workflow {

  main:
  ch_samples = channel.fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row ->
      [
        [id: row.sample, genome_size: row.genome_size],
        file(row.fastq, checkIfExists: true),
      ]
    }
  ch_files_for_qc = ch_samples
  val_checkm2_db = file(params.checkm2_db, checkIfExists: true)
  val_contaminant_fasta = [[id: 'contaminant_db'], file(params.contaminant, checkIfExists: true)]

  nanoqFilter(ch_samples)
  ch_files_for_qc = ch_files_for_qc.mix(nanoqFilter.out.fastq)

  /*
    Contamination Screening
   */

  if (!params.skip_decontamination) {
    makeBlastDb(val_contaminant_fasta)

    fastqToFasta(nanoqFilter.out.fastq)
    blastn(fastqToFasta.out.fasta, makeBlastDb.out.db_dir)
    getContaminatedIDs(blastn.out.blast_report)

    ch_samples_with_contaminated_ids = nanoqFilter.out.fastq
      .join(getContaminatedIDs.out.contaminated_ids)

    removeContaminatedIDs(ch_samples_with_contaminated_ids)

    ch_samples_clean = removeContaminatedIDs.out.fastq
    ch_files_for_qc = ch_files_for_qc.mix(removeContaminatedIDs.out.fastq)
  }
  else {
    ch_samples_clean = nanoqFilter.out.fastq
  }

  /* 
    Downsampling
   */
  rasusa(ch_samples_clean)

  /*
    Assembly and Polishing
   */
  flye(rasusa.out.fastq.map { meta, fastq -> [meta.findAll { k, _v -> k != 'step' }, fastq] })

  ch_assemblies_with_reads = flye.out.assembly
    .map { meta, assembly -> [meta, assembly] }
    .join(
      nanoqFilter.out.fastq.map { meta, fastq -> [meta.findAll { k, _v -> k != 'step' }, fastq] }
    )

  alignToAssembly(ch_assemblies_with_reads)
  postAlignToAssembly(
    alignToAssembly.out.bam
      .join(ch_assemblies_with_reads.map { meta, _assembly, fastq -> [meta, fastq] })
  )

  if (!params.skip_polishing) {
    doradoPolish(
      ch_assemblies_with_reads
        .map { meta, assembly, _fastq -> [meta, assembly] }
        .join(postAlignToAssembly.out.bam)
    )
    compareAssemblies(flye.out.assembly.join(doradoPolish.out.polished_assembly))

    flye.out.assembly | view
    doradoPolish.out.polished_assembly | view

    ch_polished_assemblies = doradoPolish.out.polished_assembly
    ch_assembly_comparisons = compareAssemblies.out.report
    ch_assemblies_collect = flye.out.assembly
      .mix(doradoPolish.out.polished_assembly)
      .collect { _meta, assembly -> assembly }
  }
  else {
    ch_polished_assemblies = channel.empty()
    ch_assembly_comparisons = channel.empty()
    ch_assemblies_collect = flye.out.assembly.collect { _meta, assembly -> assembly }
  }

  /*
    QC Reports
   */
  nanoqReport(ch_files_for_qc)
  fastQC(ch_files_for_qc)

  quast(ch_assemblies_collect)

  checkm2(ch_assemblies_collect, val_checkm2_db)

  assemblyCoverage(postAlignToAssembly.out.bam)

  ch_multiqc_reports = nanoqReport.out.report
    .mix(fastQC.out.report)
    .mix(assemblyCoverage.out.report)
    .map { _meta, report -> report }
    .mix(quast.out.report)
    .mix(checkm2.out.report)

  multiQC(ch_multiqc_reports.collect())

  publish:
  filtered_reads       = nanoqFilter.out.fastq
  subsampled_reads     = rasusa.out.fastq
  draft_assemblies     = flye.out.assembly
  polished_assemblies  = ch_polished_assemblies
  assembly_comparisons = ch_assembly_comparisons
  nanoq_reports        = nanoqReport.out.report
  fastqc_reports       = fastQC.out.report
  quast_results        = quast.out.quast_dir
  checkm2_results      = checkm2.out.checkm2_dir
  multiqc_report       = multiQC.out.report
}

output {
  filtered_reads {
    path { 'reads/postqc' }
  }
  subsampled_reads {
    path { 'reads/subsampled' }
  }
  draft_assemblies {
    path { 'assemblies/draft' }
  }
  polished_assemblies {
    path { 'assemblies/polished' }
  }
  assembly_comparisons {
    path { 'assemblies/comparison' }
  }
  nanoq_reports {
    path { 'reports/nanoq' }
  }
  fastqc_reports {
    path { 'reports/fastqc' }
  }
  quast_results {
    path { 'reports/quast' }
  }
  checkm2_results {
    path { 'reports/checkm2' }
  }
  multiqc_report {
    path { 'reports/multiqc' }
  }
}


process nanoqReport {
  label 'nanoq'
  tag "${name}"

  input:
  tuple val(meta), path(fastq)

  output:
  tuple val(meta), path("${name}.txt"), emit: report

  script:
  name = "${meta.id}" + (meta.containsKey('step') ? "_${meta.step}" : '')
  """
  nanoq \\
    --stats \\
    --input ${fastq} \\
    --report ${name}.txt \\
    -vv
  """
}

process nanoqFilter {
  label 'nanoq'
  tag "${meta.id}"

  input:
  tuple val(meta), path(fastq)

  output:
  tuple val(meta), path("${meta.id}.postqc.fastq.gz"), emit: fastq

  script:
  meta = meta + [step: 'postqc']
  """
  nanoq \\
    --min-len ${params.nanoq_min_length} \\
    --min-qual ${params.nanoq_min_quality} \\
    --input ${fastq} \\
    --output ${meta.id}.postqc.fastq.gz \\
    -vv
  """
}

process fastQC {
  label 'fastqc'
  tag "${name}"

  memory 4.GB

  input:
  tuple val(meta), path(fastq)

  output:
  tuple val(meta), path('*.{html,zip}'), emit: report

  script:
  name = "${meta.id}" + (meta.containsKey('step') ? "_${meta.step}" : '')
  """
  fastqc \\
    --memory ${task.memory.toMega()} \\
    --outdir . \\
    ${fastq}
  """
}

process multiQC {
  label 'multiqc'

  input:
  path reports, stageAs: 'input_reports/*/*'

  output:
  path "multiqc_report.html", emit: report

  script:
  """
  multiqc input_reports
  """
}

process rasusa {
  label 'rasusa'
  tag "${meta.id}"

  input:
  tuple val(meta), path(fastq)

  output:
  tuple val(meta), path("${meta.id}_subsampled.fastq.gz"), emit: fastq

  script:
  meta = meta + [step: 'subsampled']
  coverage_bases = params.rasusa_coverage && meta.containsKey('genome_size')
    ? (params.rasusa_coverage.toInteger() * meta.genome_size.toInteger()).toString()
    : params.rasusa_default_bases
  """
  rasusa reads \\
    --output ${meta.id}_subsampled.fastq.gz \\
    --bases ${coverage_bases} \\
    --seed 12345 \\
    ${fastq}
  """
}

process flye {
  label 'flye'
  tag "${meta.id}"

  cpus 8

  input:
  tuple val(meta), path(fastq)

  output:
  tuple val(meta), path("${meta.id}.fasta"), emit: assembly
  tuple val(meta), path("${meta.id}_assembly/"), emit: assembly_dir

  script:
  genome_size_opt = meta.containsKey('genome_size') ? "--genome-size ${meta.genome_size}" : ''
  """
  flye \\
    --nano-hq ${fastq} \\
    ${genome_size_opt} \\
    --out-dir ${meta.id}_assembly/ \\
    --threads ${task.cpus}

  cp ${meta.id}_assembly/assembly.fasta ${meta.id}.fasta
  """
}

process alignToAssembly {
  label 'dorado'
  tag "${meta.id}"

  cpus 4

  input:
  tuple val(meta), path(assembly), path(fastq)

  output:
  tuple val(meta), path("${meta.id}.bam"), emit: bam

  script:
  """
  dorado aligner --threads ${task.cpus} ${assembly} ${fastq} > ${meta.id}.bam
  """
}

process postAlignToAssembly {
  label 'samtools'
  tag "${meta.id}"

  cpus 2

  input:
  tuple val(meta), path(bam, stageAs: 'input/*'), path(fastq)

  output:
  tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bam

  script:
  """
  basecalling_model=\$(zcat ${fastq} | head -1 | sed -n 's/.*\\(\\(dna_r[0-9.]*_e[0-9.]*_[0-9]*bps\\|rna[0-9]*_[0-9]*bps\\)_[a-z]*@v[0-9.]*\\).*/\\1/p')

  # Fixes BAM tag needed for Dorado polishing and sorts the BAM
  samtools addreplacerg \
    -r "@RG\tID:A\tDS:basecall_model=\${basecalling_model}" \
    ${bam} \
  | samtools sort -@ ${task.cpus} -o ${meta.id}.bam

  samtools index ${meta.id}.bam
  """
}

process doradoPolish {
  label 'dorado'
  label 'gpu'
  tag "${meta.id}"

  cpus 2
  memory 16.GB

  input:
  tuple val(meta), path(assembly), path(bam), path(bai)

  output:
  tuple val(meta), path("${meta.id}_polished.fasta"), emit: polished_assembly

  script:
  """
  dorado polish \\
    --bacteria \\
    ${meta.id}.bam \\
    ${assembly} \\
    > ${meta.id}_polished.fasta
  """
}

process quast {
  label 'quast'

  input:
  path assemblies

  output:
  path 'quast_results/report.tsv', emit: report
  path 'quast_results/', emit: quast_dir

  script:
  """
  quast \\
    -o quast_results/ \\
    ${assemblies}
  """
}

process checkm2 {
  label 'checkm2'

  cpus 8

  input:
  path assemblies
  path database

  output:
  path 'checkm2_results/quality_report.tsv', emit: report
  path 'checkm2_results/', emit: checkm2_dir

  script:
  """
  checkm2 predict \\
    --input ${assemblies} \\
    --output_directory checkm2_results/ \\
    --database_path ${database} \\
    --threads ${task.cpus}
  """
}

process makeBlastDb {
  label 'blast'

  input:
  tuple val(meta), path(fasta)

  output:
  tuple val(meta), path('blast_db/'), emit: db_dir

  script:
  """
  mkdir -p blast_db/
  makeblastdb -in ${fasta} -dbtype nucl -out blast_db/${meta.id}
  """
}

process fastqToFasta {
  label 'seqkit'
  tag "${meta.id}"

  input:
  tuple val(meta), path(fastq)

  output:
  tuple val(meta), path("${meta.id}.fasta"), emit: fasta

  script:
  """
  seqkit fq2fa \\
    -o ${meta.id}.fasta \\
    ${fastq}
  """
}

process blastn {
  label 'blast'
  tag "${meta.id}"

  cpus 4

  input:
  tuple val(meta), path(fastq)
  tuple val(meta2), path(db_dir)

  output:
  tuple val(meta), path("${meta.id}_contamination.csv"), emit: blast_report

  script:
  """
  blastn \\
    -query ${fastq} \\
    -db ${db_dir}/${meta2.id} \\
    -out ${meta.id}_contamination.pre.csv \\
    -outfmt '20 qseqid sseqid pident qcovs length evalue bitscore' \\
    -subject_besthit \\
    -max_target_seqs 1 \\
    -num_threads ${task.cpus}

  # Remove duplicated headers (blastn + outfmt 20 bug?)
  cat ${meta.id}_contamination.pre.csv | awk 'NR==1 || !/^qseqid,/' > ${meta.id}_contamination.csv
  """
}

process getContaminatedIDs {
  label 'qsv'
  tag "${meta.id}"

  input:
  tuple val(meta), path(blast_report)

  output:
  tuple val(meta), path("${meta.id}_contaminated_ids.txt"), emit: contaminated_ids

  script:
  """
  qsv luau filter \\
    'tonumber(pident) >= ${params.blast_min_pident} and tonumber(qcovs) >= ${params.blast_min_qcovs}' \\
    ${blast_report} \\
  | qsv select qseqid \\
  | qsv behead \\
  | uniq \\
  > ${meta.id}_contaminated_ids.txt
  """
}

process removeContaminatedIDs {
  label 'seqkit'
  tag "${meta.id}"

  input:
  tuple val(meta), path(fastq), path(contaminated_ids)

  output:
  tuple val(meta), path("${meta.id}_decontaminated.fastq.gz"), emit: fastq

  script:
  meta = meta + [step: 'decontaminated']
  """
  seqkit grep \\
    -v -f ${contaminated_ids} \\
    -o ${meta.id}_decontaminated.fastq.gz \\
    ${fastq}
  """
}

process assemblyCoverage {
  label 'samtools'
  tag "${meta.id}"

  input:
  tuple val(meta), path(bam), path(bai)

  output:
  tuple val(meta), path("${meta.id}.txt"), emit: report

  script:
  """
  samtools coverage ${bam} > ${meta.id}.txt
  """
}

process compareAssemblies {
  label 'mummer'
  tag "${meta.id}"

  input:
  tuple val(meta), path(draft_assembly), path(polished_assembly)

  output:
  tuple val(meta), path("${meta.id}_dnadiff.report"), emit: report

  script:
  """
  dnadiff ${draft_assembly} ${polished_assembly} -p ${meta.id}_dnadiff
  """
}
