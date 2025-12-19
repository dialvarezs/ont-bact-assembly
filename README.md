# ONT Bacterial Assembly

Small Nextflow pipeline that performs:
- read qc (nanoq, fastqc)
- decontamination (blast)
- downsampling (rasusa)
- assembly (flye)
- polishing (dorado)
- assembly qc (quast, checkm2)
- comparison between draft and polished (mummer dnadiff)
- assembly coverage (samtools)
- report generation (multiqc)

## Usage

Requires Nextflow and Singularity/Apptainer.

```bash
nextflow run main.nf \
	--samplesheet input/samplesheet.csv \
	--contaminant input/contaminant/contaminant.fasta \
	--checkm2_db input/databases/checkm2_db.dmnd \
```

Check the `params` section of `nextflow.config` for additional parameters that can be set.

## Input files

### Samplesheet

```csv
sample,fastq,genome_size
ONT-S1,input/fastq/ONT-S1.fastq.gz,6000000
ONT-S2,input/fastq/ONT-S2.fastq.gz,5000000
```

### Contaminant file

A fasta file containing contaminant sequences to be removed from the input fastq files.

### CheckM2 database

Download it from https://zenodo.org/records/5571251 and decompress it, then provide the path of the `.dmnd` file to the `--checkm2_db` parameter.