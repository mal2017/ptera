# Pipeline for Transposon Expression Reanalysis (PTERA)

## Dependencies

Using mamba 0.9.2 (or conda 4.8.4) install the following:

```bash
mamba create -n ptera_v1 \
	python=3.9.6 snakemake=6.7.0 peppy=0.31.1 sra-tools=2.11.0 tabix=1.11 \
	r-tidyverse=1.3.1 r-tidymodels=0.1.3 r-sparklyr=1.7.1 r-sqldf=0.4_11 \
	xsv=0.13.0 bioconductor-deseq2=1.32.0 bioconductor-preprocesscore=1.54.0 \
	bioconductor-edaseq=2.26.0 r-wgcna=1.69 samtools=1.14 bioconductor-macsr=1.2.0
```

Make sure you run pre-install the macs env 1 time via basilisk.

```
basilisk::basiliskStart(MACSr:::env_macs)
```

Singularity must also be installed and accessible in this environment. Depending on
your system this can be done via conda. For HPCs admin privileges may be required.

Other dependencies are handled at runtime by singularity.

## Where the data comes from

Accessions are included in the pepfiles for each subworkflow. The primary data (DGRP RNA-seq) was
generated by downloading all metadata for [PRJNA483441](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA483441&o=acc_s%3Aa#)
and processing the metadata file with a custom R script (see workflow/scripts).

Other data were similarly pulled from SRA/GEO.

## Prerequisite data resources

This pipeline requires a `resources` folder which should be downloaded and placed in the
pipeline directory (level with `config/` for example).

## Setting up new subworkflows

New data can be processed in subworkflows. Each sample table should resemble
the SRA run selector export. Recommended to clean this table appropriately with a script
that can be saved in the `workflow/scripts` directory of each subworkflow.

Each table should have at a minimum the following columns:

- `sample_name`
- `LibraryLayout`
- `Experiment`
- `Run`
- `Library`

ChIP-seq data should additionally have an `input` column.

## Test Mode

Note that by default `RUN_TYPE=TEST` for the resource intensive workflows. This is useful for testing individual pipeline
components (i.e. subworkflows or the initial portions of the main workflow), but due to
the difficulty of making sure equivalent strain samples are always processed in the
test dataset for each subworkflow, test mode is not guaranteed to complete the full workflow
(main wf + subworkflows). Instead use `RUN_TYPE=FULL`.

## Running the pipeline

```bash
snakemake --profile <your profile> --use-singularity -j 999 -kp \
	--config RUN_TYPE=FULL -n
```

If disk space is at a premium, consider running with `--prioritize salmon_quant_se_vanilla `
to enforce creation of the quant output and deletion of temp files for each sample as fast as possible.

Additionally, batching may be used by specifying `--batch vanilla_salmon_terminus_collapse=1/N`
where `N` is the total number of batches the collecting rule is divided into.
This can then be repeated (`2/N` and so on) until all batches are run.
