
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Apjok-et-al-DEEPMINE-NatMicrobiol

Code deposited for paper “Characterization of antibiotic resistomes by
reprogrammed bacteriophage-enabled functional metagenomics in clinical
strains”, Apjok et al. Nature Microbiology.

## System requirements

- [Apptainer/Singularity](https://apptainer.org/) (for Illumina
  demultiplexing, only works on linux)
- [R](https://www.r-project.org/) (for all other scripts)

## Download data

Illumina reads and Nanopore contigs can be downloaded from the European
Nucleotide Archive. Study Accession: PRJEB54063 (Secondary Accession:
ERP138883)

``` r
suppressMessages(source("prepare_input_data.R"))
```

## Demultiplex Illumina reads

Download Apptainer/Singularity container for demultiplexing:

``` r
# Define a local directory that will store the apptainer/singularity container
cache_dir <- "./containers"

# convert path to absolute path
cache_dir <- normalizePath(cache_dir)

# create cache directory if it doesn't already exists
if (!dir.exists(cache_dir)) {
  dir.create(cache_dir, recursive = TRUE)
}

# download container image if it is not already downloaded
if (!file.exists(paste0(cache_dir, "/mesti90-hgttree-latest.img"))) {
  cmd <- paste0(
    "singularity pull --dir ", cache_dir,
    " mesti90-hgttree-latest.img docker://mesti90/hgttree:latest"
  )
  system(cmd)
}
```

Demultiplex Illumina reads in first run:

``` r
setwd("./run1")
cmd <- paste0(
  "singularity exec ", cache_dir,"/mesti90-hgttree-latest.img ",
  "python3 barcode_amplicon_demux_20201026.py3"
)
system(cmd, ignore.stdout = TRUE)
```

Demultiplex illumina reads in second run:

``` r
setwd("./run2")
cmd <- paste0(
  "singularity exec ", cache_dir,"/mesti90-hgttree-latest.img ",
  "python3 barcode_amplicon_demux_20210201.py3"
)
system(cmd, ignore.stdout = TRUE)
```

Demultiplex illumina reads in third run:

``` r
setwd("./run3")
cmd <- paste0(
  "singularity exec ", cache_dir,"/mesti90-hgttree-latest.img ",
  "python3 barcode_amplicon_demux_20220505.py3"
)
system(cmd, ignore.stdout = TRUE)
```

## Prepare supplementary table 6

``` r
setwd("./supplementary_table_6")
suppressMessages(source("Script - Supplementary table 6.R"))
```

## Prepare supplementary table 7

``` r
setwd("./supplementary_table_7")
suppressMessages(source("Script - Supplementary table 7.R"))
```

## Prepare supplementary table 9

``` r
setwd("./supplementary_table_9")
suppressMessages(source("Script - Supplementary table 9.R"))
```
