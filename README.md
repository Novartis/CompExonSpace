# Specificity and exon target space of splicing modifying compounds — analysis code

This repository contains the code used to analyze the data presented in the manuscript **"Specificity and exon target space of splicing modifying compounds"**.

## Data availability

The analysis requires downloading the associated dataset from Zenodo:

- DOI: `10.5281/zenodo.18416224`

Download and extract the Zenodo files into a local directory of your choice.

## Configure the data path (required)

Some scripts read the data location from the `MYDATA` environment variable (e.g., `data_dir <- Sys.getenv("MYDATA")`).

Before running the analysis, set `MYDATA` to the path where you downloaded/extracted the Zenodo dataset.

### macOS / Linux

```bash
export MYDATA=<path to zenodo downloads>
```

### Windows (PowerShell)
setx MYDATA "<path to zenodo downloads>"

After setting it, start a new terminal/session (or restart R/RStudio) so the environment variable is available.

Quick start
Clone this repository.
Download and extract the Zenodo dataset (DOI above).
Set MYDATA to the extracted dataset directory.
Run the analysis scripts (see repository scripts for entry points).


