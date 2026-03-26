# CompExonSpace  
Paper: Specificity and exon target space of splicing modifying compounds

This repository contains all R and Python code used to perform the analyses and plot figures published in the manuscript 'Specificity and exon target space of splicing modifying compounds'.

---

## Structure

```
CompExonSpace/
│
├── PCA(FigS8).Rmd
├── PCA_GA(FigS9).Rmd
├── CAGA_violin(Fig3F).Rmd
├── ESE_ESS_PSIcorrelation.ipynb
├── comparisonROC.ipynb
├── SHAPanalysis.ipynb
├── ... (other R and Python analysis scripts)
│
├── renv/                      # R environment managed by renv
├── renv.lock                  # R package versions
├── environment.yml            # Python environment
├── README.md
└── LICENSE
```

---

# Installation

## 1. R environment (renv)

Make sure you are inside the project directory, then run:

```r
install.packages("renv")
renv::restore()
```

This recreates the exact R package versions used for all `.Rmd` analyses.

---

## 2. Python environment (conda)

```bash
conda env create -f environment.yml
conda activate comp-exon-space
```

This installs all dependencies needed to run the Python notebooks.

---

# Data Availability

## 1. Raw RNA‑seq Data (SRA)

The RNA‑seq datasets used in this study are available at:

NCBI SRA BioProject: PRJNA1293306

The BioProject contains all SRR run accessions produced for this study.  
To obtain a list of all SRR runs do this:

```bash
esearch -db sra -query PRJNA1293306 | efetch -format runinfo > runinfo.csv
```

Raw FASTQ files are not needed to reproduce the figures in this repository, but can be downloaded like this:

```bash
prefetch SRRXXXXXXX
fasterq-dump SRRXXXXXXX
```

For full reproducibility of the analyses performed here, processed data from Zenodo are sufficient.

---

## 2. Processed Data (Zenodo)

Processed datasets necessary to reproduce every figure (e.g., PSI matrices, exon-level summaries, ML input tables) are available on Zenodo:

Zenodo record: https://zenodo.org/records/18416225

Download and extract:

```bash
wget https://zenodo.org/records/18416225/files/data.zip
unzip data.zip
```

The data archive contains:

```
data/
├── exon_diff/
├── psi_persample/
├── psi_pergasplicesite/
├── ML/
└── ... additional processed files
```

---

# Set the Data Path (`MYDATA`)

The scripts refer to the processed data through the environment variable `MYDATA`.

To set this do:

```r
Sys.setenv(MYDATA = "/path/to/data")
```

After the download the directory should look like:

```
$MYDATA/
   exon_diff/
   psi_persample/
   psi_pergasplicesite/
   ML/
   ...
```

Scripts such as:

```r
run_ga_pca(exp_id, data_dir = Sys.getenv("MYDATA"))
```

will automatically load data from this location.

---

## Python notebooks

Launch Jupyter inside the conda environment:

```bash
conda activate comp-exon-space
jupyter notebook
```

Open any notebook:

- `SHAPanalysis.ipynb`
- `comparisonROC.ipynb`
- `ESE_ESS_PSIcorrelation.ipynb`

---

# License

This project is licensed under the terms of the LICENSE file in the repository.
