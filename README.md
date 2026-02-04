# Exon target space paper: code

## Data directory (`MYDATA`)

Some scripts use:

```
data_dir <- Sys.getenv("MYDATA")
```

Data from Zenodo should be downloaded into a directory and the `MYDATA` 
environment variable needs to be defined

```
export MYDATA=<path to zenodo downloads>
```
