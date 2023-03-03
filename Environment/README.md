This directory contains a `.yml` file that can be used for building an environment with `conda` with all the necessary packages for running the `microSLURM_16S`. To create the environment (titled `microSLURM_16S`) use the following command:
```
$ conda env create -f microSLURM_16S.yml
```
and once the environment is built, activate it with
```
$ conda activate microSLURM_16S
```
If this environment is built, then `conda activate microSLURM_16S` should suffice for the `-p` parameter in the `microSLURM_16S.sh` script, but you will need to install the R packages `dada2`, `DECIPHER`, `phangorn`, and `phyloseq` from Bioconductor[https://www.bioconductor.org/] within the `microSLURM_16S` environment before running.