# SuperDendrix

![alt text](https://www.dropbox.com/s/fno4102uky1rhj0/fig1.pdf?dl=0)


SuperDendrix is an algorithm that uses an integer linear program (ILP) for identifying mutually exclusive sets of genomic features that are correlated with a dependency profile.
This repository includes instructions for installation and tutorials using example data for SuperDendrix.
This README is work in progress.

SuperDendrix consists of three modules:
1) Scoring differential dependencies and selecting genomic and cell-type features.
2) Finding feature sets associated with differential dependencies. 
3) Evaluating statistical significance of associations. 

## Set up

### Python and R
SuperDendrix modules are written in R and Python 3, and has some R and Python module dependencies. We suggest using Anaconda to manage the dependencies. The dependencies are listed in the `environment.yml` file in this repository.

In addition, the following modules need to be installed:
- EMMIXskew 1.0.3 (R package): https://cran.r-project.org/src/contrib/Archive/EMMIXskew/
- oncokb-annotator: https://github.com/oncokb/oncokb-annotator

To solve the ILP, SuperDendrix uses the [Gurobi Optimizer](http://www.gurobi.com/downloads/gurobi-optimizer), accessed through the `gurobi` Python module. Gurobi must be installed in order to run SuperDendrix.


## Usage

### Required data

SuperDendrix requires the following data:

1. Quantitative phenotypes (e.g. gene dependency) across samples.
2. A list of features (e.g. genomic alterations) for each sample.
3. (Optional) Categorical information (e.g. cell type) of each sample.

### Downloading required data

An example dataset for testing SuperDendrix can be downloaded using the following command.

    snakemake all

## Commands

SuperDendrix modules are run using the following commands:

### Module 1

Compute CERES z-scores and identify the six-sigma genes.

    Rscript utils/compute_CERES_zscores.R

Fitting CERES dataset with mixtures of t-distributions to find differental dependencies.

    Rscript src/fit_tmm.R

Fitting CERES dataset with mixtures of Gaussian distributions to score differential dependencies.

    python utils/fit_gmm.py -pf ${PHENOTYPE} -o ${OUTPUT_FILE}

Annotating mutations with OncoKB database of cancer mutations.

    Rscript src/oncokb_maf_annotator.R

### Modules 2 and 3

Generating randomized feature matrices using the curveball method.

    python utils/generate_null_matrices.py -m ${FEATURES} -p ${CYCLE} -o ${OUTPUT_FILE} -pre ${PREFIX} -rs ${RANDSEED}

Identifying an association between differential dependency and a set of genomic features and conducting model selection and evaluation of statistical significance.

    python src/superdendrix.py -t ${THREADS} -T ${PHENOTYPE} -Tc ${GENE} -m ${FEATURES} -gf ${FEATURELIST} -p ${CYCLE} -cp ${CP} -d ${DIRECTION} -k ${SETSIZE} -nm ${NULLMATRICES} -rs ${RANDSEED} -x -curve -o ${OUTPUT_FILE}

## Demo
A bashscript for an example analysis of BRAF dependency profile from the 20Q2 release of DepMap dataset is provided in the demo directory.

