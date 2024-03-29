# SuperDendrix

<img src="fig/overview.png" width="800">

SuperDendrix is an algorithm that uses an integer linear program (ILP) for identifying mutually exclusive sets of genomic features that are correlated with a dependency profile.
This repository includes instructions for installation and tutorials using example data for SuperDendrix.

*This README is work in progress.*

SuperDendrix consists of three modules:
1) Scoring differential dependencies and selecting genomic and cell-type features.
2) Finding feature sets associated with differential dependencies. 
3) Evaluating statistical significance of associations. 

## Set up

### Python and R
SuperDendrix modules are written in R and Python 3, and has some R and Python module dependencies. We suggest using Anaconda to manage the dependencies. The dependencies are listed in the `environment.yml` file in this repository.

In addition, code from the following repo needs to be downloaded under `superdendrix/` directory.
- oncokb-annotator: https://github.com/oncokb/oncokb-annotator

To solve the ILP, SuperDendrix uses the [Gurobi Optimizer](http://www.gurobi.com/downloads/gurobi-optimizer), accessed through the `gurobi` Python module. Gurobi must be installed in order to run SuperDendrix.


## Usage

### Required data

SuperDendrix requires the following data:

1. Quantitative phenotypes (e.g. gene dependency) across samples.
2. A list of features (e.g. genomic alterations) for each sample.
3. Categorical information (e.g. cell type) of each sample.

### Downloading required data
We provide an example dataset for testing SuperDendrix which can be downloaded using the following command.

    snakemake all

The example dataset includes a dependency scores file (`dependency_scores.csv`) containing BRAF, NRAS, and KRAS dependency profiles, a mutation file (`mutations.tsv`), and a sample information file (`sample_info.tsv`) which are part of the DepMap dataset [20Q2 release] used in our analysis.

## Commands
The following sections describes the three modules of SuperDendrix and provides instructions for testing them on the provided example dataset.

### Module 1
Module 1 of SuperDendrix scores dependencies from results of gene perturbation experiments and constructs mutation features (and optionally cancer-type features) using the OncoKB database of reported cancer mutations.

We provide a bashscript for running module 1 on the example dataset which can be executed using the following command.

    bash module_1.sh

### Modules 2 and 3
Module 2 identifies a set of approximately mutually exclusive mutation features (and optionally cancer-type features) that are associated with a differential dependency profile.

Module 3 contains two steps: first is a model selection to select the features that contribute significantly to the association identified from the second module. The second step evaluates the statistical significance of the association between the differential dependency and the selected set of features.

Modules 2 and 3 for the BRAF differential dependency profile from the example dataset can be run in a single bashscript using the following command.

    bash modules_2_3.sh

## Contacts
SuperDendrix repository is maintained by Tyler Park, a current member of the research group of prof. Ben Raphael at Princeton University.
For any related question, please email Tyler Park (typark@princeton.edu).

### License
See `LICENSE` for license information.

### Citation
The SuperDendrix manuscript will be available soon.