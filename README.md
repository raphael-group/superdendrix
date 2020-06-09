# SuperDendrix

SuperDendrix is an algorithm that uses an integer linear program (ILP) for identifying mutually exclusive sets of genomic features that are correlated with a dependency profile.

SuperDendrix consists of three modules:
1) Scoring differential dependencies and selecting genomic and cell-type features.
2) Finding feature sets associated with differential dependencies. 
3) Evaluating statistical significance of associations. 

### Set up

#### Python
SuperDendrix modules are written in R and Python 3, and has some Python module dependencies. We suggest using Anaconda to manage Python dependencies. The dependencies are listed in the `environment.yml` file in this repository.

To solve the ILP, SuperDendrix uses the [Gurobi Optimizer](http://www.gurobi.com/downloads/gurobi-optimizer), accessed through the `gurobi` Python module. Gurobi must be installed in order to run SuperDendrix.

### Usage

To run SuperDendrix on on the CERES dataset, use the following commands:

    python superdendrix.py -m ../data/examples/kim2016/events/revealer-ex1-ctnnb1-events.tsv -T  ../data/examples/kim2016/outcomes/revealer-ex1-ctnnb1-reporter.tsv -Tc Z-score -k 3
    
