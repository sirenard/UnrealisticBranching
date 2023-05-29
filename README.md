# Unrealistic Branching

This code implements a new branching strategies into the [SCIP](https://scipopt.org/) solver. 
The purpose of this branching is to have the best branching quality without regard of its complexity.

It also implements a way to generate a dataset from solving an instance. 
This dataset can be used to learn a new strategy by training a dataset. 
This model can be used itself as a branching strategy.

The dataset generation and the simple usage of the unrealistic branching strategy supports parallelization using MPI.

## Usage

New commands added to SCIP:

* `generateDataset <problemFile> <outputFile> <signA> <signB> <signC>`: generate a dataset in the csv file `outputFile` by solving the problem in the file `problemFile`. `signA` is 0 if the matrix A is composed of both negative and positive values, 1 if it is composed only of positive values, -1 if is only composed of negative values. Respectively the same value for `signB` and `signC` for the matrixes B and C.
* `trainmodel <modelName> <datasetPath> <nTrees>`: train a random forest of `nTrees` from the `datasetPath` csv file. The model is saved in the `modelname` file.
* `loadmodel <modelName> <signA> <signB> <signC>`: load a model, it will be used when the branching strategy `unrealisticTrained` is used.

To solve an instance using a trained model. The problem must be loaded before loading the model. So:

```
read problem.lp
loadmodel <modelname> ...
optimize
```

### Unrealistic Branching

The Unrealistic Branching is added to scip, it can be parametrized using the following parameters:

* `branching/unrealistic/recursiondepth`: How depth the unrealistic branching is performed with the recursion, -1 for unlimited depth. Default: 1
* `branching/unrealistic/leaftimelimit`: Time limit allowed to the leaf of the recursion tree (not solved using the unrealistic branching). -1 for no limit.
* `branching/unrealisticTrained/threshold`: Above which threshold does the predicted score is relevant (for the 1v1 approach)

### Dataset Generation Configuration

* `dataset/scoremethod`: How the score must be computed, can take the following values: 
    - 'c'ontinuous mapping between 0 and 1. 
    - 'a'lpha: the score is 1 if variable SB score is <= (1+alpha)bestUbScore
    - '1'vs1: apply the 1vs1 method. The feature vector is the descriptio of 2 candidates and the output value is 1 is the first is better than the other 
    - 'n'one: no score computed for each features vectorthe current SB score, the minimum and maximum SB score are reported
* `dataset/alpha`: If the scoremethod is 'a', threshold parameter in [0,1]. By default: 0.2
*  `dataset/epsilon`: Exploration branching is performed with probability epsilon at each step during dataset generation


## Automatic Testing

A python script [test.py](/test.py) allows to perform automatic tests.

```
usage: test.py [-h] -e EXECUTABLE -d DIR [-m MODELS [MODELS ...]] [-o OUTPUT]
               [-f FROMS] [-t TO] [-r RECURSIONDEPTH [RECURSIONDEPTH ...]] -a
               ABC ABC ABC [-md MAXDEPTH [MAXDEPTH ...]]
               [-th THRESHOLDS [THRESHOLDS ...]]
               [-bs BRANCHINGSTRATEGIES [BRANCHINGSTRATEGIES ...]]
               [-eps EPSILON] [-exp EXPLOREMETHOD] [-tl TIMELIMIT]
               [-ltl LEAFTIMELIMIT] [-nl NODELIMIT] [-of OUTPUTFILE]

optional arguments:
  -h, --help            show this help message and exit
  -e EXECUTABLE, --executable EXECUTABLE
                        path to scip learning executable
  -d DIR, --dir DIR     path to lp instances
  -m MODELS [MODELS ...], --models MODELS [MODELS ...]
                        path to models to use
  -o OUTPUT, --output OUTPUT
                        path to output directory
  -f FROMS, --froms FROMS
                        index from
  -t TO, --to TO        index to
  -r RECURSIONDEPTH [RECURSIONDEPTH ...], --recursiondepth RECURSIONDEPTH [RECURSIONDEPTH ...]
                        list of recursion (in addition to 0)
  -a ABC ABC ABC, --abc ABC ABC ABC
                        list a b c
  -md MAXDEPTH [MAXDEPTH ...], --maxdepth MAXDEPTH [MAXDEPTH ...]
                        list of max depths tested
  -th THRESHOLDS [THRESHOLDS ...], --thresholds THRESHOLDS [THRESHOLDS ...]
                        threshold values
  -bs BRANCHINGSTRATEGIES [BRANCHINGSTRATEGIES ...], --branchingstrategies BRANCHINGSTRATEGIES [BRANCHINGSTRATEGIES ...]
                        Branching strategies used
  -eps EPSILON, --epsilon EPSILON
                        Epsilon for exploration
  -exp EXPLOREMETHOD, --exploremethod EXPLOREMETHOD
                        Name of the branching scheme used for exploration
  -tl TIMELIMIT, --timelimit TIMELIMIT
                        if solve take more than that time, don't do it
  -ltl LEAFTIMELIMIT, --leaftimelimit LEAFTIMELIMIT
                        leaftime limit (coefficient of the running time)
  -nl NODELIMIT, --nodelimit NODELIMIT
                        if solve take more than in node number that time,
                        don't do it
  -of OUTPUTFILE, --outputfile OUTPUTFILE
```