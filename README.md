# Tensor Belief Propagation 

[![Build Status](https://travis-ci.org/akxlr/tbp.svg?branch=master)](https://travis-ci.org/akxlr/tbp)

[Tensor Belief Propagation](http://proceedings.mlr.press/v70/wrigley17a/wrigley17a.pdf) (TBP) is an experimental algorithm for
approximate inference in discrete graphical models [1]. It takes a factor graph in [.uai](#other-file-formats) or [.fg](#other-file-formats) format and outputs approximate marginals for
each variable.

[1] [Wrigley, Andrew, Wee Sun Lee, and Nan Ye. "Tensor Belief Propagation." International Conference on Machine Learning. 2017.](http://proceedings.mlr.press/v70/wrigley17a/wrigley17a.pdf)

## Requirements

 * Linux or OSX
 * Python 3.6+

## Installation

Install libDAI prerequisites:
```
# Linux
$ sudo apt-get install g++ make doxygen graphviz libboost-dev libboost-graph-dev libboost-program-options-dev libboost-test-dev libgmp-dev cimg-dev

# OSX
$ brew install boost gmp doxygen graphviz
```

Install tbp with the Python package manager `pip`:
```bash
$ pip install tbp
...
Successfully installed tbp-X.X.X
```
This will take a while as libDAI must be compiled.

## Usage
TBP takes a factor graph in either [.fg](#other-file-formats) or [.uai](#other-file-formats) format as input, and
outputs the approximate marginal distribution of each variable in [.MAR format](#other-file-formats).
This involves two steps — first, all potential functions in the graph must be decomposed
into sums of rank-1 tensors yielding a decomposed factor graph (.dfg). Then, the
message passing procedure must be run on the decomposed graph to give approximate marginals.

### Command line
After installation, the command line utility `tbp` is available to do either or both of these steps. For usage
instructions, run `tbp --help`. 

#### Examples

Decompose the factor graph `ising_8x8.fg` and find marginals:
```bash
$ tbp tests/ising_8x8.fg
64 2 0.594961 0.405039 2 ... 0.608573 0.391427
```

Decompose input potentials into 3 rank-1 components and save the resulting decomposed graph (but don't find marginals):
```bash
$ tbp tests/ising_8x8.fg -r 3 -o tests/ising_8x8.dfg --verbosity 2
Reading graph tests/ising_8x8.fg (libDAI format)...
Decomposing input graph (r=3 terms per factor)...
Successfully saved decomposed graph to tests/ising_8x8.dfg.
```

Decompose the factor graph `Promedus_11.uai` after applying some evidence, find marginals using TBP with sample size of 1000, and save the output
to `out.MAR`:
```bash
$ tbp tests/uai/MAR_prob/Promedus_11.uai -e tests/uai/MAR_prob/Promedus_11.uai.evid -k 1000 -o out.MAR --verbosity 2
Reading graph tests/uai/MAR_prob/Promedus_11.uai (UAI format)...
Applying evidence file tests/uai/MAR_prob/Promedus_11.uai.evid...
Decomposing input graph (r=4 terms per factor)...
Running TBP with sample size K=1000...
Successfully saved marginals to out.MAR.
```


### Python library
The `tbp` package can also be used directly from Python, for example:

```python
import tbp

# Load a factor graph in .uai format
g = tbp.load_uai_graph('tests/uai/MAR_prob/linkage_11.uai')

# Apply evidence (fixed variable assignments)
g.apply_evidence('tests/uai/MAR_prob/linkage_11.uai.evid')

# Decompose each factor into a weighted sum of 4 rank-1 tensors
dg = g.decompose(r=4)

# Run TBP to find marginals with sample size of 10000
mar = dg.tbp_marg(K=10000)

```

### Troubleshooting
#### Installing into a virtual environment
If `pip install` has issues with dependencies or version conflicts, you can install the necessary
 packages into a virtual environment (a project-specific folder rather than globally on your system):
```bash
$ sudo pip3 install virtualenv  # pip or pip3, depending on your system
$ virtualenv -p python3 venv    # create venv folder to store packages
$ source venv/bin/activate      # activate virtual environment
$ pip install tbp               # install tbp into venv folder
```
Now when you invoke `tbp`, the local versions will be used.

#### Building from GitHub clone
To use the `tbp` Python package from source without installation via `pip install`, libDAI must first be compiled:
 
```bash
$ git clone git@github.com:akxlr/tbp.git
$ cd tbp/libdai
$ cp Makefile.<platform> Makefile.conf  # Choose <platform> according to your platform
$ make
...
libDAI built successfully!
```
This produces a utility `libdai/utils/dfgmarg` which is symlinked from `tbp/dfgmarg` and used during inference. See [libDAI README](libdai/README) for full installation instructions.


## Using MATLAB for the decomposition
The decomposition of potential functions uses the non-negative CP decomposition algorithm in the Tensorly tensor 
library. As an alternative to TensorLy, the [MATLAB Tensor Toolbox](http://www.sandia.gov/~tgkolda/TensorToolbox) can be used
(this was what we used in [1]). To use this instead of Tensorly:

 * Install MATLAB
 * Install the
 [MATLAB Python API](https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html)
 * Install the [MATLAB Tensor Toolbox](http://www.sandia.gov/~tgkolda/TensorToolbox)
 
You can now replace `method='tensorly'` with `method='matlab'` when calling decomposition functions in [core.py](python/tbp/core.py).


## File formats

### .dfg (decomposed factor graph)
We created the `.dfg` file format based on
[libDAI's .fg file format](https://staff.fnwi.uva.nl/j.m.mooij/libDAI/doc/fileformats.html)
to represent decomposed factor graphs. A decomposed factor graph is a
factor graph with all factors represented as sums of rank-1 tensors rather than multidimensional tables.

The first line of a `.dfg` file contains the number of factors in the graph, followed by a blank line. Then, factors
are described in turn by blocks separated by a single blank line. Each factor block is structured as follows:

     1. n_terms
     2. <weights>
     3. n_variables
     4. <variable indices>
     5. <variable cardinalities>
     6. n_nonzero_1
     7. 1 0.5
     8. 3 0.1
     9. 4 0.1
    10. ...
    11. n_nonzero_2
    12. 1 0.5
    13. 3 0.1
    14. 4 0.3
    15. ...

In the header section of the factor block (lines 1-5), `n_terms` is the number of terms in the decomposition and
`<weights>`, `<variable indices>` and `<variable cardinalities>` are self-explanatory space-separated lists of length `n_terms`,
`n_variables` and `n_variables` respectively. 

The remainder of the factor block (line 6 onwards) describes
a series of `n_variables` 2D matrices that together describe the `n_terms` rank-1 tensors.
Each matrix corresponds to a single variable and has shape `(cardinality, n_terms)`, where `cardinality` is
the cardinality of the variable and `n_terms` is the number of rank-1 terms in the decomposition (constant
for all variables). Each matrix begins with the
number of nonzero values in the matrix, followed by a series of `index value` pairs describing the nonzero
entries of the matrix in column-major order. See
[libDAI's documentation](https://staff.fnwi.uva.nl/j.m.mooij/libDAI/doc/fileformats.html) for examples of how to
reshape these lists back into matrices.

 The *i*th rank-1 tensor is constructed by taking the outer product of the *i*th columns of
all matrices. The complete factor is then reconstructed by adding up these rank-1 tensors and weighting
according to `<weights>`.

### Other file formats
Other file formats used in this project are:

 * `.fg` (libDAI factor graph): https://staff.fnwi.uva.nl/j.m.mooij/libDAI/doc/fileformats.html
 * `.uai` (UAI factor graph): http://www.hlt.utdallas.edu/~vgogate/uai14-competition/modelformat.html
 * `.MAR` (marginals): http://www.hlt.utdallas.edu/~vgogate/uai14-competition/resformat.html
 * `.evid` (evidence): http://www.hlt.utdallas.edu/~vgogate/uai14-competition/evidformat.html

## To do
 * ICML experiments - finish cleaning code used for experiments (see `icml17.py` for partial code)
 * Rewrite code that loads .uai files to handle all problems (currently breaks on some)
 * Deal with Z <= 0 warning from C++ code
 * Clean up C++ code and compiler warnings
 * Add more tests

## Feedback
Bug reports, suggestions and comments are welcome. Please email [andrew@wrigley.io](mailto:andrew@wrigley.io) or use the issue tracker.

## License

See [LICENSE.txt](LICENSE.txt) (MIT).

## Acknowledgments

* [libDAI](https://staff.fnwi.uva.nl/j.m.mooij/libDAI/) (included in [libdai](libdai) folder with modifications; libDAI's junction tree implementation is used for the message passing step)
* [Eigen](http://eigen.tuxfamily.org/) (version 3.3.4 included in [libdai/vendor/include](libdai/vendor/include) folder)
* [TensorLy](https://github.com/tensorly/tensorly) (used to perform initial non-negative CP decomposition of potential functions)
* [MATLAB Tensor Toolbox](http://www.sandia.gov/~tgkolda/TensorToolbox)





