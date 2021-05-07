# IEnKS
A C++ implementation of the IEnKS-Q and variants using Armadillo for efficient
linear algebra. Python controls the executable and post-process its results for
ease of use.

Lorenz95 and Lorenz63 dynamical models are provided in ```C/lorenz95.h```.

## Compilation

```
cd Debug
make
```

This compiles the C++ and produces the ```IEnKS``` executable. It relies on the
libraries Armadillo and cnpy. These library are provided for convenience.

## Usage

The IEnKS executable is controlled with a Python3 script.
Such a script is provided in ```Python/test.py``` as an example.

This script takes a list of list of parameters, launch the Data Assimilation
algorithms in parallel over each possible parameters combination and plot scores
as function of the parameters.



