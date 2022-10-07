# GridMatch
A tool for rubber-sheet distorsion removal in python. These routines can be used to:
1) remove seeing distorsions from a burst of images where the object does not evolve
2) remove redisual seeing distorsions of a time series, where the object is eveolving.
3) to morph one image into another assuming that they are relatively similar.

These tools are a python port of the tools originally developed in ANA/C by R. A. Shine (LMSAL).

Coded in Python by J. de la Cruz Rodriguez (ISP-SU,2022).

## Dependencies
These routines make use of Numpy and Numba. They can be used without Numba by commenting out all the numba function decorators and by replacing numba.prange by range. But they become extremely slow.


## Usage
We have included example code showing how to use every function, please check the examply.py file in the example folder.


## Acknowledgements
I greatefulyl acknowledge funding from the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme (SUNMAG, grant agreement 759548).
