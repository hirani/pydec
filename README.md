## **PyDEC: A Python Library for Discretizations of Exterior Calculus.**

PyDEC is a Python library implementing Discrete Exterior Calculus (DEC) and lowest order finite element exterior calculus (FEEC) i.e., Whitney forms. The main functionality implemented:
- simplicial complexes of dimension n embeded in dimension N >= n
- abstract simplicial complexes (no embedding needed)
- cubical complexes in all dimensions
- boundary operators for all the above type of complexes
- discrete exterior derivative (i.e., coboundary operator)
- DEC discrete Hodge star (i.e., primal-dual diagonal mass matrix)
- FEEC mass matrix for Whitney forms
  
The code and companion paper include examples for numerically solving PDEs and computing cohomology: ACM Transactions on Mathematical Software, Vol. 39, No. 1, pp. 3:1-3:41, 2012, [DOI: 10.1145/2382585.2382588](http://dx.doi.org/10.1145/2382585.2382588).

Installation:
- `cd` to the folder where you cloned PyDEC
- `pip install .`

Note about installation:
- There is another package called `pydec` which has nothing to do with this package. If you simply do `pip install pydec` then you will get that other package. PyDEC has existed since about 2008, but we were too slow in grabbing `pydec` name on pip :-(
- So for now, please install by first cloning or downloading the source and then using the installation instructions above.

