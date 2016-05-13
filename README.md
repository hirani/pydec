## **PyDEC: A Python Library for Discretizations of Exterior Calculus.**

Exterior calculus is the generalization of vector calculus to manifolds. PyDEC is a Python library for computations related to the discretization of exterior calculus which includes numerical solution of partial differential equations. It is also useful for purely topological computations. Thus PyDEC facilitates inquiry into both physical problems on manifolds as well as purely topological problems on abstract complexes. It uses efficient algorithms for constructing the operators and objects and related topological problems. Our algorithms are formulated in terms of high-level matrix operations which extend to arbitrary dimension. As a result, our implementations map well to the facilities of numerical libraries such as NumPy and SciPy. The availability of such libraries makes Python suitable for prototyping numerical methods. The code and the companion paper includes examples where we demonstrate how PyDEC is used to solve physical and topological problems.

A companion paper, with many examples, is in the ACM Transactions on Mathematical Software, Vol. 39, No. 1, pp. 3:1-3:41, 2012, [DOI: 10.1145/2382585.2382588](http://dx.doi.org/10.1145/2382585.2382588). If you use PyDEC in your work please consider sending us a note. The package can be cited using the following BibTeX entry:
```
@webpage{BeHi2008,
    Author = {Nathan Bell and Anil N. Hirani},
    Note = {Software made available on Google Code website.},
    Title = {PyDEC: A {P}ython Library for {D}iscrete {E}xterior {C}alculus},
    Url = {https://github.com/hirani/pydec},
```
