# tmmnlay
Cython implementarion of the Transfer Matrix Method

This extension implements the Transfer Matrix Method, which is
used in optics and acoustics to analyze the propagation of
electromagnetic or acoustic waves through a multilayer medium.

## Credits
tmmnlay has been forked from [PyTMM](https://kitchenknif.github.io/PyTMM).

# Using tmmnlay

## Compiling the Code:
To compile the Cython extension you need Cython, a C compiler, [NumPy](http://www.numpy.org/):

 - **cython (>=0.23.4)**
 - **python-numpy (>= 1.0.3)**
 - **python-scipy (>=0.5.2)**
 - **python-all-dev (any version)**
 - **python-numpy-dev (any version)**

And to compile the Debian package you need some additional tools:

 - **debhelper (>=7.0.0)**
 - **dh-python (any version)**
 - **cdbs (>= 0.4.49)**

Compilation options

 - **make source** - Create source package for Cython extension
 - **make install** - Install Cython extension on local system
 - **make ext** - Create Cython extension in place
 - **make deb** - Generate a deb package for Cython extension
 - **make rpm** - Generate a rpm package for Cython extension
 - **make clean** - Delete temporal files

## Using tmmnlay
  
  ```python
from tmmnlay import TransferMatrix, solvePropagation
...
a = TransferMatrix.boundingLayer(1, n1)

R, T = solvePropagation(a)
...
  ```

## License

GPL v3+
