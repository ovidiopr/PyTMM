# tmmnlay
Cython implementation of the Transfer Matrix Method

This extension implements the Transfer Matrix Method, which is
used in optics and acoustics to analyze the propagation of
electromagnetic or acoustic waves through a multilayer medium.
For details of the implementation, check the [algorithm](https://github.com/ovidiopr/tmmnlay/wiki/Algorithm).

## Credits
tmmnlay has been forked from [PyTMM](https://kitchenknif.github.io/PyTMM).
However, I have rewritten all the code from scratch in order to (i) increase
the calculation speed and (ii) make it more 'Pythonic' (e.g., with the new code
you can run the calculation for all wavelengths in a single step). In the process
of rewriting it I have also implemented functionality that was missing in PyTMM,
like the calculation of the electric field in the multilayer.

# Using tmmnlay

## Compiling the Code:
To compile the Cython extension you need Cython, a C compiler, [NumPy](http://www.numpy.org/):

 - **cython (>=0.23.4)**
 - **python-numpy (>= 1.0.3)**
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
from tmmnlay import MultiLayer
...
b = MultiLayer(n=(1.0, n2, n1), d=(0.0, d, 0.0), wvl=wavelength)
r, t = b.rt_TE
...
  ```

## License

GPL v3+
