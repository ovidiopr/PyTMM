PYTHON=`which python`
CYTHON=`which cython`
DESTDIR=/
PROJECT=python-tmmnlay
VERSION=1.0.0

all:
	@echo "make source - Create source package for Cython extension"
	@echo "make install - Install Cython extension on local system"
	@echo "make ext - Create Cython extension in place"
	@echo "make deb - Generate a deb package for Cython extension"
	@echo "make rpm - Generate a rpm package for Cython extension"
	@echo "make clean - Delete temporal files"

source:
	# build the source package in the parent directory
	$(PYTHON) setup.py sdist $(COMPILE) --dist-dir=../ --prune
	rm -rf MANIFEST

install:
	$(PYTHON) setup.py install --root $(DESTDIR) $(COMPILE)

ext: $(CURDIR)/tmmnlay.pyx
	$(PYTHON) setup.py build_ext --inplace

deb:
	# build the source package in the parent directory
	# then rename it to project_version.orig.tar.gz
	$(PYTHON) setup.py sdist $(COMPILE) --dist-dir=../ --prune
	rename -f 's/$(PROJECT)-(.*)\.tar\.gz/$(PROJECT)_$$1\.orig\.tar\.gz/' ../*
	# build the package
	dpkg-buildpackage -i -I -rfakeroot

rpm:
	$(PYTHON) setup.py bdist_rpm --dist-dir=../

clean:
	$(PYTHON) setup.py clean
	$(MAKE) -f $(CURDIR)/debian/rules clean
	rm -rf build/ dist/ MANIFEST
	find . -name '*.pyc' -delete
	find . -name '*.c' -delete
	find . -name '*.o' -delete
	find . -name '*.so' -delete
