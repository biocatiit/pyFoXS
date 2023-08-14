# pyFoXS
Python version of the FoXS software

This is a python port of the open source program FoXS
(https://github.com/salilab/imp/tree/develop/modules/foxs) in the
Integrated Modelling Platform (https://integrativemodeling.org/). It is used
to determine theoretical small angle x-ray scattering (SAXS) profiles
from atomic models.

The python version is still in development, so use with caution!

## Install and run

First, you need to install the scipy and numpy library:
```
pip install scipy numpy
```

You can then run the foxs.py file like a normal python file:
```
python foxs.py examples/nup133/3KFO.pdb examples/nup133/23922_merge.dat
```

## IMP Code and build

This is not needed, it is only if you want to access the C++ code.

First, you need to download the IMP library. Two possibilities:

- Download the source code tarball from our download page, then extract it with something like:
`tar -xvzf ../imp-<version>.tar.gz`
- Alternatively you can use git to get the code directly from our GitHub repository with something like:
```
git clone -b main https://github.com/salilab/imp.git
(cd imp && git submodule update --init && ./setup_git.py)
```

Then, you need to build/install it for python:

To build IMP source found in `path/to/imp-source` and install it in
`path_to_install` do:

1. `mkdir build && cd build`
2. `cmake path/to/imp-source -DCMAKE_INSTALL_PREFIX=path_to_install`
3. `make -j4`
4. `make install`

Everything about IMP and how to download/compile it is here: [https://integrativemodeling.org/nightly/doc/manual/installation.html](https://integrativemodeling.org/nightly/doc/manual/installation.html).

Then you can run pyFoXS/foxs.py like a normal python file.
