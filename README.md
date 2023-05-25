# pyFoXS
Python version of the FoXS software

\brief Determine small angle X-ray (SAXS) profiles.

## Installation

First, you need to install the IMP library:
```
conda install -c conda-forge imp
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

## foxs {#foxs_bin}

Determine small angle X-ray (SAXS) profiles.
The IMP.saxs module contains functions that, given an atomic protein structure,
can calculate its SAXS profile using the Debye formula, and then fit this
profile against the experimentally determined one. FoXS is a simple command
line interface to this functionality which takes as input a number of PDB
files and/or SAXS profiles. There is also a \salilab{foxs/,web server}
available.

_Examples_:
 - [Determination of a Nup133 structure](@ref foxs_nup133)

## Info

_Author(s)_: Dina Schneidman

_Maintainer_: `duhovka`

_License_: [LGPL](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html)
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

_Publications_:
 - Dina Schneidman-Duhovny, Michael Hammel, Andrej Sali, \quote{FoXS: A Web server for Rapid Computation and Fitting of SAXS Profiles}, <em>Nucleic Acids Research</em>, 2010.
 - Dina Schneidman-Duhovny, Michael Hammel, John A. Tainer, Andrej Sali, \quote{Accurate SAXS profile computation and its assessment by contrast variation experiments}, <em> Biophysical Journal </em>, 2013.
