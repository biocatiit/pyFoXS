"""
\file IMP/saxs/FitParameters.h \brief

Copyright 2007-2022 IMP Inventors. All rights reserved.
"""

import sys

class FitParameters:
    def __init__(self, chi_square=0.0, c1=0.0, c2=0.0,
                 c=0.0, o=0.0, default_chi_square=0.0):
        self.chi_square = chi_square
        self.c1 = c1
        self.c2 = c2
        self.c = c
        self.o = o
        self.default_chi_square = default_chi_square
        self.profile_file_name = None
        self.pdb_file_name = None
        self.mol_index = 0

    def show(self, stream=sys.stdout):
        if self.pdb_file_name is not None:
            stream.write(self.pdb_file_name + " ")
        if self.profile_file_name is not None:
            stream.write(self.profile_file_name + " ")
        stream.write("Chi^2 = " + str(self.chi_square) + " c1 = " + str(self.c1) +
                     " c2 = " + str(self.c2) + " default chi^2 = " +
                     str(self.default_chi_square) + "\n")

    def __lt__(self, other):
        return self.chi_square < other.chi_square
