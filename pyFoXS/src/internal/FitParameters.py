"""
\file IMP/saxs/FitParameters.h \brief

Copyright 2007-2022 IMP Inventors. All rights reserved.
"""

from typing import Optional
import sys

class FitParameters:
    def __init__(self, chi_square: float = 0.0, c1: float = 0.0, c2: float = 0.0,
                 c: float = 0.0, o: float = 0.0, default_chi_square: float = 0.0):
        self.chi_square = chi_square
        self.c1 = c1
        self.c2 = c2
        self.c = c
        self.o = o
        self.default_chi_square = default_chi_square
        self.profile_file_name: Optional[str] = None
        self.pdb_file_name: Optional[str] = None
        self.mol_index: int = 0

    def get_score(self) -> float:
        return self.chi_square

    def get_chi_square(self) -> float:
        return self.chi_square

    def get_c1(self) -> float:
        return self.c1

    def get_c2(self) -> float:
        return self.c2

    def get_scale(self) -> float:
        return self.c

    def get_offset(self) -> float:
        return self.o

    def get_default_chi_square(self) -> float:
        return self.default_chi_square

    def get_profile_file_name(self) -> Optional[str]:
        return self.profile_file_name

    def get_pdb_file_name(self) -> Optional[str]:
        return self.pdb_file_name

    def get_mol_index(self) -> int:
        return self.mol_index

    def set_score(self, score: float):
        self.chi_square = score

    def set_chi_square(self, chi_square: float):
        self.chi_square = chi_square

    def set_default_chi_square(self, chi_square: float):
        self.default_chi_square = chi_square

    def set_profile_file_name(self, file_name: str):
        self.profile_file_name = file_name

    def set_pdb_file_name(self, file_name: str):
        self.pdb_file_name = file_name

    def set_mol_index(self, index: int):
        self.mol_index = index

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
