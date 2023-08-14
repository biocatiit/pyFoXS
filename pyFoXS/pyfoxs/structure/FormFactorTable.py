"""
\file IMP/saxs/FormFactorTable.h   \brief A class for computation of
atomic and residue level form factors for SAXS calculations

Copyright 2007-2022 IMP Inventors. All rights reserved.
"""

import enum
import math

from .Residue import ResidueType, residue_to_string

class FormFactorType(enum.Enum):
    ALL_ATOMS = 0
    HEAVY_ATOMS = 1
    CA_ATOMS = 2
    RESIDUES = 3

class FormFactor:
    def __init__(self, ff=0, vacuum_ff=0, dummy_ff=0):
        self.ff_ = ff
        self.vacuum_ff_ = vacuum_ff
        self.dummy_ff_ = dummy_ff

class FormFactorTable:
    # electron density of solvent - default=0.334 e/A^3 (H2O)
    rho_ = 0.334

    class FormFactorAtomType(enum.Enum):
        H = 0
        He = 1
        Li = 2
        Be = 3
        B = 4
        C = 5
        N = 6
        O = 7
        F = 8
        Ne = 9
        Na = 10
        Mg = 11
        Al = 12
        Si = 13
        P = 14
        S = 15
        Cl = 16
        Ar = 17
        K = 18
        Ca = 19
        Cr = 20
        Mn = 21
        Fe = 22
        Co = 23
        Ni = 24
        Cu = 25
        Zn = 26
        Se = 27
        Br = 28
        Ag = 29
        I = 30
        Ir = 31
        Pt = 32
        Au = 33
        Hg = 34
        ALL_ATOM_SIZE = 35
        CH = 35
        CH2 = 36
        CH3 = 37
        NH = 38
        NH2 = 39
        NH3 = 40
        OH = 41
        OH2 = 42
        SH = 43
        HEAVY_ATOM_SIZE = 44
        UNK = 45

    element_ff_type_map_ = {member.name.upper(): member for member in FormFactorAtomType}
    residue_type_form_factor_map_ = {
        ResidueType.ALA : FormFactor(9.037, 37.991, 28.954),
        ResidueType.ARG : FormFactor(23.289, 84.972, 61.683),
        ResidueType.ASP : FormFactor(20.165, 58.989, 38.824),
        ResidueType.ASN : FormFactor(19.938, 59.985, 40.047),
        ResidueType.CYS : FormFactor(18.403, 53.991, 35.588),
        ResidueType.GLN : FormFactor(19.006, 67.984, 48.978),
        ResidueType.GLU : FormFactor(19.233, 66.989, 47.755),
        ResidueType.GLY : FormFactor(10.689, 28.992, 18.303),
        ResidueType.HIS : FormFactor(21.235, 78.977, 57.742),
        ResidueType.ILE : FormFactor(6.241, 61.989, 55.748),
        ResidueType.LEU : FormFactor(6.241, 61.989, 55.748),
        ResidueType.LYS : FormFactor(10.963, 70.983, 60.020),
        ResidueType.MET : FormFactor(16.539, 69.989, 53.450),
        ResidueType.PHE : FormFactor(9.206, 77.986, 68.7806),
        ResidueType.PRO : FormFactor(8.613, 51.9897, 43.377),
        ResidueType.SER : FormFactor(13.987, 45.991, 32.004),
        ResidueType.THR : FormFactor(13.055, 53.99, 40.935),
        ResidueType.TYR : FormFactor(14.156, 85.986, 71.83),
        ResidueType.TRP : FormFactor(14.945, 98.979, 84.034),
        ResidueType.VAL : FormFactor(7.173, 53.9896, 46.817),
        ResidueType.POP : FormFactor(45.616, 365.99, 320.41),
        ResidueType.UNK : FormFactor(9.037, 37.991, 28.954)
    }

    zero_form_factors_ = [
        -0.720147, -0.720228,
        #   H       He - periodic table line 1
        1.591, 2.591, 3.591, 0.50824, 6.16294, 4.94998, 7.591, 6.993,
        # Li     Be      B     C       N        O       F      Ne - line 2
        7.9864, 8.9805, 9.984, 10.984, 13.0855, 9.36656, 13.984, 16.591,
        #  Na      Mg        Al       Si        P        S       Cl    Ar - line 3
        15.984, 14.9965, 20.984, 21.984, 20.9946, 23.984,
        # K       Ca2+       Cr      Mn      Fe2+      Co - line 4
        24.984, 25.984, 24.9936, 30.9825, 31.984, 43.984, 49.16,
        # Ni     Cu          Zn2+      Se       Br       Ag      I
        70.35676, 71.35676, 72.324, 73.35676,
        # Ir         Pt      Au      Hg
        -0.211907, -0.932054, -1.6522, 5.44279, 4.72265, 4.0025, 4.22983, 3.50968,
        8.64641
        #  CH        CH2        CH3     NH       NH2       NH3     OH       OH2
        # SH
    ]

    vacuum_zero_form_factors_ = [
        #   H       He - periodic table line 1
        0.999953, 0.999872,
        # Li  Be    B     C       N       O       F     Ne - line 2
        2.99, 3.99, 4.99, 5.9992, 6.9946, 7.9994, 8.99, 9.999,
        #  Na     Mg     Al     Si      P        S        Cl     Ar - line 3
        10.9924, 11.9865, 12.99, 13.99, 14.9993, 15.9998, 16.99, 17.99,
        # K    Ca2+     Cr     Mn     Fe2+     Co - line 4
        18.99, 18.0025, 23.99, 24.99, 24.0006, 26.99,
        # Ni   Cu      Zn2+     Se     Br - line 4 cont.
        27.99, 28.99, 27.9996, 33.99, 34.99,
        # Ag    I       Ir     Pt      Au     Hg - some elements from lines 5, 6
        46.99, 52.99, 76.99, 77.99, 78.9572, 79.99,
        # CH      CH2     CH3     NH       NH2       NH3     OH      OH2      SH
        6.99915, 7.99911, 8.99906, 7.99455, 8.99451, 9.99446, 8.99935, 9.9993,
        16.9998
    ]

    dummy_zero_form_factors_ = [
        1.7201, 1.7201, 1.399, 1.399, 1.399, 5.49096, 0.83166, 3.04942,
        1.399, 3.006,
        #  H     He     Li?    Be?    B?       C        N        O      F?     Ne
        3.006, 3.006, 3.006, 3.006, 1.91382, 6.63324, 3.006, 1.399,
        # Na     Mg    Al?    Si?      P        S      Cl?    Ar?
        3.006, 3.006, 3.006, 3.006, 3.006, 3.006,
        # K?   Ca2+    Cr?    Mn?   Fe2+   Co?
        3.006, 3.006, 3.006, 3.006, 3.006,
        # Ni?   Cu?   Zn2+    Se     Br?
        3.006, 3.83, 6.63324, 6.63324, 6.63324, 6.63324,
        # Ag?   I?       Ir?      Pt?       Au      Hg
        7.21106, 8.93116, 10.6513, 2.55176, 4.27186, 5.99196, 4.76952, 6.48962,
        8.35334
        #  CH       CH2      CH3     NH       NH2       NH3     OH       OH2   SH
    ]

    # form_factor_type_key_ = IntKey()

    form_factors_coefficients_ = []
    form_factors_ = []
    vacuum_form_factors_ = []
    dummy_form_factors_ = []
    min_q_ = 0.0
    max_q_ = 0.0
    delta_q_ = 0.0
    # warn_context_ = WarningContext()

    class AtomFactorCoefficients:
        def __init__(self, atom_type, a, c, b, excl_vol):
            self.atom_type_ = atom_type
            self.a_ = a
            self.c_ = c
            self.b_ = b
            self.excl_vol_ = excl_vol

        def __str__(self):
            return (
                f"Atom Type: {self.atom_type_}\n"
                f"a values: {self.a_}\n"
                f"c value: {self.c_}\n"
                f"b values: {self.b_}\n"
                f"Excluded Volume: {self.excl_vol_}"
            )

    def __init__(self, table_name="", min_q=0.0, max_q=0.0, delta_q=0.0):
        if table_name:
            self.read_form_factor_table(table_name)
        self.min_q_ = min_q
        self.max_q_ = max_q
        self.delta_q_ = delta_q
        self.dummy_form_factors_ = self.dummy_zero_form_factors_
        self.vacuum_form_factors_ = self.vacuum_zero_form_factors_
        self.zero_form_factors_ = self.zero_form_factors_
        self.form_factor_type_key_ = "form_factor_type_key"


    def get_form_factor(self, p, ff_type=FormFactorType.HEAVY_ATOMS):
        if ff_type in (FormFactorType.CA_ATOMS, FormFactorType.RESIDUES):  # residue level form factors
            residue_type = p.residue.residue_type
            return self.get_form_factor_r(residue_type)

        # atomic form factor, initialization by request
        if self.form_factor_type_key_ in p.cache_attributes:
            return self.zero_form_factors_[p.cache_attributes[self.form_factor_type_key_].value]

        ff_atom_type = self.get_form_factor_atom_type(p, ff_type)
        if ff_atom_type.value >= self.FormFactorAtomType.HEAVY_ATOM_SIZE.value:
            print("Can't find form factor for particle",
                  p.atom_type,
                  "using default")
            ff_atom_type = self.FormFactorAtomType.N

        form_factor = self.zero_form_factors_[ff_atom_type.value]
        p.cache_attributes[self.form_factor_type_key_] = ff_atom_type
        return form_factor

    def get_form_factor_r(self, rt):
        if rt in self.residue_type_form_factor_map_:
            return self.residue_type_form_factor_map_[rt].ff_
        else:
            print("Can't find form factor for residue", rt,
                  "using default value of ALA")
            return self.residue_type_form_factor_map_[self.FormFactorAtomType.UNK].ff_

    def get_vacuum_form_factor(self, p, ff_type):
        if ff_type == FormFactorType.CA_ATOMS:  # residue level form factors
            residue_type = p.residue.residue_type if p.residue else None
            if residue_type:
                return self.get_vacuum_form_factor_r(residue_type)

        if ff_type == FormFactorType.RESIDUES:  # residue level form factors
            residue_type = p.residue.residue_type if p.residue else None
            if residue_type:
                return self.get_form_factor_r(residue_type)

        if self.form_factor_type_key_ in p.cache_attributes:
            return self.vacuum_zero_form_factors_[p.cache_attributes[self.form_factor_type_key_].value]

        ff_atom_type = self.get_form_factor_atom_type(p, ff_type)
        form_factor = self.vacuum_zero_form_factors_[ff_atom_type.value]
        p.cache_attributes[self.form_factor_type_key_] = ff_atom_type
        return form_factor

    def get_vacuum_form_factor_r(self, rt):
        if rt in self.residue_type_form_factor_map_:
            return self.residue_type_form_factor_map_[rt].vacuum_ff_
        else:
            print("Can't find form factor for residue", rt, "using default value of ALA")
            return self.residue_type_form_factor_map_[ResidueType.UNK].vacuum_ff_

    def get_dummy_form_factor(self, p, ff_type):
        if ff_type == FormFactorType.CA_ATOMS:
            # Residue level form factors
            residue_type = p.residue.residue_type if p.residue else None
            if residue_type:
                return self.get_dummy_form_factor_r(residue_type)

        if ff_type == FormFactorType.RESIDUES:
            # Residue level form factors
            residue_type = p.residue.residue_type if p.residue else None
            if residue_type:
                return self.get_form_factor_r(residue_type)

        if self.form_factor_type_key_ in p.cache_attributes:
            return self.dummy_zero_form_factors_[p.cache_attributes[self.form_factor_type_key_].value]

        ff_atom_type = self.get_form_factor_atom_type(p, ff_type)
        form_factor = self.dummy_zero_form_factors_[ff_atom_type.value]
        p.cache_attributes[self.form_factor_type_key_] = ff_atom_type
        return form_factor

    def get_dummy_form_factor_r(self, rt):
        if rt in self.residue_type_form_factor_map_:
            return self.residue_type_form_factor_map_[rt].dummy_ff_
        else:
            print("Can't find form factor for residue ", rt,
                    " using default value of ALA")
            return self.residue_type_form_factor_map_[ResidueType.UNK].dummy_ff_

    def get_form_factor_atom_type(self, p, ff_type):
        ad = p
        residue_type = ad.residue.residue_type if ad.residue is not None else ""
        residue_type = residue_to_string(residue_type)
        atom_type = ad.atom_name
        # Find FormFactorAtomType
        ret_type = self.element_ff_type_map_[ad.atom_type] if ad.atom_type in self.element_ff_type_map_ else self.FormFactorAtomType.UNK

        if ff_type == FormFactorType.HEAVY_ATOMS:
            if ret_type == self.FormFactorAtomType.C:
                ret_type = self.get_carbon_atom_type(atom_type, residue_type)
            elif ret_type == self.FormFactorAtomType.N:
                ret_type = self.get_nitrogen_atom_type(atom_type, residue_type)
            elif ret_type == self.FormFactorAtomType.O:
                ret_type = self.get_oxygen_atom_type(atom_type, residue_type)
            elif ret_type == self.FormFactorAtomType.S:
                ret_type = self.get_sulfur_atom_type(atom_type, residue_type)

        if ret_type.value >= self.FormFactorAtomType.HEAVY_ATOM_SIZE.value:
            print("Can't find form factor for particle "
                    + ad.atom_type
                    + " using default value of nitrogen\n")
            ret_type = self.FormFactorAtomType.N

        return ret_type

    def get_water_form_factor(self):
        return self.zero_form_factors_[self.FormFactorAtomType.OH2.value]

    def get_vacuum_water_form_factor(self):
        return self.vacuum_zero_form_factors_[self.FormFactorAtomType.OH2.value]

    def get_dummy_water_form_factor(self):
        return self.dummy_zero_form_factors_[self.FormFactorAtomType.OH2.value]

    def get_form_factors(self, p, ff_type=FormFactorType.HEAVY_ATOMS):
        return self.get_form_factors(p, ff_type)

    def get_vacuum_form_factors(self, p, ff_type=FormFactorType.HEAVY_ATOMS):
        return self.get_vacuum_form_factors(p, ff_type)

    def get_dummy_form_factors(self, p, ff_type=FormFactorType.HEAVY_ATOMS):
        return self.get_dummy_form_factors(p, ff_type)

    def get_water_form_factors(self):
        return self.form_factors_[self.FormFactorAtomType.OH2.value]

    def get_water_vacuum_form_factors(self):
        return self.vacuum_form_factors_[self.FormFactorAtomType.OH2.value]

    def get_water_dummy_form_factors(self):
        return self.dummy_form_factors_[self.FormFactorAtomType.OH2.value]

    def get_radius(self, p, ff_type):
        # dummy_zero_form_factor = volume * rho
        # volume = 4/3 * pi * r^3
        # r^3 = 3*dummy_zero_form_factor / 4*pi*rho
        one_third = 1.0 / 3
        c = 3.0 / (4 * math.pi * self.rho_)
        form_factor = self.get_dummy_form_factor(p, ff_type)
        return math.pow(c * form_factor, one_third)

    def get_volume(self, p, ff_type):
        # dummy_zero_form_factor = volume * rho
        form_factor = self.get_dummy_form_factor(p, ff_type)
        return form_factor / self.rho_

    def show(self, out=print, prefix=""):
        self.show(out, prefix)

    def get_carbon_atom_type(self, atom_type, residue_type):
        # protein atoms
        # CH
        if atom_type == "CH":
            return self.FormFactorAtomType.CH
        # CH2
        if atom_type == "CH2":
            return self.FormFactorAtomType.CH2
        # CH3
        if atom_type == "CH3":
            return self.FormFactorAtomType.CH3
        # C
        if atom_type == "C":
            return self.FormFactorAtomType.C

        # CA
        if atom_type == "CA":
            if residue_type == "GLY":
                return self.FormFactorAtomType.CH2  # Glycine has 2 hydrogens
            return self.FormFactorAtomType.CH
        # CB
        if atom_type == "CB":
            if (residue_type == "ILE" or residue_type == "THR" or
                    residue_type == "VAL"):
                return self.FormFactorAtomType.CH
            if residue_type == "ALA":
                return self.FormFactorAtomType.CH3
            return self.FormFactorAtomType.CH2
        # CG1
        if atom_type == "CG":
            if (residue_type == "ASN" or residue_type == "ASP" or
                    residue_type == "HIS" or residue_type == "PHE" or
                    residue_type == "TRP" or residue_type == "TYR"):
                return self.FormFactorAtomType.C
            if residue_type == "LEU":
                return self.FormFactorAtomType.CH
            return self.FormFactorAtomType.CH2
        # CG1
        if atom_type == "CG1":
            if residue_type == "ILE":
                return self.FormFactorAtomType.CH2
            if residue_type == "VAL":
                return self.FormFactorAtomType.CH3
        # CG2 - only VAL, ILE, and THR
        if atom_type == "CG2":
            return self.FormFactorAtomType.CH3
        # CD
        if atom_type == "CD":
            if residue_type == "GLU" or residue_type == "GLN":
                return self.FormFactorAtomType.C
            return self.FormFactorAtomType.CH2
        # CD1
        if atom_type == "CD1":
            if residue_type == "LEU" or residue_type == "ILE":
                return self.FormFactorAtomType.CH3
            if (residue_type == "PHE" or residue_type == "TRP" or
                    residue_type == "TYR"):
                return self.FormFactorAtomType.CH
            return self.FormFactorAtomType.C
        # CD2
        if atom_type == "CD2":
            if residue_type == "LEU":
                return self.FormFactorAtomType.CH3
            if (residue_type == "PHE" or residue_type == "HIS" or
                    residue_type == "TYR"):
                return self.FormFactorAtomType.CH
            return self.FormFactorAtomType.C
        # CE
        if atom_type == "CE":
            if residue_type == "LYS":
                return self.FormFactorAtomType.CH2
            if residue_type == "MET":
                return self.FormFactorAtomType.CH3
            return self.FormFactorAtomType.C
        # CE1
        if atom_type == "CE1":
            if (residue_type == "PHE" or residue_type == "HIS" or
                    residue_type == "TYR"):
                return self.FormFactorAtomType.CH
            return self.FormFactorAtomType.C
        # CE2
        if atom_type == "CE2":
            if residue_type == "PHE" or residue_type == "TYR":
                return self.FormFactorAtomType.CH
            return self.FormFactorAtomType.C
        # CZ
        if atom_type == "CZ":
            if residue_type == "PHE":
                return self.FormFactorAtomType.CH
            return self.FormFactorAtomType.C
        # CZ2, CZ3, CE3
        if (atom_type == "CZ2" or atom_type == "CZ3" or
                atom_type == "CE3"):
            if residue_type == "TRP":
                return self.FormFactorAtomType.CH
            return self.FormFactorAtomType.C

        # DNA/RNA atoms
        # C5'
        if atom_type == "C5p":
            return self.FormFactorAtomType.CH2
        # C1', C2', C3', C4'
        if (atom_type == "C4p" or atom_type == "C3p" or
                atom_type == "C2p" or atom_type == "C1p"):
            return self.FormFactorAtomType.CH
        # C2
        if atom_type == "C2":
            if (residue_type == "DADE" or residue_type == "ADE"):
                return self.FormFactorAtomType.CH
            return self.FormFactorAtomType.C
        # C4
        if atom_type == "C4":
            return self.FormFactorAtomType.C
        # C5
        if atom_type == "C5":
            if (residue_type == "DCYT" or residue_type == "CYT" or
                    residue_type == "DURA" or residue_type == "URA"):
                return self.FormFactorAtomType.CH
            return self.FormFactorAtomType.C
        # C6
        if atom_type == "C6":
            if (residue_type == "DCYT" or residue_type == "CYT" or
                    residue_type == "DURA" or residue_type == "URA" or
                    residue_type == "DTHY" or residue_type == "THY"):
                return self.FormFactorAtomType.CH
            return self.FormFactorAtomType.C
        # C7
        if atom_type == "C7":
            return self.FormFactorAtomType.CH3
        # C8
        if atom_type == "C8":
            return self.FormFactorAtomType.CH

        print("Carbon atom not found, using default C form factor for "
                    + atom_type + " " + residue_type)
        return self.FormFactorAtomType.C

    def get_nitrogen_atom_type(self, atom_type, residue_type):
        # protein atoms
        # N
        if atom_type == "N":
            if residue_type == "PRO":
                return self.FormFactorAtomType.N
            return self.FormFactorAtomType.NH
        # ND1
        if atom_type == "ND1":
            if residue_type == "HIS":
                return self.FormFactorAtomType.NH
            return self.FormFactorAtomType.N
        # ND2
        if atom_type == "ND2":
            if residue_type == "ASN":
                return self.FormFactorAtomType.NH2
            return self.FormFactorAtomType.N
        # NH1, NH2
        if atom_type == "NH1" or atom_type == "NH2":
            if residue_type == "ARG":
                return self.FormFactorAtomType.NH2
            return self.FormFactorAtomType.N
        # NE
        if atom_type == "NE":
            if residue_type == "ARG":
                return self.FormFactorAtomType.NH
            return self.FormFactorAtomType.N
        # NE1
        if atom_type == "NE1":
            if residue_type == "TRP":
                return self.FormFactorAtomType.NH
            return self.FormFactorAtomType.N
        # NE2
        if atom_type == "NE2":
            if residue_type == "GLN":
                return self.FormFactorAtomType.NH2
            return self.FormFactorAtomType.N
        # NZ
        if atom_type == "NZ":
            if residue_type == "LYS":
                return self.FormFactorAtomType.NH3
            return self.FormFactorAtomType.N

        # DNA/RNA atoms
        # N1
        if atom_type == "N1":
            if residue_type == "DGUA" or residue_type == "GUA":
                return self.FormFactorAtomType.NH
            return self.FormFactorAtomType.N
        # N2, N4, N6
        if atom_type == "N2" or atom_type == "N4" or atom_type == "N6":
            return self.FormFactorAtomType.NH2
        # N3
        if atom_type == "N3":
            if residue_type == "DURA" or residue_type == "URA":
                return self.FormFactorAtomType.NH
            return self.FormFactorAtomType.N
        # N7, N9
        if atom_type == "N7" or atom_type == "N9":
            return self.FormFactorAtomType.N

        print(f"Nitrogen atom not found, using default N form factor for {atom_type} {residue_type}")

        return self.FormFactorAtomType.N

    def get_oxygen_atom_type(self, atom_type, residue_type):
        # O OE1 OE2 OD1 OD2 O1A O2A OXT OT1 OT2
        if atom_type == "O" or atom_type == "OE1" or \
                atom_type == "OE2" or atom_type == "OD1" or \
                atom_type == "OD2" or atom_type == "OXT":
            return self.FormFactorAtomType.O
        # OG
        if atom_type == "OG":
            if residue_type == "SER":
                return self.FormFactorAtomType.OH
            return self.FormFactorAtomType.O
        # OG1
        if atom_type == "OG1":
            if residue_type == "THR":
                return self.FormFactorAtomType.OH
            return self.FormFactorAtomType.O
        # OH
        if atom_type == "OH":
            if residue_type == "TYR":
                return self.FormFactorAtomType.OH
            return self.FormFactorAtomType.O

        # DNA/RNA atoms
        # O1P, O3', O2P, O2',O4',05', O2,O4,O6
        if atom_type == "OP1" or atom_type == "O3p "or \
                atom_type == "OP2" or atom_type == "O4p" or \
                atom_type == "O5p" or atom_type == "O2" or \
                atom_type == "O4" or atom_type == "O6":
            return self.FormFactorAtomType.O
        # O2'
        if atom_type == "O2p":
            return self.FormFactorAtomType.OH

        # water molecule
        if residue_type == "HOH":
            return self.FormFactorAtomType.OH2

        print(f"Oxygen atom not found, using default O form factor for {atom_type} {residue_type}")

        return self.FormFactorAtomType.O


    def get_sulfur_atom_type(self, atom_type, residue_type):
        # SD
        if atom_type == "SD":
            return self.FormFactorAtomType.S
        # SG
        if atom_type == "SG":
            if residue_type == "CYS":
                return self.FormFactorAtomType.SH
            return self.FormFactorAtomType.S

        print(f"Sulfur atom not found, using default S form factor for {atom_type} {residue_type}")

        return self.FormFactorAtomType.S


def get_default_form_factor_table():
    return FormFactorTable()

