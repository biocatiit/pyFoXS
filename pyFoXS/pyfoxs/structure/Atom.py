"""
\file Atom.cpp   \brief Simple atoms decorator.

Copyright 2007-2022 IMP Inventors. All rights reserved.
"""

import numpy as np
from enum import IntEnum
from scipy.spatial import distance

from .Residue import Residue

atom_entry_type_field_ = 0
atom_number_field_ = 6
atom_type_field_ = 12
atom_alt_loc_field_ = 16
atom_res_name_field_ = 17
atom_chain_id_field_ = 21
atom_res_number_field_ = 22
atom_res_insertion_field_ = 26
atom_xcoord_field_ = 30
atom_ycoord_field_ = 38
atom_zcoord_field_ = 46
atom_occupancy_field_ = 54
atom_temp_factor_field_ = 60
atom_element_field_ = 76
model_index_field_ = 6

class Element(IntEnum):
    UNKNOWN_ELEMENT = 0
    OH = -1
    H2O = -2
    H = 1
    He = 2
    Li = 3
    Be = 4
    B = 5
    C = 6
    N = 7
    O = 8
    F = 9
    Ne = 10
    Na = 11
    Mg = 12
    Al = 13
    Si = 14
    P = 15
    S = 16
    Cl = 17
    Ar = 18
    K = 19
    Ca = 20
    Sc = 21
    Ti = 22
    V = 23
    Cr = 24
    Mn = 25
    Fe = 26
    Co = 27
    Ni = 28
    Cu = 29
    Zn = 30
    Ga = 31
    Ge = 32
    As = 33
    Se = 34
    Br = 35
    Kr = 36
    Rb = 37
    Sr = 38
    Y = 39
    Zr = 40
    Nb = 41
    Mo = 42
    Tc = 43
    Ru = 44
    Rh = 45
    Pd = 46
    Ag = 47
    Cd = 48
    In = 49
    Sn = 50
    Sb = 51
    Te = 52
    I = 53
    Xe = 54
    Cs = 55
    Ba = 56
    La = 57
    Ce = 58
    Pr = 59
    Nd = 60
    Pm = 61
    Sm = 62
    Eu = 63
    Gd = 64
    Tb = 65
    Dy = 66
    Ho = 67
    Er = 68
    Tm = 69
    Yb = 70
    Lu = 71
    Hf = 72
    Ta = 73
    W = 74
    Re = 75
    Os = 76
    Ir = 77
    Pt = 78
    Au = 79
    Hg = 80
    Tl = 81
    Pb = 82
    Bi = 83
    Po = 84
    At = 85
    Rn = 86
    Fr = 87
    Ra = 88
    Ac = 89
    Th = 90
    Pa = 91
    U = 92
    Np = 93
    Pu = 94
    Am = 95
    Cm = 96
    Bk = 97
    Cf = 98
    Es = 99
    Fm = 100
    Md = 101
    No = 102
    Lr = 103
    Db = 104
    Jl = 105
    Rf = 106

ElementType = {Element.UNKNOWN_ELEMENT : 'UNKNOWN',
                Element.H : 'H',
                Element.C : 'C',
                Element.N : 'N',
                Element.O : 'O',
                Element.S : 'S',
                Element.Se : 'SE',
                Element.P : 'P'}

AtomType = {'N' : Element.N,
            'H' : Element.H,
            'C' : Element.C,
            'O' : Element.O,
            'S' : Element.S,
            'SE' : Element.Se,
            'OXT' : Element.O,
            'OT1' : Element.O,
            'OT2' : Element.O,
            'CH3' : Element.C,
            'CH' : Element.C,
            'CA' : Element.C,
            'HA' : Element.H,
            'HA1' : Element.H,
            'HA2' : Element.H,
            'HA3' : Element.H,
            'CB' : Element.C,
            'HB' : Element.H,
            'HB1' : Element.H,
            'HB2' : Element.H,
            'HB3' : Element.H,
            'CG' : Element.C,
            'CG1' : Element.C,
            'CG2' : Element.C,
            'HG' : Element.H,
            'HG1' : Element.H,
            'HG2' : Element.H,
            'HG3' : Element.H,
            'HG11' : Element.H,
            'HG21' : Element.H,
            'HG31' : Element.H,
            'HG12' : Element.H,
            'HG13' : Element.H,
            'HG22' : Element.H,
            'HG23' : Element.H,
            'HG32' : Element.H,
            'OG' : Element.O,
            'OG1' : Element.O,
            'SG' : Element.S,
            'CD' : Element.C,
            'CD1' : Element.C,
            'CD2' : Element.C,
            'HD' : Element.H,
            'HD1' : Element.H,
            'HD2' : Element.H,
            'HD3' : Element.H,
            'HD11' : Element.H,
            'HD21' : Element.H,
            'HD31' : Element.H,
            'HD12' : Element.H,
            'HD13' : Element.H,
            'HD22' : Element.H,
            'HD23' : Element.H,
            'HD32' : Element.H,
            'SD' : Element.S,
            'OD1' : Element.O,
            'OD2' : Element.O,
            'ND1' : Element.N,
            'ND2' : Element.N,
            'CE' : Element.C,
            'CE1' : Element.C,
            'CE2' : Element.C,
            'CE3' : Element.C,
            'HE' : Element.H,
            'HE1' : Element.H,
            'HE2' : Element.H,
            'HE3' : Element.H,
            'HE21' : Element.H,
            'HE22' : Element.H,
            'OE1' : Element.O,
            'OE2' : Element.O,
            'NE' : Element.N,
            'NE1' : Element.N,
            'NE2' : Element.N,
            'CZ' : Element.C,
            'CZ2' : Element.C,
            'CZ3' : Element.C,
            'NZ' : Element.N,
            'HZ' : Element.H,
            'HZ1' : Element.H,
            'HZ2' : Element.H,
            'HZ3' : Element.H,
            'CH2' : Element.C,
            'NH1' : Element.N,
            'NH2' : Element.N,
            'OH' : Element.O,
            'HH' : Element.H,
            'HH11' : Element.H,
            'HH21' : Element.H,
            'HH2' : Element.H,
            'HH12' : Element.H,
            'HH22' : Element.H,
            'HH23' : Element.H,
            'HH33' : Element.H,
            'HH13' : Element.H,
            'P' : Element.P,
            'OP1' : Element.O,
            'OP2' : Element.O,
            'OP3' : Element.O,
            'O5p' : Element.O,
            'C5p' : Element.C,
            'H5pp' : Element.H,
            'C4p' : Element.C,
            'H4p' : Element.H,
            'H5p' : Element.H,
            'O4p' : Element.O,
            'C1p' : Element.C,
            'H1p' : Element.H,
            'C3p' : Element.C,
            'H3p' : Element.H,
            'O3p' : Element.O,
            'C2p' : Element.C,
            'H2p' : Element.H,
            'H2pp' : Element.H,
            'O2p' : Element.O,
            'HO2p' : Element.H,
            'N9' : Element.N,
            'C8' : Element.C,
            'H8' : Element.H,
            'N7' : Element.N,
            'C5' : Element.C,
            'C4' : Element.C,
            'N3' : Element.N,
            'C2' : Element.C,
            'N1' : Element.N,
            'C6' : Element.C,
            'N6' : Element.N,
            'H61' : Element.H,
            'H62' : Element.H,
            'O6' : Element.O,
            'N2' : Element.N,
            'NT' : Element.N,
            'H21' : Element.H,
            'H22' : Element.H,
            'H6' : Element.H,
            'H5' : Element.H,
            'O2' : Element.O,
            'N4' : Element.N,
            'H41' : Element.H,
            'H42' : Element.H,
            'O4' : Element.O,
            'C7' : Element.C,
            'H71' : Element.H,
            'H72' : Element.H,
            'H73' : Element.H,
            'O1A' : Element.O,
            'O2A' : Element.O,
            'O3A' : Element.O,
            'O1B' : Element.O,
            'O2B' : Element.O,
            'O3B' : Element.O,
            'CAY' : Element.C,
            'CY' : Element.C,
            'OY' : Element.O,
            'CAT' : Element.C,
            'NO2' : Element.N,
            'UNKNOWN' : Element.UNKNOWN_ELEMENT}

def add_atom_type(name, e):
    if name in AtomType:
        raise Exception("An AtomType with that name already exists: " + name)
    AtomType[name] = AtomType[e]
    return name

def get_element_for_atom_type(at):
    if at not in AtomType:
        raise Exception("Invalid AtomType index: " + str(at))
    return ElementType[AtomType[at]]

def get_atom_type_exists(name):
    return name in AtomType

def get_residue(d, nothrow=False):
    mhd = d
    while mhd is not None:
        mhd = mhd.parent
        if isinstance(mhd, Residue):
            return mhd
    if nothrow:
        return Residue()
    else:
        print("Atom is not the child of a residue: " + str(d))

class Atom:
    def __init__(self, name='', type='', element='', atom_index=0, coord=(0,0,0), occupancy=None, temp_factor=None):
        self.atom_type = type
        # self.element = element
        self.atom_index = atom_index
        self.coordinates = coord
        self.occupancy = occupancy
        self.temp_factor = temp_factor
        self.atom_name = name
        self.parent = None
        self.radius = None
        self.residue = None
        self.cache_attributes = {}

    def setup_particle(self, radius):
        self.radius = radius

    @staticmethod
    def get_atom_type_key():
        return "atom_type"

    @staticmethod
    def get_element_key():
        return "element"

    @staticmethod
    def get_input_index_key():
        return "pdb_atom_index"

    @staticmethod
    def get_occupancy_key():
        return "occupancy"

    @staticmethod
    def get_temperature_factor_key():
        return "tempFactor"

class Chain:
    def __init__(self, chain_id):
        self.name = "Chain"
        self.chain_id = chain_id
        self.parent = None
        self.residues = []

def compute_max_distance(particles):
    max_dist2 = 0.0
    coordinates = [particle.coordinates for particle in particles]

    for i in range(len(coordinates)):
        for j in range(i + 1, len(coordinates)):
            dist2 = distance.sqeuclidean(coordinates[i], coordinates[j])
            if dist2 > max_dist2:
                max_dist2 = dist2

    return np.sqrt(max_dist2)

def compute_max_distance_between_particles(particles1, particles2):
    max_dist2 = 0.0
    coordinates1 = [particle.coordinates for particle in particles1]
    coordinates2 = [particle.coordinates for particle in particles2]

    for i in range(len(coordinates1)):
        for j in range(len(coordinates2)):
            dist2 = distance.sqeuclidean(coordinates1[i], coordinates2[j])
            if dist2 > max_dist2:
                max_dist2 = dist2

    return np.sqrt(max_dist2)
