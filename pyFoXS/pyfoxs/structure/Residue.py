"""
\file Residue.cpp   \brief Simple residue decorator.

Copyright 2007-2022 IMP Inventors. All rights reserved.
"""

from enum import IntEnum
import math

class ResidueType(IntEnum):
    UNK = 0
    GLY = 1
    ALA = 2
    VAL = 3
    LEU = 4
    ILE = 5
    SER = 6
    THR = 7
    CYS = 8
    MET = 9
    PRO = 10
    ASP = 11
    ASN = 12
    GLU = 13
    GLN = 14
    LYS = 15
    ARG = 16
    HIS = 17
    PHE = 18
    TYR = 19
    TRP = 20
    ACE = 21
    NH2 = 22
    MSE = 23
    ADE = 24
    URA = 25
    CYT = 26
    GUA = 27
    THY = 28
    DADE = 29
    DURA = 30
    DCYT = 31
    DGUA = 32
    DTHY = 33
    HOH = 34
    HEME = 35
    POP = 36

residue_types = {
        "UNK": ResidueType.UNK,
        "GLY": ResidueType.GLY,
        "ALA": ResidueType.ALA,
        "VAL": ResidueType.VAL,
        "LEU": ResidueType.LEU,
        "ILE": ResidueType.ILE,
        "SER": ResidueType.SER,
        "THR": ResidueType.THR,
        "CYS": ResidueType.CYS,
        "MET": ResidueType.MET,
        "PRO": ResidueType.PRO,
        "ASP": ResidueType.ASP,
        "ASN": ResidueType.ASN,
        "GLU": ResidueType.GLU,
        "GLN": ResidueType.GLN,
        "LYS": ResidueType.LYS,
        "ARG": ResidueType.ARG,
        "HIS": ResidueType.HIS,
        "PHE": ResidueType.PHE,
        "TYR": ResidueType.TYR,
        "TRP": ResidueType.TRP,
        "ACE": ResidueType.ACE,
        "NH2": ResidueType.NH2,
        "MSE": ResidueType.MSE,
        "ADE": ResidueType.ADE,
        "URA": ResidueType.URA,
        "CYT": ResidueType.CYT,
        "GUA": ResidueType.GUA,
        "THY": ResidueType.THY,
        "DADE": ResidueType.DADE,
        "DURA": ResidueType.DURA,
        "DCYT": ResidueType.DCYT,
        "DGUA": ResidueType.DGUA,
        "DTHY": ResidueType.DTHY,
        "A" : ResidueType.ADE,
        "U" : ResidueType.URA,
        "C" : ResidueType.CYT,
        "G" : ResidueType.GUA,
        "T" : ResidueType.THY,
        "DA": ResidueType.DADE,
        "DU": ResidueType.DURA,
        "DC": ResidueType.DCYT,
        "DG": ResidueType.DGUA,
        "DT": ResidueType.DTHY,
        "HOH": ResidueType.HOH,
        "HEME": ResidueType.HEME,
        "POP": ResidueType.POP
    }
def get_residue_type(code):
    return residue_types[code]

residue_string_types = {
        ResidueType.UNK : "UNK",
        ResidueType.GLY : "GLY",
        ResidueType.ALA : "ALA",
        ResidueType.VAL : "VAL",
        ResidueType.LEU : "LEU",
        ResidueType.ILE : "ILE",
        ResidueType.SER : "SER",
        ResidueType.THR : "THR",
        ResidueType.CYS : "CYS",
        ResidueType.MET : "MET",
        ResidueType.PRO : "PRO",
        ResidueType.ASP : "ASP",
        ResidueType.ASN : "ASN",
        ResidueType.GLU : "GLU",
        ResidueType.GLN : "GLN",
        ResidueType.LYS : "LYS",
        ResidueType.ARG : "ARG",
        ResidueType.HIS : "HIS",
        ResidueType.PHE : "PHE",
        ResidueType.TYR : "TYR",
        ResidueType.TRP : "TRP",
        ResidueType.ACE : "ACE",
        ResidueType.NH2 : "NH2",
        ResidueType.MSE : "MSE",
        ResidueType.ADE : "ADE",
        ResidueType.URA : "URA",
        ResidueType.CYT : "CYT",
        ResidueType.GUA : "GUA",
        ResidueType.THY : "THY",
        ResidueType.DADE : "DADE",
        ResidueType.DURA : "DURA",
        ResidueType.DCYT : "DCYT",
        ResidueType.DGUA : "DGUA",
        ResidueType.DTHY : "DTHY",
        ResidueType.HOH : "HOH",
        ResidueType.HEME : "HEME",
        ResidueType.POP : "POP"
    }
def residue_to_string(code):
    return residue_string_types[code]

one_letter_codes = {
        ResidueType.UNK: "X",
        ResidueType.GLY: "G",
        ResidueType.ALA: "A",
        ResidueType.VAL: "V",
        ResidueType.LEU: "L",
        ResidueType.ILE: "I",
        ResidueType.SER: "S",
        ResidueType.THR: "T",
        ResidueType.CYS: "C",
        ResidueType.MET: "M",
        ResidueType.PRO: "P",
        ResidueType.ASP: "D",
        ResidueType.ASN: "N",
        ResidueType.GLU: "E",
        ResidueType.GLN: "Q",
        ResidueType.LYS: "K",
        ResidueType.ARG: "R",
        ResidueType.HIS: "H",
        ResidueType.PHE: "F",
        ResidueType.TYR: "Y",
        ResidueType.TRP: "W",
        ResidueType.ACE: "ACE",
        ResidueType.NH2: "NH2",
        ResidueType.MSE: "MSE",
        ResidueType.ADE: "A",
        ResidueType.URA: "U",
        ResidueType.CYT: "C",
        ResidueType.GUA: "G",
        ResidueType.THY: "T",
        ResidueType.DADE: "DA",
        ResidueType.DURA: "DU",
        ResidueType.DCYT: "DC",
        ResidueType.DGUA: "DG",
        ResidueType.DTHY: "DT",
        ResidueType.HOH: "HOH",
        ResidueType.HEME: "HEME",
        ResidueType.POP: "POP"
    }

def get_one_letter_code(residue_type):
    return one_letter_codes.get(residue_type)

residue_masses = {
        ResidueType.ALA: 71.079,
        ResidueType.ARG: 156.188,
        ResidueType.ASP: 115.089,
        ResidueType.ASN: 114.104,
        ResidueType.CYS: 103.144,
        ResidueType.GLN: 128.131,
        ResidueType.GLU: 129.116,
        ResidueType.GLY: 57.052,
        ResidueType.HIS: 137.142,
        ResidueType.ILE: 113.160,
        ResidueType.LEU: 113.160,
        ResidueType.LYS: 128.174,
        ResidueType.MET: 131.198,
        ResidueType.PHE: 147.177,
        ResidueType.PRO: 97.117,
        ResidueType.SER: 87.078,
        ResidueType.THR: 101.105,
        ResidueType.TYR: 163.170,
        ResidueType.TRP: 186.213,
        ResidueType.VAL: 99.133,
        ResidueType.UNK: 113.160,
        ResidueType.ADE: 507.2,
        ResidueType.URA: 484.2,
        ResidueType.CYT: 483.2,
        ResidueType.GUA: 523.2,
        ResidueType.DADE: 491.2,
        ResidueType.DTHY: 482.2,
        ResidueType.DCYT: 467.2,
        ResidueType.DGUA: 507.2,
        ResidueType.POP: 686.97
    }

def get_mass(residue_type):
    return residue_masses.get(residue_type, math.nan)

class Residue:
    def __init__(self, residue_type, index=-1, insertion_code=" "):
        self.name = "Residue" + str(index)
        self.residue_type = residue_type
        self.index = index
        self.insertion_code = insertion_code
        self.parent = None
        self.atoms = []

    def add_child(self, atom):
        self.atoms.append(atom)

    @property
    def residue_type(self):
        return self._residue_type

    @residue_type.setter
    def residue_type(self, value):
        if isinstance(value, ResidueType):
            self._residue_type = value
        elif isinstance(value, str):
            try:
                self._residue_type = get_residue_type(value)
            except KeyError:
                residue_types[value] = value
                residue_string_types[value] = value
                self._residue_type = get_residue_type(value)

        else:
            raise ValueError("Invalid residue type.")

    @property
    def index(self):
        return self._index

    @index.setter
    def index(self, value):
        self._index = int(value)
