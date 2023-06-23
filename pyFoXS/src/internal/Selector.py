import string
from abc import abstractmethod

class PDBSelector:
    def __init__(self, name="PDBSelector%1%"):
        self.name = name

    @abstractmethod
    def get_is_selected(self, pdb_line):
        pass

class NonAlternativePDBSelector(PDBSelector):
    def __init__(self, name="NonAlternativePDBSelector%1%"):
        PDBSelector.__init__(self, name)

    def get_is_selected(self, pdb_line):
        alt_loc_indicator = pdb_line[16]
        return alt_loc_indicator == ' ' or alt_loc_indicator == 'A'

class WaterPDBSelector(NonAlternativePDBSelector):
    def __init__(self, name="WaterPDBSelector%1%"):
        NonAlternativePDBSelector.__init__(self, name)

    def get_is_selected(self, pdb_line):
        if not NonAlternativePDBSelector.get_is_selected(self, pdb_line):
            return False
        res_name = pdb_line[17:20].strip()
        return (res_name[:3] == 'HOH' or res_name[:3] == 'DOD')

class HydrogenPDBSelector(NonAlternativePDBSelector):
    def __init__(self, name="HydrogenPDBSelector%1%"):
        NonAlternativePDBSelector.__init__(self, name)

    def is_hydrogen(self, pdb_line):
        if not NonAlternativePDBSelector.get_is_selected(self, pdb_line):
            return False
        elem = pdb_line[76:78].strip()
        # Determine if the line is a hydrogen atom
        if len(elem) == 1 and elem == 'H':
            return True
        if len(elem) == 2 and elem[0] == 'H' and elem[1] in ['E', 'e', 'O', 'o', 'F', 'f', 'G', 'g']:
            return False
        atom_name = pdb_line[12:16].strip()
        if (atom_name[0] == ' ' or atom_name[0] in string.digits) and (atom_name[1] == 'H' or atom_name[1] == 'D'):
            return True
        if atom_name[0] == 'H' or atom_name[0] == 'D':
            return True
        return False

    def get_is_selected(self, pdb_line):
        if not NonAlternativePDBSelector.get_is_selected(self, pdb_line):
            return False
        return self.is_hydrogen(pdb_line)

class NonWaterPDBSelector(NonAlternativePDBSelector):
    def __init__(self, name="NonWaterPDBSelector%1%"):
        super().__init__(name)
        self.ws_ = WaterPDBSelector()

    def get_is_selected(self, pdb_line):
        if not super().get_is_selected(pdb_line):
            return False
        return not self.ws_.get_is_selected(pdb_line)

class NonHydrogenPDBSelector(NonAlternativePDBSelector):
    def __init__(self, name="NonHydrogenPDBSelector%1%"):
        super().__init__(name)
        self.hs_ = HydrogenPDBSelector()

    def get_is_selected(self, pdb_line):
        if not super().get_is_selected(pdb_line):
            return False
        return not self.hs_.get_is_selected(pdb_line)

class NonWaterNonHydrogenPDBSelector(NonAlternativePDBSelector):
    def __init__(self, name="NonWaterNonHydrogenPDBSelector"):
        NonAlternativePDBSelector.__init__(self, name)
        self.ws_ = WaterPDBSelector()
        self.hs_ = HydrogenPDBSelector()

    def get_is_selected(self, pdb_line):
        if not NonAlternativePDBSelector.get_is_selected(self, pdb_line):
            return False
        return not (self.ws_.get_is_selected(pdb_line) or
                    self.hs_.get_is_selected(pdb_line))

class CAlphaPDBSelector(NonAlternativePDBSelector):
    def __init__(self, name="CAlphaPDBSelector%1%"):
        super().__init__(name)

    def get_is_selected(self, pdb_line):
        if not super().get_is_selected(pdb_line):
            return False
        atom_type = pdb_line[12:16].strip()
        return atom_type == 'CA' # (atom_type == 'C' and atom_type[2] == 'A' and atom_type[3] == ' ')

def get_default_pdb_selector():
    return NonWaterPDBSelector()
