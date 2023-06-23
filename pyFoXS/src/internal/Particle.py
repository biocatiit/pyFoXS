from .Atom import get_element_for_atom_type, get_atom_type_exists, add_atom_type, Atom, Chain
from .Residue import Residue

StructureProvenance = {}

class Particle:
    def __init__(self, model, name="", id=0, element=None, elname="", coordinates=(0,0,0), radius=None, atom=None, chain=None, molecule=None, residue=None):
        self.name = name
        self.id = id
        self.element = element
        self.elname = elname
        self.coordinates = coordinates
        self.model = model
        self.was_used = False
        self.radius = radius
        self.atom = atom
        self.chain = chain
        self.molecule = molecule
        self.residue = residue
        self.parent = None
        self.children = []

    def setup_particle(self, radius):
        self.radius = radius
    
    def add_child(self, child):
        self.children.append(child)
        child.parent = self

def atom_particle(model, line):
    atom_type = line[12:16].strip()
    element = line[76:78].strip()
    is_hetatm = line[:6] == "HETATM"
    atom_index = int(line[6:11])
    residue_index = int(line[22:26])
    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
    occupancy = float(line[54:60])
    temp_factor = float(line[60:66])
    atom_name = None
    string_name = atom_type

    # Determine AtomType
    if is_hetatm:
        string_name = "HET:" + string_name
        if not get_atom_type_exists(string_name):
            atom_name = add_atom_type(string_name, element)
        else:
            atom_name = string_name
    else:  # ATOM line
        string_name = string_name.strip()
        if string_name == "":
            string_name = "UNK"
        if not get_atom_type_exists(string_name):
            print("ATOM record type not found: \""
                            + string_name + "\" in PDB file " + "\n")
            atom_name = add_atom_type(string_name, element)
        else:
            atom_name = string_name

    # Create new particle
    p = Particle(model)
    # p.add_attribute(get_pdb_index_key(), atom_index)
    v = (x, y, z)

    # Atom decorator
    p.atom = Atom(name=atom_name, type=atom_type, element=element, atom_index=atom_index, coord=v, occupancy=occupancy, temp_factor=temp_factor)
    oss = "Atom " + atom_name + " of residue " + str(residue_index)
    p.name = oss

    # Check if the element matches
    e2 = get_element_for_atom_type(atom_name)
    if element != e2:
        print("AtomType element and PDB line elements don't match. AtomType "
                        + str(e2) + " vs. determined from PDB " + str(element) + "\n")

    return p

def chain_particle(m, chain_id, filename):
    p = Particle(m)
    p.chain = Chain(chain_id)
    p.name = "Chain " + str(chain_id)
    # atom.Molecule.setup_particle(p)

    # Set provenance of this chain
    sp = StructureProvenance[Particle(m)] = (filename, chain_id)
    # atom.core.add_provenance(m, p.get_index(), sp)

    return p

def residue_particle(model, line):
    residue_name = line[17:20].strip()
    residue_index = int(line[22:26])
    ins_code = line[26]

    p = Particle(model)

    if not residue_name:
        residue_name = "UNK"

    p.residue = Residue(residue_name, residue_index, ord(ins_code))
    p.name = residue_name

    return p


"""
class StructureProvenance(core.Provenance):
    @staticmethod
    def setup_particle(m, pi,
                          filename: str, chain_id: str,
                          residue_offset: int = 0) -> None:
        core.Provenance.setup_particle(m, pi)
        assert filename, "The filename cannot be empty."
        m.add_attribute(StructureProvenance.get_filename_key(), pi,
                        filename)
        m.add_attribute(StructureProvenance.get_chain_key(), pi, chain_id)
        m.add_attribute(StructureProvenance.get_residue_offset_key(), pi,
                        residue_offset)

    @staticmethod
    def do_setup_particle(m, pi, o):
        StructureProvenance.do_setup_particle(m, pi, o.get_filename(),
                                              o.get_chain_id(),
                                              o.get_residue_offset())

    @staticmethod
    def get_filename_key() -> str:
        return "filename"

    @staticmethod
    def get_chain_key() -> str:
        return "chain_id"

    @staticmethod
    def get_residue_offset_key() -> str:
        return "residue_offset"

    @classmethod
    def get_is_setup(cls, m, pi) -> bool:
        return (m.get_has_attribute(cls.get_filename_key(), pi) and
                m.get_has_attribute(cls.get_chain_key(), pi) and
                m.get_has_attribute(cls.get_residue_offset_key(), pi))

    def set_filename(self, filename: str) -> None:
        assert filename, "The filename cannot be empty"
        self.get_model().set_attribute(self.get_filename_key(),
                                        self.get_particle_index(),
                                        core.get_absolute_path(filename))

    def get_filename(self) -> str:
        return self.get_model().get_attribute(self.get_filename_key(),
                                              self.get_particle_index())

    def set_chain_id(self, chain_id: str) -> None:
        self.get_model().set_attribute(self.get_chain_key(),
                                        self.get_particle_index(),
                                        chain_id)

    def get_chain_id(self) -> str:
        return self.get_model().get_attribute(self.get_chain_key(),
                                              self.get_particle_index())

    def set_residue_offset(self, residue_offset: int) -> None:
        self.get_model().set_attribute(self.get_residue_offset_key(),
                                        self.get_particle_index(),
                                        residue_offset)

    def get_residue_offset(self) -> int:
        return self.get_model().get_attribute(self.get_residue_offset_key(),
                                              self.get_particle_index())
"""