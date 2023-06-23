from .Particle import Particle

class HierarchyTraits:
    def __init__(self):
        self.children_ = None
        self.parent_ = None

    def __init__(self, name):
        self.children_ = None
        self.parent_ = None
        self.name = name

    def get_children_key(self):
        return self.children_

    def get_parent_key(self):
        return self.parent_
    
    def get_traits(self):
        return self

    def __eq__(self, o):
        return self.parent_ == o.parent_

    def __repr__(self):
        return f"HierarchyTraits(parent={self.parent_})"

"""
class Hierarchy:
    def __init__(self, model, particle_index):
        super().__init__(model, particle_index, self.get_traits())

    @staticmethod
    def setup_particle(model, particle_index,
                       children=[]):
        H.setup_particle(model, particle_index, Hierarchy.get_traits())
        ret = Hierarchy(model, particle_index)
        for child_index in children:
            if not Hierarchy.get_is_setup(model, child_index):
                Hierarchy.setup_particle(model, child_index)
            ret.add_child(Hierarchy(model, child_index))
        return ret

    @staticmethod
    def get_is_setup(model, particle_index):
        return H.get_is_setup(model, particle_index, Hierarchy.get_traits())

    def get_is_valid(self, print_info: bool = False):
        return super().get_is_valid(print_info)

    def add_child(self, child):
        if child != self:
            super().add_child(child)

    def get_child(self, i: int):
        return Hierarchy(super().get_child(i))

    def get_children(self):
        return [Hierarchy(child) for child in super().get_children()]

    def get_parent(self):
        parent = super().get_parent()
        if parent == Hierarchy():
            return Hierarchy()
        else:
            return Hierarchy(parent)
"""
class Hierarchy:
    def __init__(self):
        self.atoms = []
        self.residues = {}
        self.chains = {}
        self.models = {}

    def add_atom(self, atom):
        self.atoms.append(atom)

    def add_residue(self, residue):
        residue_id = residue.get_id()
        chain_id = residue.get_parent().get_id()
        model_id = residue.get_parent().get_parent().get_id()

        if chain_id not in self.chains:
            self.chains[chain_id] = {}
        if residue_id not in self.chains[chain_id]:
            self.chains[chain_id][residue_id] = residue

        if residue_id not in self.residues:
            self.residues[residue_id] = residue

        if chain_id not in self.models:
            self.models[chain_id] = model_id

    def get_atom_count(self):
        return len(self.atoms)

    def get_residue_count(self):
        return len(self.residues)

    def get_chain_count(self):
        return len(self.chains)

    def get_model_count(self):
        return len(self.models)
    
def setup_particle(model, particle_index, children=[]):
    particle = Particle(model)
    # particle.id = particle_index

    for child in children:
        Particle(model, child).add_child(particle_index)

    return particle


def get_by_type(mhd, t):
    type_names = [
        "ATOM_TYPE", "RESIDUE_TYPE", "CHAIN_TYPE", "MOLECULE_TYPE",
        "DOMAIN_TYPE", "FRAGMENT_TYPE", "XYZ_TYPE", "XYZR_TYPE",
        "MASS_TYPE", "STATE_TYPE"
    ]
    type_name = type_names[t]
    return [Hierarchy(p) for p in Hierarchy.get_by_type(mhd, type_name)]

def get_residue(mhd, index):
    return Hierarchy(Hierarchy.get_residue(mhd, index))

def create_fragment(particles):
    fragment = Hierarchy.create_fragment([Hierarchy(p) for p in particles])
    return Hierarchy(fragment)

def get_internal_bonds(mhd):
    return [Bond(b) for b in Hierarchy.get_internal_bonds(mhd)]

def get_root(h):
    while h.get_parent():
        h = h.get_parent()
    return h

def get_leaves(h):
    return [Hierarchy(leaf) for leaf in IMP.core.get_leaves(h)]
