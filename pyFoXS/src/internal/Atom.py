from enum import IntEnum
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

def get_atom(rd, at):
    mhd = Hierarchy(rd.get_particle())
    for i in range(mhd.get_number_of_children()):
        a = Atom(mhd.get_child(i))
        if a.get_atom_type() == at:
            return a
    print("Atom not found " + str(at) + "\n")
    return Atom()

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
    def do_setup_particle(m, pi, t):
        m.add_attribute(Atom.get_atom_type_key(), pi, t.get_index())
        if not Hierarchy.get_is_setup(m, pi):
            Hierarchy.setup_particle(m, pi)
        m.add_attribute(Atom.get_element_key(), pi, Element.UNKNOWN_ELEMENT)

        ret = Atom(m, pi)
        if not Mass.get_is_setup(m, pi):
            Mass.setup_particle(m, pi, 0)
        ret.set_atom_type(t)

    @staticmethod
    def do_setup_particle_from_atom(m, pi, o):
        Atom.do_setup_particle(m, pi, o.get_atom_type())

    # def __init__(self, m, pi):
    #     self.m = m
    #     self.pi = pi

    def show(self, out):
        if self.get_input_index() != -1:
            out.write("#" + str(self.get_input_index()) + " ")
        out.write(self.get_atom_type())
        out.write(" (" + get_element_table().get_name(self.get_element()) + ")")

    def set_atom_type(self, t):
        self.get_particle().set_value(Atom.get_atom_type_key(), t.get_index())
        e = get_element_for_atom_type(t)
        if e != Element.UNKNOWN_ELEMENT:
            self.set_element(e)
        self.get_model().set_trigger_updated(Residue.get_type_changed_key())

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

    def set_element(self, e):
        self.get_particle().set_value(Atom.get_element_key(), e)
        Mass(self.get_particle()).set_mass(get_element_table().get_mass(e))


def ihm_category_new(reader, name, data_callback, end_frame_callback,
                     finalize_callback, data, free_func):
    category = {
        'name': name,
        'data_callback': data_callback,
        'end_frame_callback': end_frame_callback,
        'finalize_callback': finalize_callback,
        'data': data,
        'free_func': free_func,
        'keyword_map': ihm_mapping_new(ihm_keyword_free)
    }
    reader['category_map'][category['name']] = category
    return category


class Category:
    def __init__(self, reader, name, callback):
        self.c_ = ihm_category_new(reader, name, callback, None, None, self, None)

class Keyword:
    def __init__(self, c, name):
        self.k_ = ihm_keyword_new(c, name.encode())

    def data(self):
        return self.k_.data

    def as_str(self):
        if self.k_.omitted or self.k_.unknown or not self.k_.in_file:
            return ""
        else:
            return self.k_.data

    def as_float(self, default_value=0.0):
        if self.k_.omitted or self.k_.unknown or not self.k_.in_file:
            return default_value
        else:
            return float(self.k_.data)

    def as_int(self, default_value=0):
        if self.k_.omitted or self.k_.unknown or not self.k_.in_file:
            return default_value
        else:
            return int(self.k_.data)

class AtomSiteCategory(Category):
    def __init__(self, reader, name, filename, model, hiers, selector,
                 read_all_models, honor_model_num):
        super().__init__(reader, "_atom_site", self.callback)
        self.name_ = name
        self.filename_ = filename
        self.model_ = model
        self.selector_ = selector
        self.read_all_models_ = read_all_models
        self.honor_model_num_ = honor_model_num
        self.atom_name_ = Keyword(c_, "label_atom_id")
        self.residue_name_ = Keyword(c_, "label_comp_id")
        self.chain_ = Keyword(c_, "label_asym_id")
        self.auth_chain_ = Keyword(c_, "auth_asym_id")
        self.element_ = Keyword(c_, "type_symbol")
        self.seq_id_ = Keyword(c_, "label_seq_id")
        self.group_ = Keyword(c_, "group_pdb")
        self.id_ = Keyword(c_, "id")
        self.occupancy_ = Keyword(c_, "occupancy")
        self.temp_factor_ = Keyword(c_, "b_iso_or_equiv")
        self.ins_code_ = Keyword(c_, "pdbx_pdb_ins_code")
        self.x_ = Keyword(c_, "cartn_x")
        self.y_ = Keyword(c_, "cartn_y")
        self.z_ = Keyword(c_, "cartn_z")
        self.model_num_ = Keyword(c_, "pdbx_pdb_model_num")
        self.auth_seq_id_ = Keyword(c_, "auth_seq_id")
        self.alt_loc_id_ = Keyword(c_, "label_alt_id")
        self.cp_ = None
        self.rp_ = None
        self.root_p_ = None
        self.hiers_ = hiers
        self.curr_chain_ = ""
        self.curr_seq_id_ = 0
        self.curr_auth_seq_id_ = 0
        self.curr_model_num_ = 0
        self.curr_residue_icode_ = ""
        self.hetatm_ = "HETATM"
        self.chain_map_ = {}
        self.root_map_ = {}
        self.pdb_line_ = ""

    def callback(self, reader, data, error):
        self.handle()

    def set_root_particle_name(self, model_num):
        oss = StringIO()
        oss.write(str(model_num))
        root_name = oss.getvalue()
        root_p_.set_name(name_ + ": " + root_name)

    def get_root_particle(self, model_num):
        if self.root_p_ is None or model_num != self.curr_model_num_:
            if not self.read_all_models_ and self.root_p_ is not None:
                return False
            self.curr_model_num_ = model_num
            if model_num not in self.root_map_:
                root_p_ = Particle(model_)
                self.set_root_particle_name(model_num)
                self.hiers_.append(Hierarchy.setup_particle(root_p_))
                self.root_map_[model_num] = root_p_
            else:
                self.root_p_ = self.root_map_[model_num]
            self.cp_ = None
        return True

    def get_chain_particle(self, chain):
        if self.cp_ is None or chain != self.curr_chain_:
            self.curr_chain_ = chain
            root_chain = (root_p_, chain)
            if root_chain not in self.chain_map_:
                cp_ = internal.chain_particle(model_, chain, filename_)
                Hierarchy(root_p_).add_child(Chain(cp_))
                self.chain_map_[root_chain] = cp_
            else:
                self.cp_ = self.chain_map_[root_chain]
            self.rp_ = None

    def replace(self, dest, pos, maxlen, repl):
        len_ = min(maxlen, len(repl))
        if len_ > 0:
            dest[pos:pos + len_] = repl[:len_]

    def get_is_selected(self):
        pdb_line_ = " " * 80
        self.replace(pdb_line_, atom_entry_type_field_, 6,
                     self.group_.as_str())
        if (len(self.atom_name_.as_str()) >= 4 or
                len(self.element_.as_str()) == 2):
            self.replace(pdb_line_, atom_type_field_, 4,
                         self.atom_name_.as_str())
        else:
            self.replace(pdb_line_, atom_type_field_ + 1, 3,
                         self.atom_name_.as_str())
        self.replace(pdb_line_, atom_alt_loc_field_, 1,
                     self.alt_loc_id_.as_str())
        self.replace(pdb_line_, atom_res_name_field_, 3,
                     self.residue_name_.as_str())
        self.replace(pdb_line_, atom_chain_id_field_, 1,
                     self.chain_.as_str())
        start = self.auth_seq_id_.as_str()
        endptr = ctypes.POINTER(ctypes.c_char)()
        int(start, ctypes.byref(endptr), 10)
        pdb_line_[atom_res_insertion_field_] = endptr.contents or ' '
        self.replace(pdb_line_, atom_element_field_, 2,
                     self.element_.as_str())
        return self.selector_.get_is_selected(pdb_line_)

    def handle(self):
        if not self.get_is_selected():
            return
        if not self.get_root_particle(self.curr_model_num_ if self.honor_model_num_ else 1):
            return

        e = get_element_table().get_element(self.element_.as_str())
        seq_id = self.seq_id_.as_int(1)
        residue_icode = self.ins_code_.as_str()

        if len(self.auth_chain_.as_str()) > 0:
            self.get_chain_particle(self.auth_chain_.as_str())
        else:
            self.get_chain_particle(self.chain_.as_str())

        if (self.rp_ is None or seq_id != self.curr_seq_id_ or
                residue_icode != self.curr_residue_icode_):
            self.curr_seq_id_ = seq_id
            self.curr_residue_icode_ = residue_icode
            si = self.auth_seq_id_.as_str()
            start = si
            endptr = ctypes.POINTER(ctypes.c_char)()
            auth_seq_id = int(start, ctypes.byref(endptr), 10)
            if not endptr.contents:
                auth_seq_id = seq_id
            one_icode = 32
            if endptr.contents:
                one_icode = endptr.contents
            self.curr_auth_seq_id_ = auth_seq_id
            self.rp_ = internal.residue_particle(
                model_, auth_seq_id, one_icode, self.residue_name_.as_str()
            )
            Chain(cp_).add_child(Residue(rp_))

        ap = internal.atom_particle(
            model_, self.atom_name_.as_str(), e,
            self.group_.as_str() == "HETATM", self.id_.as_int(),
            self.curr_auth_seq_id_, self.x_.as_float(), self.y_.as_float(),
            self.z_.as_float(), self.occupancy_.as_float(),
            self.temp_factor_.as_float()
        )
        Residue(rp_).add_child(Atom(ap))

class Chain:
    def __init__(self, chain_id):
        self.name = "Chain" + str(chain_id)
        self.chain_id = chain_id
        self.residues = []

    def add_child(self, residue):
        self.residues.append(residue)

class Molecule:
    def __init__(self, molecule_id):
        self.name = "Molecule" + str(molecule_id)
        self.molecule_id = molecule_id
        self.chains = []

    def add_chain(self, chain):
        self.chains.append(chain)