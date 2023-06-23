
def get_atom(rd, at):
    mhd = Hierarchy(rd.get_particle())
    for i in range(mhd.get_number_of_children()):
        a = Atom(mhd.get_child(i))
        if a.get_atom_type() == at:
            return a
    print("Atom not found " + str(at) + "\n")
    return Atom()

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
