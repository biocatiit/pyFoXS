"""
This is the program for SAXS profile computation and fitting.
see FOXS for webserver (salilab.org/foxs)
"""

import os
import numpy as np
from scipy.spatial import distance

from .SolventAccessibleSurface import SolventAccessibleSurface
from .Profile import Profile
from .Particle import Particle
from .Residue import Residue
from .Atom import Atom, get_atom_type_exists, add_atom_type
from .Selector import NonAlternativePDBSelector, NonWaterNonHydrogenPDBSelector, NonHydrogenPDBSelector, NonWaterPDBSelector, CAlphaPDBSelector, get_default_pdb_selector

def compute_profile(particles, min_q, max_q, delta_q, ft, ff_type, hydration_layer, fit, reciprocal, ab_initio, vacuum, beam_profile_file):
    profile = Profile(qmin=min_q, qmax=max_q, delta=delta_q, constructor=0)
    # profile = Profile(min_q, max_q, delta_q)
    if reciprocal:
        profile.ff_table_ = ft
    if len(beam_profile_file) > 0:
        profile.beam_profile_ = beam_profile_file

    surface_area = []
    s = SolventAccessibleSurface()
    average_radius = 0.0
    if hydration_layer:
        for particle in particles:
            radius = ft.get_radius(particle, ff_type)
            particle.radius = radius
            average_radius += radius
        surface_area = s.get_solvent_accessibility(particles)
        average_radius /= len(particles)
        profile.average_radius_ = average_radius

    if not fit:
        if ab_initio:
            profile.calculate_profile_constant_form_factor(particles, ft)
        elif vacuum:
            profile.calculate_profile_partial(particles, surface_area, ff_type)
            profile.sum_partial_profiles(0.0, 0.0)
        else:
            profile.calculate_profile(particles, ff_type, reciprocal)
    else:
        if reciprocal:
            profile.calculate_profile_reciprocal_partial(particles, surface_area, ff_type)
        else:
            # default use
            profile.calculate_profile_partial(particles, surface_area, ff_type)

    return profile

def read_pdb(model, file, pdb_file_names, particles_vec, residue_level, heavy_atoms_only, multi_model_pdb, explicit_water):
    mhds = []
    selector = None
    if residue_level:
        selector = CAlphaPDBSelector()
    elif heavy_atoms_only:
        if explicit_water:
            selector = NonHydrogenPDBSelector()
        else:
            # default path
            selector = NonWaterNonHydrogenPDBSelector()
    else:
        if explicit_water:
            selector = NonAlternativePDBSelector()
        else:
            selector = NonWaterPDBSelector()

    if multi_model_pdb == 2:
        mhds = read_multimodel_pdb_or_mmcif(file, model, selector, True)
    else:
        if multi_model_pdb == 3:
            mhd = read_pdb_or_mmcif(file, model, selector, False, True)
            mhds.append(mhd)
        else:
            # Default path
            mhd = read_pdb_or_mmcif(file, model, selector, True, True)
            mhds.append(mhd)

    for h_index in range(len(mhds)):
        # ps = get_by_type(mhds[h_index], Atom)
        ps = mhds[h_index]
        if len(ps) > 0:
            pdb_id = file
            if len(mhds) > 1:
                pdb_id = trim_extension(file) + "_m" + str(h_index + 1) + ".pdb"
            pdb_file_names.append(pdb_id)
            particles_vec.append(ps)
            if len(mhds) > 1:
                print(str(len(ps)) + " atoms were read from PDB file " + file + " MODEL " + str(h_index + 1))
            else:
                print(str(len(ps)) + " atoms were read from PDB file " + file)
    return pdb_file_names, particles_vec

def get_by_type(root_particle, particle_type):
    particles = []
    visited = set()

    def depth_first_search(particle):
        visited.add(particle)
        if isinstance(particle, particle_type):
            particles.append(particle)

        if isinstance(particle, Particle):
            for child in particle.children:
                if child not in visited:
                    depth_first_search(child)

    depth_first_search(root_particle)
    return particles

def read_multimodel_pdb_or_mmcif(input, model, selector=get_default_pdb_selector(), noradii=False):
    if input.endswith('.cif'):
        return read_multimodel_mmcif(input, model, selector, noradii)
    else:
        return read_multimodel_pdb(input, model, selector, noradii)

def read_multimodel_mmcif(in_text, model, selector, noradii):
    fname = os.path.split(in_text)[1]
    nicename = os.path.splitext(fname)[0]
    ret = read_mmcif_hierarchy(in_text, nicename, fname,
                     model, selector, True, True, noradii)
    if len(ret) == 0:
        raise ValueError("No molecule read from file " + in_text)
    return ret

def read_mmcif_hierarchy(in_text, model, selector, select_first_model, noradii):
    fname = os.path.split(in_text)[1]
    nicename = os.path.splitext(fname)[0]
    ret = read_mmcif(in_text, nicename, fname,
                     model, selector, False, select_first_model, noradii)
    if len(ret) == 0:
        raise ValueError("No molecule read from file " + in_text)
    return ret[0]

def read_mmcif(in_stream, name, filename, model, selector, read_all_models,
               honor_model_num, noradii):
    print("read_mmcif not yet implemented")
    return
    err = None
    fh = ihm_file_new(read_callback, in_stream, None)
    r = ihm_reader_new(fh)
    ret = Hierarchies()

    asc = AtomSiteCategory(r, name, filename, model, ret, selector,
                           read_all_models, honor_model_num)

    more_data = 1
    if not ihm_read_file(r, more_data, err):
        errmsg = err.msg.decode('utf-8')
        ihm_error_free(err)
        ihm_reader_free(r)
        raise IOError(errmsg)
    ihm_reader_free(r)

    if not noradii:
        internal.add_pdb_radii(ret)

    return ret

def read_multimodel_pdb(in_stream, model, selector, noradii):
    fname = os.path.split(in_text)[1]
    nicename = os.path.splitext(fname)[0]
    ret = read_pdb(in_stream, nicename, fname,
                   model, selector, False, True, noradii)
    if len(ret) == 0:
        raise ValueError("No molecule read from file " + in_stream)
    return ret


def read_pdb_hierarchy(in_text, model, selector, select_first_model, no_radii):
    fname = os.path.split(in_text)[1]
    nicename = os.path.splitext(fname)[0]
    ret = read_pdb_func(in_text, nicename, fname, model, selector, select_first_model, False, no_radii)
    if len(ret) == 0:
        raise ValueError("No molecule read from file " + in_text)
    return ret[0]

# class Atom:
#     def __init__(self, atom_name, atom_type, coordinates):
#         self.atom_name = atom_name
#         self.atom_type = atom_type
#         self.coordinates = coordinates


# class Residue:
#     def __init__(self, residue_id, residue_name):
#         self.residue_id = residue_id
#         self.residue_name = residue_name
#         self.atoms = []


class Chain:
    def __init__(self, chain_id):
        self.name = "Chain"
        self.chain_id = chain_id
        self.parent = None
        self.residues = []


class Molecule:
    def __init__(self, molecule_id):
        self.molecule_id = molecule_id
        self.parent = None
        self.chains = []


# class Hierarchy:
#     def __init__(self):
#         self.molecules = []
#         self.parent = None
#     
#     def add_molecule(self, molecule):
#         self.molecules.append(molecule)


# class Particle:
#     def __init__(self, name):
#         self.name = name
#         self.children = []
# 
#     def add_child(self, child):
#         self.children.append(child)
#         child.parent = self

"""
def read_pdb_func(in_stream, name, filename, model, selector, select_first_model, split_models, noradii):
    hierarchy = Hierarchy()
    ret = []
    current_molecule = None
    current_chain = None
    current_residue = None
    current_particle = None
    root_particle = None
    has_atom = False

    f = open(in_stream, 'r')
    line = f.readline()
    while line:
        if line.startswith(('ATOM', 'HETATM')):
            if not selector.get_is_selected(line):
                print("Selector rejected line " + line)
                line = f.readline()
                continue
            atom_name = line[12:16].strip()
            atom_type = line[76:78].strip()
            coordinates = (float(line[30:38]), float(line[38:46]), float(line[46:54]))


            atom = Atom(atom_name, atom_type, coordinates)
            if current_residue is not None:
                current_residue.atoms.append(atom)
            
            has_atom = True

        elif line.startswith('TER'):
            current_residue = None

        elif line.startswith('ENDMDL'):
            current_chain = None
            if split_models:
                current_molecule = None

        elif line.startswith('MODEL'):
            model_number = int(line[10:14])
            if model_number == model or (select_first_model and model_number == 1):
                current_molecule = Molecule(model_number)
                hierarchy.add_molecule(current_molecule)

        elif line.startswith('CHAIN'):
            chain_id = line[21]
            if not selector or chain_id in selector:
                current_chain = Chain(chain_id)
                current_molecule.chains.append(current_chain)

        if line.startswith(('ATOM', 'HETATM')):
            residue_id = int(line[22:26])
            residue_name = line[17:20].strip()
            if current_chain is not None:
                if current_residue is None or residue_id != current_residue.residue_id:
                    current_residue = Residue(residue_id, residue_name)
                    current_chain.residues.append(current_residue)

            if root_particle is None:
                root_particle = Particle(residue_name)
                ret.append(root_particle)
            root_particle.add_child(atom)


        line = f.readline()

    if not has_atom:
        print("No atoms were read from " + filename +
                 "; perhaps it is not a PDB file.\n")
        return []

    return ret

def read_pdb_func(in_stream, name, filename, model, selector, select_first_model,
                split_models, noradii):
    ret = []
    root_name = ""
    root_p = None
    cp = None
    rp = None

    curr_residue_icode = '-'
    curr_chain = '-'
    chain_name_set = False
    first_model_read = False
    has_atom = False

    f = open(in_stream, 'r')
    line = f.readline()
    while line:
        line = line.strip()
        if line == "":
            line = f.readline()
            continue

        # handle MODEL reading
        if line[:5] == "MODEL":
            if first_model_read and select_first_model:
                break
            if split_models:
                root_name = line
                root_p = None
                cp = None
                rp = None
            first_model_read = True

        # check if line is an HETATM or ATOM record and if selector accepts the line
        # if true, create a new Particle using the line and add it to the Model
        if line[:4] == "ATOM" or line[:6] == "HETATM":
            if not selector.get_is_selected(line):
                print("Selector rejected line " + line)
                line = f.readline()
                continue

            residue_index = int(line[22:26])
            residue_icode = line[17:20].strip()
            chain = line[21]

            # create atom particle
            ap = atom_particle(model, line)
            # make sure all children have coordinates (no residues without valid atoms)
            if ap:
                # check if new chain
                if root_p is None:
                    root_p = Particle(model)
                    ret.append(setup_particle(model, root_p))
                    if root_name != "" or name != "":
                        root_p.name = name + ": " + root_name

                if cp is None or chain != curr_chain:
                    curr_chain = chain
                    # create new chain particle
                    cp = chain_particle(model, chain, filename)
                    chain_name_set = False
                    root_p.add_child(cp.chain)
                    rp = None  # make sure we get a new residue

                # check if new residue
                if rp is None or residue_index != rp.id or residue_icode != curr_residue_icode:
                    curr_residue_icode = residue_icode
                    # create new residue particle
                    rp = residue_particle(model, line)
                    cp.chain.add_child(rp.residue)

                # set chain name (protein/nucleotide/other) according to residue name
                if not chain_name_set:
                    cp.chain.name = rp.residue.name
                    chain_name_set = True

                rp.residue.add_child(ap.atom)
                has_atom = True

        line = f.readline()

    if not has_atom:
        print("No atoms were read from " + filename +
                 "; perhaps it is not a PDB file.\n")
        return []

    # TODO
    # if not noradii:
    #     internal.add_pdb_radii(ret)

    return ret

def read_pdb_func(in_stream, name, filename, model, selector, select_first_model,
                split_models, noradii):
    ret = []
    current_chain = None
    current_residue = None
    root_particle = None
    has_atom = False
    root_name = ""

    curr_residue_icode = '-'
    curr_chain = '-'
    chain_name_set = False
    first_model_read = False

    f = open(in_stream, 'r')
    line = f.readline()
    while line:
        line = line.strip()
        if line == "":
            line = f.readline()
            continue

        # handle MODEL reading
        if line[:5] == "MODEL":
            if first_model_read and select_first_model:
                break
            if split_models:
                root_name = line
                root_particle = None
                current_chain = None
                current_residue = None
            first_model_read = True

        # check if line is an HETATM or ATOM record and if selector accepts the line
        # if true, create a new Particle using the line and add it to the Model
        if line[:4] == "ATOM" or line[:6] == "HETATM":
            if not selector.get_is_selected(line):
                print("Selector rejected line " + line)
                line = f.readline()
                continue

            atom_name = line[12:16].strip()
            atom_type = line[76:78].strip()
            coordinates = (float(line[30:38]), float(line[38:46]), float(line[46:54]))

            atom = Atom(atom_name, atom_type, coordinates)
            
            residue_index = int(line[22:26])
            residue_icode = line[17:20].strip()
            chain = line[21]

            # # create atom particle
            # make sure all children have coordinates (no residues without valid atoms)
            if atom:
                # check if new chain
                if root_particle is None:
                    root_particle = Particle(model)
                    ret.append(setup_particle(model, root_particle))
                    if root_name != "" or name != "":
                        root_particle.name = name + ": " + root_name

                if current_chain is None or chain != curr_chain:
                    curr_chain = chain
                    current_chain = Chain(chain)
                    # current_molecule.chains.append(current_chain)
                    # create new chain particle
                    chain_name_set = False
                    root_particle.add_child(current_chain)
                    current_residue = None  # make sure we get a new residue

                # check if new residue
                if current_residue is None or residue_index != current_residue.index or residue_icode != curr_residue_icode:
                    curr_residue_icode = residue_icode
                    # create new residue particle
                    current_residue = Residue(residue_icode, residue_index)
                    current_chain.residues.append(current_residue)

                # set chain name (protein/nucleotide/other) according to residue name
                if not chain_name_set:
                    current_chain.name = current_residue.name
                    chain_name_set = True

                current_residue.add_child(atom)
                has_atom = True

        line = f.readline()

    if not has_atom:
        print("No atoms were read from " + filename +
                 "; perhaps it is not a PDB file.\n")
        return []

    # TODO
    # if not noradii:
    #     internal.add_pdb_radii(ret)
    return ret
    
"""

def read_pdb_func(in_stream, name, filename, model, selector, select_first_model,
                split_models, noradii):
    ret = []
    current_chain = None
    current_residue = None
    root_particle = None
    has_atom = False
    root_name = ""

    curr_residue_icode = '-'
    curr_chain = '-'
    chain_name_set = False
    first_model_read = False

    f = open(in_stream, 'r')
    line = f.readline()
    while line:
        line = line.strip()
        if line == "":
            line = f.readline()
            continue

        # handle MODEL reading
        if line[:5] == "MODEL":
            if first_model_read and select_first_model:
                break
            if split_models:
                root_name = line
                root_particle = None
                current_chain = None
                current_residue = None
            first_model_read = True

        # check if line is an HETATM or ATOM record and if selector accepts the line
        # if true, create a new Particle using the line and add it to the Model
        if line[:4] == "ATOM" or line[:6] == "HETATM":
            if not selector.get_is_selected(line):
                print("Selector rejected line " + line)
                line = f.readline()
                continue

            string_name = atom_name = line[12:16].strip()
            atom_type = line[76:78].strip()
            coordinates = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            residue_index = int(line[22:26])
            residue_icode = line[17:20].strip()
            chain = line[21]
            is_hetatm = line[:6] == "HETATM"
            residue_index = int(line[22:26])
            occupancy = float(line[54:60])
            temp_factor = float(line[60:66])
            # Determine AtomType
            if is_hetatm:
                string_name = "HET:" + string_name
                if not get_atom_type_exists(string_name):
                    atom_name = add_atom_type(string_name, atom_type)
                else:
                    atom_name = string_name
            else:  # ATOM line
                string_name = string_name.strip()
                if string_name == "":
                    string_name = "UNK"
                if not get_atom_type_exists(string_name):
                    print("ATOM record type not found: \""
                                    + string_name + "\" in PDB file " + "\n")
                    atom_name = add_atom_type(string_name, atom_type)
                else:
                    atom_name = string_name
            
            atom = Atom(name=atom_name, type=atom_type, coord=coordinates, occupancy=occupancy, temp_factor=temp_factor)

            # # create atom particle
            # make sure all children have coordinates (no residues without valid atoms)
            if atom:
                # check if new chain
                if root_particle is None:
                    root_particle = Particle(model)
                    
                    if root_name != "" or name != "":
                        root_particle.name = name + ": " + root_name

                if current_chain is None or chain != curr_chain:
                    curr_chain = chain
                    current_chain = Chain(chain)
                    # current_molecule.chains.append(current_chain)
                    # create new chain particle
                    chain_name_set = False
                    root_particle.add_child(current_chain)
                    current_residue = None  # make sure we get a new residue

                # check if new residue
                if current_residue is None or residue_index != current_residue.index or residue_icode != curr_residue_icode:
                    curr_residue_icode = residue_icode
                    # create new residue particle
                    current_residue = Residue(residue_icode, residue_index)
                atom.residue = current_residue

                # set chain name (protein/nucleotide/other) according to residue name
                if not chain_name_set:
                    current_chain.name = current_residue.name
                    chain_name_set = True

                current_residue.add_child(atom)
                ret.append(atom)
                has_atom = True

        line = f.readline()

    if not has_atom:
        print("No atoms were read from " + filename +
                 "; perhaps it is not a PDB file.\n")
        return []

    return [ret]

def read_pdb_or_mmcif(input, model, selector=get_default_pdb_selector(), select_first_model=True, noradii=False):
    if input.endswith(".cif"):
        return read_mmcif(input, model, selector, select_first_model, noradii)
    else:
        return read_pdb_hierarchy(input, model, selector, select_first_model, noradii)

# import ihm.reader
# from .Atom import AtomSiteCategory, add_pdb_radii
# from .Hierarchy import Hierarchy
# 
# def read_mmcif(in_stream, name, filename, model, selector,
#                read_all_models, honor_model_num, noradii):
#     ret = Hierarchy()
# 
#     def read_callback(buffer_size, user_data):
#         return in_stream.read(buffer_size)
# 
#     err = None
#     fh = ihm.reader.file_new(read_callback, in_stream, None)
#     r = ihm.reader.new(fh)
# 
#     try:
#         asc = AtomSiteCategory(r, name, filename, model, ret, selector,
#                                read_all_models, honor_model_num)
#         more_data = True
#         while more_data:
#             more_data, err = ihm.reader.read_file(r)
#         if err:
#             errmsg = err.msg.decode()
#             raise errmsg
#     finally:
#         ihm.reader.free(r)
# 
#     if not noradii:
#         add_pdb_radii(ret)
# 
#     return ret


def read_files(m, files, pdb_file_names, dat_files, particles_vec, exp_profiles, residue_level, heavy_atoms_only, multi_model_pdb, explicit_water, max_q, units):
    for file in files:
        # check if file exists
        if not os.path.exists(file):
            print("Can't open file " + file)
            return
        # 1. try as pdb or mmcif
        try:
            read_pdb(m, file, pdb_file_names, particles_vec, residue_level, heavy_atoms_only, multi_model_pdb, explicit_water)
        except ValueError:  # not a pdb file
            # 2. try as a dat profile file
            profile = Profile(file_name=file, fit_file=False, max_q=max_q, units=units, constructor=1)
            # profile = Profile(file, False, max_q, units)
            if profile.size() == 0:
                print("Can't parse input file " + file)
                return
            else:
                dat_files.append(file)
                exp_profiles.append(profile)
                print("Profile read from file " + file + " size = " + str(profile.size()))

def trim_extension(file_name):
    # if file_name[-4:] == '.':
    #     return file_name[:-4]
    return os.path.splitext(file_name)[0]

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
