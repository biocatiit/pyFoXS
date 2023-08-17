"""
This is the program for SAXS profile computation and fitting.
see FOXS for webserver (salilab.org/foxs)
"""

import os

from .Profile import Profile
from .Selector import NonAlternativePDBSelector, NonWaterNonHydrogenPDBSelector, NonHydrogenPDBSelector, NonWaterPDBSelector, CAlphaPDBSelector, get_default_pdb_selector
from ..structure.SolventAccessibleSurface import SolventAccessibleSurface
from ..structure.Residue import Residue
from ..structure.Atom import Atom, Chain, get_atom_type_exists, add_atom_type

def compute_profile(particles, min_q, max_q, delta_q, ft, ff_type,
    hydration_layer, fit, reciprocal, ab_initio, vacuum, beam_profile_file):
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

def read_files(m, files, pdb_file_names, dat_files, particles_vec, exp_profiles,
    residue_level, heavy_atoms_only, multi_model_pdb, explicit_water, max_q, units):
    for file in files:
        # check if file exists
        if not os.path.exists(file):
            print("Can't open file " + file)
            return
        # 1. try as pdb or mmcif
        try:
            read_pdb(m, file, pdb_file_names, particles_vec, residue_level,
                heavy_atoms_only, multi_model_pdb, explicit_water)
        except ValueError:  # not a pdb file
            # 2. try as a dat profile file
            profile = Profile(file_name=file, fit_file=False, max_q=max_q,
                units=units, constructor=1)
            # profile = Profile(file, False, max_q, units)
            if profile.size() == 0:
                print("Can't parse input file " + file)
                return
            else:
                dat_files.append(file)
                exp_profiles.append(profile)
                print("Profile read from file " + file + " size = " + str(profile.size()))

def read_pdb(model, file, pdb_file_names, particles_vec, residue_level,
    heavy_atoms_only, multi_model_pdb, explicit_water):
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

def read_pdb_or_mmcif(in_text, model, selector=get_default_pdb_selector(),
    select_first_model=True, noradii=False):
    if in_text.endswith(".cif"):
        return read_mmcif(in_text, model, selector, select_first_model, noradii)
    else:
        return read_pdb_hierarchy(in_text, model, selector, select_first_model,
            noradii)

def read_multimodel_pdb_or_mmcif(in_text, model,
    selector=get_default_pdb_selector(), noradii=False):
    if in_text.endswith('.cif'):
        return read_multimodel_mmcif(in_text, model, selector, noradii)
    else:
        return read_multimodel_pdb(in_text, model, selector, noradii)

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
    return []
    # err = None
    # fh = ihm_file_new(read_callback, in_stream, None)
    # r = ihm_reader_new(fh)
    # ret = Hierarchies()

    # asc = AtomSiteCategory(r, name, filename, model, ret, selector,
    #                        read_all_models, honor_model_num)

    # more_data = 1
    # if not ihm_read_file(r, more_data, err):
    #     errmsg = err.msg.decode('utf-8')
    #     ihm_error_free(err)
    #     ihm_reader_free(r)
    #     raise IOError(errmsg)
    # ihm_reader_free(r)

    # if not noradii:
    #     internal.add_pdb_radii(ret)

    # return ret

def read_multimodel_pdb(in_stream, model, selector, noradii):
    fname = os.path.split(in_stream)[1]
    nicename = os.path.splitext(fname)[0]
    ret = read_pdb_func(in_stream, nicename, fname,
                   model, selector, False, True, noradii)
    if len(ret) == 0:
        raise ValueError("No molecule read from file " + in_stream)
    return ret

def read_pdb_hierarchy(in_text, model, selector, select_first_model, no_radii):
    fname = os.path.split(in_text)[1]
    nicename = os.path.splitext(fname)[0]
    ret = read_pdb_func(in_text, nicename, fname, model, selector,
        select_first_model, False, no_radii)
    if len(ret) == 0:
        raise ValueError("No molecule read from file " + in_text)
    return ret[0]

def read_pdb_func(in_stream, name, filename, model, selector, select_first_model,
                split_models, noradii):
    ret = []
    current_chain = None
    current_residue = None
    root_particle = None
    has_atom = False

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
                root_particle = None
                current_chain = None
                current_residue = None
            first_model_read = True

        # check if line is an HETATM or ATOM record and if selector accepts the line
        # if true, create a new Particle using the line and add it to the Model
        if line[:4] == "ATOM" or line[:6] == "HETATM":
            if not selector.get_is_selected(line):
                # print("Selector rejected line " + line)
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

            atom = Atom(name=atom_name, type=atom_type, coord=coordinates,
                occupancy=occupancy, temp_factor=temp_factor)

            # # create atom particle
            # make sure all children have coordinates (no residues without valid atoms)
            if atom:
                # check if new chain
                if root_particle is None:
                    root_particle = [] # Particle(model)

                if current_chain is None or chain != curr_chain:
                    curr_chain = chain
                    current_chain = Chain(chain)
                    # current_molecule.chains.append(current_chain)
                    # create new chain particle
                    chain_name_set = False
                    root_particle.append(current_chain)
                    current_residue = None  # make sure we get a new residue

                # check if new residue
                if (current_residue is None
                    or residue_index != current_residue.index
                    or residue_icode != curr_residue_icode):
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

def trim_extension(file_name):
    # if file_name[-4:] == '.':
    #     return file_name[:-4]
    return os.path.splitext(file_name)[0]
