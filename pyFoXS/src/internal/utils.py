"""
This is the program for SAXS profile computation and fitting.
see FOXS for webserver (salilab.org/foxs)
"""

import IMP
import IMP.saxs
import os

def compute_profile(particles, min_q, max_q, delta_q, ft, ff_type, hydration_layer, fit, reciprocal, ab_initio, vacuum, beam_profile_file):
    profile = IMP.saxs.Profile(min_q, max_q, delta_q)
    if reciprocal:
        profile.set_ff_table(ft)
    if len(beam_profile_file) > 0:
        profile.set_beam_profile(beam_profile_file)

    surface_area = []
    s = IMP.saxs.SolventAccessibleSurface()
    average_radius = 0.0
    if hydration_layer:
        for particle in particles:
            radius = ft.get_radius(particle, ff_type)
            IMP.core.XYZR.setup_particle(particle, radius)
            average_radius += radius
        surface_area = s.get_solvent_accessibility(IMP.core.XYZRs(particles))
        average_radius /= len(particles)
        profile.set_average_radius(average_radius)

    if not fit:
        if ab_initio:
            profile.calculate_profile_constant_form_factor(particles)
        elif vacuum:
            profile.calculate_profile_partial(particles, surface_area, ff_type)
            profile.sum_partial_profiles(0.0, 0.0)
        else:
            profile.calculate_profile(particles, ff_type, reciprocal)
    else:
        if reciprocal:
            profile.calculate_profile_reciprocal_partial(particles, surface_area, ff_type)
        else:
            profile.calculate_profile_partial(particles, surface_area, ff_type)

    return profile

def read_pdb(model, file, pdb_file_names, particles_vec, residue_level, heavy_atoms_only, multi_model_pdb, explicit_water):
    mhds = []
    selector = None
    if residue_level:
        selector = IMP.atom.CAlphaPDBSelector()
    elif heavy_atoms_only:
        if explicit_water:
            selector = IMP.atom.NonHydrogenPDBSelector()
        else:
            selector = IMP.atom.NonWaterNonHydrogenPDBSelector()
    else:
        if explicit_water:
            selector = IMP.atom.NonAlternativePDBSelector()
        else:
            selector = IMP.atom.NonWaterPDBSelector()

    if multi_model_pdb == 2:
        mhds = read_multimodel_pdb_or_mmcif(file, model, selector, True)
    else:
        if multi_model_pdb == 3:
            mhd = IMP.atom.read_pdb_or_mmcif(file, model, selector, False, True)
            mhds.append(mhd)
        else:
            mhd = IMP.atom.read_pdb_or_mmcif(file, model, selector, True, True)
            mhds.append(mhd)

    for h_index in range(len(mhds)):
        ps = IMP.atom.get_by_type(mhds[h_index], IMP.atom.ATOM_TYPE)
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

def read_files(m, files, pdb_file_names, dat_files, particles_vec, exp_profiles, residue_level, heavy_atoms_only, multi_model_pdb, explicit_water, max_q, units):
    for file in files:
        # check if file exists
        if not os.path.exists(file):
            IMP.IMP_WARN("Can't open file " + file)
            return
        # 1. try as pdb or mmcif
        try:
            read_pdb(m, file, pdb_file_names, particles_vec, residue_level, heavy_atoms_only, multi_model_pdb, explicit_water)
        except IMP.ValueException:  # not a pdb file
            # 2. try as a dat profile file
            profile = IMP.saxs.Profile(file, False, max_q, units)
            if profile.size() == 0:
                print("Can't parse input file " + file)
                return
            else:
                dat_files.append(file)
                exp_profiles.append(profile)
                print("Profile read from file " + file + " size = " + str(profile.size()))

def trim_extension(file_name):
    if file_name[-4:] == '.':
        return file_name[:-4]
    return file_name

def read_multimodel_mmcif(in_file, model, selector, noradii):
    sp = IMP.PointerMember(selector)
    ret = read_mmcif(in_file, cif_nicename(in_file.get_name()), in_file.get_name(), model, selector, True, True, noradii)
    if not ret:
        raise ValueError(f"No molecule read from file {in_file.get_name()}")
    return ret

def read_mmcif(in_file, model, selector, select_first_model, noradii):
    sp = IMP.PointerMember(selector)
    ret = read_mmcif(in_file, cif_nicename(in_file.get_name()), in_file.get_name(), model, selector, False, select_first_model, noradii)
    if not ret:
        raise ValueError(f"No molecule read from file {in_file.get_name()}")
    return ret[0]
