"""
This is the program for SAXS profile computation and fitting.
see FOXS for webserver (salilab.org/foxs)
"""

import os
import sys
import numpy as np

pyfoxs_path = os.path.abspath(os.path.join('.', __file__, '..', '..'))
if pyfoxs_path not in os.sys.path:
    os.sys.path.append(pyfoxs_path)
from pyfoxs import __version__
from pyfoxs.utils.utils import compute_profile, read_files, trim_extension
from pyfoxs.utils.Profile import Profile
from pyfoxs.utils.ProfileFitter import ProfileFitter
from pyfoxs.score.ChiScoreLog import ChiScoreLog
from pyfoxs.score.ChiFreeScore import ChiFreeScore
from pyfoxs.score.RatioVolatilityScore import RatioVolatilityScore
from pyfoxs.utils.FitParameters import FitParameters
from pyfoxs.utils.Distribution import RadialDistributionFunction
from pyfoxs.structure.FormFactorTable import get_default_form_factor_table, FormFactorTable, FormFactorType
from pyfoxs.structure.Atom import compute_max_distance

def pyfoxs(files, profile_size=500, max_q=0.5, min_c1=0.99, max_c1=1.05,
    min_c2=-2.0, max_c2=4.0, hydrogens=False, residues=False,
    background_adjustment_q=0.0, use_offset=False, write_partial_profile=False,
    multi_model_pdb=1, units=1, vr_score=False, score_log=False,
    explicit_water=False, form_factor_table_file="", beam_profile_file="",
    ab_initio=False, vacuum=False, chi_free=0, pr_dmax=0.0, write_output=True,
    random_seed=None):
    """
    Main function to run pyFoXS
    """
    if random_seed is not None:
        np.random.seed(random_seed)

    heavy_atoms_only = not hydrogens
    residue_level = residues


    if not ab_initio:
        fit = True
    else:
        fit = False

    pdb_files = []
    dat_files = []

    if not files:
        print("WARNING: You need to specify a file to the program.\n")
        return

    if multi_model_pdb not in (1, 2, 3):
        print(f"Incorrect option for multi_model_pdb {multi_model_pdb}")
        print("Use 1 to read first MODEL only")
        print("    2 to read each MODEL into a separate structure,")
        print("    3 to read all models into a single structure")
        print("Default value of 1 is used")
        multi_model_pdb = 1

    if units not in (1, 2, 3):
        print(f"Incorrect option for units {units}")
        print("Use 1 for unknown units, 2 for 1/A, 3 for 1/nm")
        print("Default value of 1 is used")
        units = 1

    # determine form factor type
    ff_type = FormFactorType.HEAVY_ATOMS

    if not heavy_atoms_only:
        ff_type = FormFactorType.ALL_ATOMS

    if residue_level:
        ff_type = FormFactorType.CA_ATOMS

    # 1. read pdbs and profiles, prepare particles
    particles_vec = []
    exp_profiles = []
    m = None # Model()

    read_files(m, files, pdb_files, dat_files, particles_vec, exp_profiles,
            residue_level, heavy_atoms_only, multi_model_pdb, explicit_water,
            max_q, units)

    if background_adjustment_q > 0.0:
        for profile in exp_profiles:
            profile.background_adjust(background_adjustment_q)

    if len(exp_profiles) == 0 and not write_partial_profile:
        fit = False

    if max_q == 0.0:
        if len(exp_profiles) > 0:
            for profile in exp_profiles:
                if profile.max_q_ > max_q:
                    max_q = profile.max_q_
        else:
            max_q = 0.5

    delta_q = max_q / profile_size

    # read in or use default form factor table
    reciprocal = False
    ft = None

    if len(form_factor_table_file) > 0:
        # reciprocal space calculation, requires form factor file
        ft = FormFactorTable(form_factor_table_file, 0.0, max_q, delta_q)
        reciprocal = True
    else:
        ft = get_default_form_factor_table()

    # 2. compute profiles for input pdbs
    profiles = []
    fit_profiles = []
    fps = []
    prs = []

    for i, part in enumerate(particles_vec):
        print("Computing profile for", pdb_files[i], len(part), "atoms")
        profile = compute_profile(part, 0.0, max_q, delta_q, ft, ff_type,
                                not explicit_water, fit, reciprocal, ab_initio, vacuum,
                                beam_profile_file)

        # save the profile
        profiles.append(profile)
        # write profile file
        profile_file_name = pdb_files[i] + ".dat"
        if write_partial_profile and write_output:
            profile.write_partial_profiles(profile_file_name)
        else:  # write normal profile
            profile.add_errors()

            if write_output:
                profile.write_SAXS_file(profile_file_name)

        # calculate P(r)
        if pr_dmax > 0.0:
            pr = RadialDistributionFunction(0.5)
            profile.profile_2_distribution(pr, pr_dmax)
            pr.normalize()
            pr_file_name = pdb_files[i] + ".pr"
            prs.append(pr)
            if write_output:
                with open(pr_file_name, "w") as pr_file:
                    pr_file.write("Distance distribution\n")
                    for item in pr.distribution:
                        pr_file.write(str(item)+"\n")

        # 3. fit experimental profiles
        for j, dat_file in enumerate(dat_files):
            exp_saxs_profile = exp_profiles[j]

            if write_output:
                fit_file_name2 = trim_extension(pdb_files[i]) + "_" + \
                    trim_extension(os.path.basename(dat_file)) + ".dat"
            else:
                fit_file_name2 = ""

            fp = FitParameters()
            if score_log:
                pf = ProfileFitter(exp_saxs_profile, scoring_function=ChiScoreLog())
                fp, fit_profile = pf.fit_profile(profile, min_c1, max_c1, min_c2, max_c2,
                                    use_offset, fit_file_name2)
            else:
                if vr_score:
                    pf = ProfileFitter(exp_saxs_profile, scoring_function=RatioVolatilityScore())
                    fp, fit_profile = pf.fit_profile(profile, min_c1, max_c1, min_c2, max_c2,
                                        use_offset, fit_file_name2)
                else:
                    # default path
                    pf = ProfileFitter(exp_saxs_profile) # scoring_function=ChiScore() by default
                    fp, fit_profile = pf.fit_profile(profile, min_c1, max_c1, min_c2, max_c2,
                                        use_offset, fit_file_name2)
                    if chi_free > 0:
                        dmax = compute_max_distance(particles_vec[i])
                        ns = int(round(exp_saxs_profile.max_q_ * dmax / np.pi))
                        K = chi_free
                        cfs = ChiFreeScore(K, ns)
                        resampled_profile = Profile(qmin=exp_saxs_profile.min_q_, qmax=exp_saxs_profile.max_q_,
                                                    delta=exp_saxs_profile.delta_q_, constructor=0)

                        pf.resample(profile, resampled_profile)
                        chi_free, fit_profile = cfs.compute_score(exp_saxs_profile, resampled_profile, use_offset)
                        fp.chi_square = chi_free

            fit_profiles.append(fit_profile)

            fp.pdb_file_name = pdb_files[i]
            fp.profile_file_name = dat_file
            fp.mol_index = i
            fp.show(sys.stdout)
            fps.append(fp)

    fps.sort(key=lambda x: x.chi_square)

    return profiles, fit_profiles, fps, prs
