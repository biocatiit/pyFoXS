"""
This is the program for SAXS profile computation and fitting.
see FOXS for webserver (salilab.org/foxs)
"""

import os
import sys
import argparse
import numpy as np

from src import __version__
from src.utils.utils import compute_profile, read_files, trim_extension
from src.utils.Profile import Profile
from src.utils.ProfileFitter import ProfileFitter
from src.score.ChiScoreLog import ChiScoreLog
from src.score.ChiFreeScore import ChiFreeScore
from src.score.RatioVolatilityScore import RatioVolatilityScore
from src.utils.FitParameters import FitParameters
from src.utils.Distribution import RadialDistributionFunction
from src.structure.FormFactorTable import get_default_form_factor_table, FormFactorTable, FormFactorType
from src.structure.Atom import compute_max_distance

def main():
    """
    Main function to run pyFoXS in the terminal
    """
    np.random.seed(42)
    profile_size = 500
    max_q = 0.5 # change after read
    min_c1 = 0.99
    max_c1 = 1.05
    min_c2 = -2.0
    max_c2 = 4.0
    heavy_atoms_only = True
    residue_level = False
    background_adjustment_q = 0.0
    use_offset = False
    write_partial_profile = False
    multi_model_pdb = 1
    units = 1 # determine automatically
    vr_score = False
    score_log = False
    # gnuplot_script = False
    explicit_water = False
    desc_prefix = ""
    form_factor_table_file = ""
    beam_profile_file = ""
    ab_initio = False
    vacuum = False
    chi_free = 0
    pr_dmax = 0.0

    hidden = argparse.ArgumentParser(add_help=False)
    hidden.add_argument("input_files", nargs="*", help="input PDB and profile files")
    hidden.add_argument("--form_factor_table", "-f", help="ff table name", default=form_factor_table_file)
    hidden.add_argument("--explicit_water", help="use waters from input PDB (default = False)", action="store_true")
    hidden.add_argument("--beam_profile", help="beam profile file name for desmearing", default=beam_profile_file)
    hidden.add_argument("--ab_initio", help="compute profile for a bead model with constant form factor (default = False)", action="store_true")
    hidden.add_argument("--vacuum", help="compute profile in vacuum (default = False)", action="store_true")
    hidden.add_argument("--chi_free", help="compute chi-free instead of chi, specify iteration number (default = 0)", type=int, default=chi_free)
    hidden.add_argument("--pr_dmax", help="Dmax value for P(r) calculation. P(r) is calculated only if pr_dmax > 0", type=float, default=pr_dmax)

    parser = argparse.ArgumentParser(description=desc_prefix + ".\n\nOptions", parents=[hidden])
    parser.add_argument("--version", help="Show version info and exit.", action="store_true")
    parser.add_argument("--profile_size", "-s", help="number of points in the profile", type=int, default=profile_size)
    parser.add_argument("--max_q", "-q", help="max q value", type=float, default=max_q)
    parser.add_argument("--min_c1", help="min c1 value", type=float, default=min_c1)
    parser.add_argument("--max_c1", help="max c1 value", type=float, default=max_c1)
    parser.add_argument("--min_c2", help="min c2 value", type=float, default=min_c2)
    parser.add_argument("--max_c2", help="max c2 value", type=float, default=max_c2)
    parser.add_argument("--hydrogens", "-hyd", help="explicitly consider hydrogens in PDB files (default = False)", action="store_true")
    parser.add_argument("--residues", "-r", help="fast coarse grained calculation using CA atoms only (default = False)", action="store_true")
    parser.add_argument("--background_q", "-b", help="background adjustment, not used by default. if enabled, recommended q value is 0.2", type=float, default=background_adjustment_q)
    parser.add_argument("--offset", "-o", help="use offset in fitting (default = False)", action="store_true")
    parser.add_argument("--write-partial-profile", "-p", help="write partial profile file (default = False)", action="store_true")
    parser.add_argument("--multi-model-pdb", "-m", help="1 - read the first MODEL only (default), 2 - read each MODEL into a separate structure, 3 - read all models into a single structure", type=int, default=multi_model_pdb)
    parser.add_argument("--units", "-u", help="1 - unknown --> determine automatically (default) 2 - q values are in 1/A, 3 - q values are in 1/nm", type=int, default=units)
    parser.add_argument("--volatility_ratio", "-v", help="calculate volatility ratio score (default = False)", action="store_true")
    parser.add_argument("--score_log", "-l", help="use log(intensity) in fitting and scoring (default = False)", action="store_true")
    # parser.add_argument("--gnuplot_script", "-g", help="print gnuplot script for gnuplot viewing (default = False)", action="store_true")

    args = parser.parse_args()

    if args.version:
        print(f"pyFoXS Version: {__version__}")
        return

    print("Usage: <pdb_file1> <pdb_file2> ... <profile_file1> <profile_file2> ...\n"
      "\nAny number of input PDBs and profiles is supported.\n"
      "Each PDB will be fitted against each profile.\n\n"
      "This program is part of IMP, the Integrative Modeling Platform.\n")

    fit = True
    files = []
    pdb_files = []
    dat_files = []

    if args.input_files:
        files = args.input_files

    if not files:
        print("WARNING: You need to specify a file to the program.\n")
        return

    if args.hydrogens:
        heavy_atoms_only = False

    if args.residues:
        residue_level = True

    if args.offset:
        use_offset = True

    if args.write_partial_profile:
        write_partial_profile = True

    if args.score_log:
        score_log = True

    if args.explicit_water:
        explicit_water = True

    # no water layer or fitting in ab initio mode for now
    if args.ab_initio:
        ab_initio = True
        fit = False

    if args.vacuum:
        vacuum = True

    if args.volatility_ratio:
        vr_score = True

    form_factor_table_file = args.form_factor_table
    beam_profile_file = args.beam_profile
    chi_free = args.chi_free
    pr_dmax = args.pr_dmax

    profile_size = args.profile_size
    max_q = args.max_q
    min_c1 = args.min_c1
    max_c1 = args.max_c1
    min_c2 = args.min_c2
    max_c2 = args.max_c2
    background_adjustment_q = args.background_q
    multi_model_pdb = args.multi_model_pdb
    units = args.units

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
    fps = []

    for i, part in enumerate(particles_vec):
        print("Computing profile for", pdb_files[i], len(part), "atoms")
        profile = compute_profile(part, 0.0, max_q, delta_q, ft, ff_type,
                                not explicit_water, fit, reciprocal, ab_initio, vacuum,
                                beam_profile_file)

        # save the profile
        profiles.append(profile)
        # write profile file
        profile_file_name = pdb_files[i] + ".dat"
        if write_partial_profile:
            profile.write_partial_profiles(profile_file_name)
        else:  # write normal profile
            profile.add_errors()
            profile.write_SAXS_file(profile_file_name)

        # calculate P(r)
        if pr_dmax > 0.0:
            pr = RadialDistributionFunction(0.5)
            profile.profile_2_distribution(pr, pr_dmax)
            pr.normalize()
            pr_file_name = pdb_files[i] + ".pr"
            with open(pr_file_name, "w") as pr_file:
                pr_file.write("Distance distribution\n")
                for item in pr.distribution:
                    pr_file.write(str(item)+"\n")

        # 3. fit experimental profiles
        for j, dat_file in enumerate(dat_files):
            exp_saxs_profile = exp_profiles[j]
            fit_file_name2 = trim_extension(pdb_files[i]) + "_" + \
                trim_extension(os.path.basename(dat_file)) + ".dat"
            fp = FitParameters()
            if score_log:
                pf = ProfileFitter(exp_saxs_profile, scoring_function=ChiScoreLog())
                fp = pf.fit_profile(profile, min_c1, max_c1, min_c2, max_c2,
                                    use_offset, fit_file_name2)
            else:
                if vr_score:
                    pf = ProfileFitter(exp_saxs_profile, scoring_function=RatioVolatilityScore())
                    fp = pf.fit_profile(profile, min_c1, max_c1, min_c2, max_c2,
                                        use_offset, fit_file_name2)
                else:
                    # default path
                    pf = ProfileFitter(exp_saxs_profile) # scoring_function=ChiScore() by default
                    fp = pf.fit_profile(profile, min_c1, max_c1, min_c2, max_c2,
                                        use_offset, fit_file_name2)
                    if chi_free > 0:
                        dmax = compute_max_distance(particles_vec[i])
                        ns = int(round(exp_saxs_profile.max_q_ * dmax / np.pi))
                        K = chi_free
                        cfs = ChiFreeScore(K, ns)
                        resampled_profile = Profile(qmin=exp_saxs_profile.min_q_, qmax=exp_saxs_profile.max_q_,
                                                    delta=exp_saxs_profile.delta_q_, constructor=0)

                        pf.resample(profile, resampled_profile)
                        chi_free = cfs.compute_score(exp_saxs_profile, resampled_profile, use_offset)
                        fp.chi_square = chi_free
            fp.pdb_file_name = pdb_files[i]
            fp.profile_file_name = dat_file
            fp.mol_index = i
            fp.show(sys.stdout)
            fps.append(fp)

    fps.sort(key=lambda x: x.chi_square)
    return

if __name__ == "__main__":
    main()
