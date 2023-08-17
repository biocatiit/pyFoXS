"""
This is the program for SAXS profile computation and fitting.
see FOXS for webserver (salilab.org/foxs)
"""

import os
import sys
import argparse
import numpy as np

from pyfoxs import __version__
# from pyfoxs.utils.utils import compute_profile, read_files, trim_extension
# from pyfoxs.utils.Profile import Profile
# from pyfoxs.utils.ProfileFitter import ProfileFitter
# from pyfoxs.score.ChiScoreLog import ChiScoreLog
# from pyfoxs.score.ChiFreeScore import ChiFreeScore
# from pyfoxs.score.RatioVolatilityScore import RatioVolatilityScore
# from pyfoxs.utils.FitParameters import FitParameters
# from pyfoxs.utils.Distribution import RadialDistributionFunction
# from pyfoxs.structure.FormFactorTable import get_default_form_factor_table, FormFactorTable, FormFactorType
# from pyfoxs.structure.Atom import compute_max_distance
from pyfoxs import pyfoxs_api

def main():
    """
    Main function to run pyFoXS in the terminal
    """
    # np.random.seed(42)
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
    parser.add_argument("--random_seed", "-rs", help="Input an integer random seed to use.", default=None)
    # parser.add_argument("--gnuplot_script", "-g", help="print gnuplot script for gnuplot viewing (default = False)", action="store_true")

    args = parser.parse_args()

    if args.version:
        print(f"pyFoXS Version: {__version__}")
        return

    print("Usage: <pdb_file1> <pdb_file2> ... <profile_file1> <profile_file2> ...\n"
      "\nAny number of input PDBs and profiles is supported.\n"
      "Each PDB will be fitted against each profile.\n")
    files = []

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

    if args.vacuum:
        vacuum = True

    if args.volatility_ratio:
        vr_score = True

    random_seed = args.random_seed
    if random_seed is not None:
        random_seed = int(random_seed)

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

    hydrogens = not heavy_atoms_only

    pyfoxs_api.pyfoxs(files, profile_size=profile_size, max_q=max_q,
        min_c1=min_c1, max_c1=max_c1, min_c2=min_c2, max_c2=max_c2,
        hydrogens=hydrogens, residues=residue_level,
        background_adjustment_q=background_adjustment_q, use_offset=use_offset,
        write_partial_profile=write_partial_profile,
        multi_model_pdb=multi_model_pdb, units=units, vr_score=vr_score,
        score_log=score_log, explicit_water=explicit_water,
        form_factor_table_file=form_factor_table_file,
        beam_profile_file=beam_profile_file, ab_initio=ab_initio,
        vacuum=vacuum, chi_free=chi_free, pr_dmax=pr_dmax, write_output=True,
        random_seed=random_seed)

if __name__ == "__main__":
    main()
