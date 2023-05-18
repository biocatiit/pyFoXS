"""
This is the program for SAXS profile computation and fitting.
see FOXS for webserver (salilab.org/foxs)
"""

import sys
import argparse
import IMP.saxs
from src.internal.Gnuplot import *
from src.internal.JmolWriter import *
from src.internal.utils import *

def main():
    __version__ = "0.1"
    profile_size = 500
    max_q = 0.0 # change after read
    min_c1 = 0.99
    max_c1 = 1.05
    min_c2 = -0.5
    max_c2 = 2.0
    heavy_atoms_only = True
    residue_level = False
    background_adjustment_q = 0.0
    desc_prefix = ""
    use_offset = False
    write_partial_profile = False
    multi_model_pdb = 1
    units = 1 # determine automatically
    vr_score = False
    score_log = False
    gnuplot_script = False
    explicit_water = False
    form_factor_table_file = ""
    beam_profile_file = ""
    ab_initio = False
    vacuum = False
    javascript = False
    chi_free = 0
    pr_dmax = 0.0
    print("Usage: <pdb_file1> <pdb_file2> ... <profile_file1> <profile_file2> ...\n"
      "\nAny number of input PDBs and profiles is supported.\n"
      "Each PDB will be fitted against each profile.\n\n"
      "This program is part of IMP, the Integrative Modeling Platform.")

    hidden = argparse.ArgumentParser(add_help=False)
    hidden.add_argument("input_files", nargs="*", help="input PDB and profile files")
    hidden.add_argument("--form_factor_table", "-f", help="ff table name", default=form_factor_table_file)
    hidden.add_argument("--explicit_water", help="use waters from input PDB (default = False)", action="store_true")
    hidden.add_argument("--beam_profile", help="beam profile file name for desmearing", default=beam_profile_file)
    hidden.add_argument("--ab_initio", help="compute profile for a bead model with constant form factor (default = False)", action="store_true")
    hidden.add_argument("--vacuum", help="compute profile in vacuum (default = False)", action="store_true")
    hidden.add_argument("--javascript", help="output JavaScript for browser viewing of the results (default = False)", action="store_true")
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
    parser.add_argument("--gnuplot_script", "-g", help="print gnuplot script for gnuplot viewing (default = False)", action="store_true")

    args = parser.parse_args()

    if args.version:
        print(f"Version: \"{__version__}\"")
        return 0

    fit = True
    files = []
    pdb_files = []
    dat_files = []

    if args.input_files:
        files = args.input_files

    if not files:
        print(parser)
        return 0

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

    if args.gnuplot_script:
        gnuplot_script = True

    if args.explicit_water:
        explicit_water = True

    # no water layer or fitting in ab initio mode for now
    if args.ab_initio:
        ab_initio = True
        fit = False

    if args.vacuum:
        vacuum = True

    if args.javascript:
        javascript = True

    if args.volatility_ratio:
        vr_score = True

    if multi_model_pdb != 1 and multi_model_pdb != 2 and multi_model_pdb != 3:
        print(f"Incorrect option for multi_model_pdb {multi_model_pdb}")
        print("Use 1 to read first MODEL only")
        print("    2 to read each MODEL into a separate structure,")
        print("    3 to read all models into a single structure")
        print("Default value of 1 is used")
        multi_model_pdb = 1

    if units != 1 and units != 2 and units != 3:
        print(f"Incorrect option for units {units}")
        print("Use 1 for unknown units, 2 for 1/A, 3 for 1/nm")
        print("Default value of 1 is used")
        units = 1

    # IMP::benchmark::Profiler pp("prof_out");

    # determine form factor type
    ff_type = IMP.saxs.HEAVY_ATOMS

    if not heavy_atoms_only:
        ff_type = IMP.saxs.ALL_ATOMS

    if residue_level:
        ff_type = IMP.saxs.CA_ATOMS

    # 1. read pdbs and profiles, prepare particles
    particles_vec = []
    exp_profiles = []
    m = IMP.Model()

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
                if profile.get_max_q() > max_q:
                    max_q = profile.get_max_q()
        else:
            max_q = 0.5

    delta_q = max_q / profile_size

    # read in or use default form factor table
    reciprocal = False
    ft = None

    if len(form_factor_table_file) > 0:
        # reciprocal space calculation, requires form factor file
        ft = IMP.saxs.FormFactorTable(form_factor_table_file, 0.0, max_q, delta_q)
        reciprocal = True
    else:
        ft = IMP.saxs.get_default_form_factor_table()

    # 2. compute profiles for input pdbs
    profiles = []
    fps = []

    for i in range(len(particles_vec)):
        print("Computing profile for", pdb_files[i], len(particles_vec[i]), "atoms")
        profile = compute_profile(particles_vec[i], 0.0, max_q, delta_q, ft, ff_type,
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
            if gnuplot_script:
                Gnuplot.print_profile_script(pdb_files[i])

        # calculate P(r)
        if pr_dmax > 0.0:
            pr = IMP.saxs.RadialDistributionFunction(0.5)
            profile.profile_2_distribution(pr, pr_dmax)
            pr.normalize()
            pr_file_name = pdb_files[i] + ".pr"
            with open(pr_file_name, "w") as pr_file:
                pr.show(pr_file)

        # 3. fit experimental profiles
        for j in range(len(dat_files)):
            exp_saxs_profile = exp_profiles[j]
            fit_file_name2 = trim_extension(pdb_files[i]) + "_" + \
                trim_extension(os.path.basename(dat_files[j])) + ".dat"

            fp = IMP.saxs.FitParameters()
            if score_log:
                pf = IMP.saxs.ProfileFitterChiLog(exp_saxs_profile)
                fp = pf.fit_profile(profile, min_c1, max_c1, min_c2, max_c2,
                                    use_offset, fit_file_name2)
            else:
                if vr_score:
                    pf = IMP.saxs.ProfileFitterRatioVolatility(exp_saxs_profile)
                    fp = pf.fit_profile(profile, min_c1, max_c1, min_c2, max_c2,
                                        use_offset, fit_file_name2)
                else:
                    pf = IMP.saxs.ProfileFitterChi(exp_saxs_profile)
                    fp = pf.fit_profile(profile, min_c1, max_c1, min_c2, max_c2,
                                        use_offset, fit_file_name2)
                    if chi_free > 0:
                        dmax = IMP.saxs.compute_max_distance(particles_vec[i])
                        ns = int(round(exp_saxs_profile.get_max_q() * dmax / math.pi))
                        K = chi_free
                        cfs = IMP.saxs.ChiFreeScore(ns, K)
                        cfs.set_was_used(True)
                        resampled_profile = IMP.saxs.Profile(exp_saxs_profile.get_min_q(), exp_saxs_profile.get_max_q(),
                                                    exp_saxs_profile.get_delta_q())
                        pf.resample(profile, resampled_profile)
                        chi_free = cfs.compute_score(exp_saxs_profile, resampled_profile)
                        fp.set_chi_square(chi_free)
            fp.set_pdb_file_name(pdb_files[i])
            fp.set_profile_file_name(dat_files[j])
            fp.set_mol_index(i)
            fp.show(sys.stdout)
            if gnuplot_script:
                Gnuplot.print_fit_script(fp)
            fps.append(fp)

    fps.sort(key=lambda x: x.get_score())

    if len(pdb_files) > 1 and gnuplot_script:
        Gnuplot.print_profile_script(pdb_files)
        if len(dat_files) > 0:
            Gnuplot.print_fit_script(fps)

    if javascript:
        if len(dat_files) > 0:
            Gnuplot.print_canvas_script(fps, JmolWriter.MAX_DISPLAY_NUM_)
            JmolWriter.prepare_jmol_script(fps, particles_vec, "jmoltable")
        else:
            Gnuplot.print_canvas_script(pdb_files, JmolWriter.MAX_DISPLAY_NUM_)
            JmolWriter.prepare_jmol_script(pdb_files, particles_vec, "jmoltable")

    return 0

if __name__ == "__main__":
    main()