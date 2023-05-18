"""
\file IMP/foxs/Gnuplot.h   \brief A class for printing gnuplot scripts
for profile viewing

Copyright 2007-2022 IMP Inventors. All rights reserved.
"""

from .ColorCoder import ColorCoder
import IMP.saxs

class Gnuplot:
    @staticmethod
    def print_profile_script(pdb):
        # file names
        profile_file_name = pdb + ".dat"
        plt_file_name = trim_extension(pdb) + ".plt"
        png_file_name = trim_extension(pdb) + ".png"
        # output script
        with open(plt_file_name, "w") as plt_file:
            plt_file.write("set terminal png enhanced;set output '{}';\n".format(png_file_name))
            plt_file.write("set ylabel '';\nset format y '';\nset xtics nomirror;\n")
            plt_file.write("set ytics nomirror; set border 3;\n")
            plt_file.write("set style line 11 lc rgb '#808080' lt 1;\n")
            plt_file.write("set border 3 back ls 11;\n")
            plt_file.write("plot '{}' u 1:(log($2)) t 'FoXS' w lines lw 2.5 lc rgb '#e26261'\n".format(profile_file_name))
            plt_file.write("reset\n")

    @staticmethod
    def print_profile_script(pdbs):
        # ColorCoder.set_number(len(pdbs))
        hex_color = "#ZZZZZZ"
        with open("profiles.plt", "w") as plt_file:
            plt_file.write("set terminal png enhanced;set output \"profiles.png\"\n")
            plt_file.write("set ylabel '';\nset format y '';\nset xtics nomirror;\n")
            plt_file.write("set ytics nomirror; set border 3;\n")
            plt_file.write("set style line 11 lc rgb '#808080' lt 1;\n")
            plt_file.write("set border 3 back ls 11;\n")
            plt_file.write("plot ")
            for i, pdb in enumerate(pdbs):
                ColorCoder.html_hex_color(hex_color, i)
                profile_file_name = pdb + ".dat"
                plt_file.write("'{}' u 1:(log($2)) t \"{}\" w lines lw 2 lt {}"
                               .format(profile_file_name, trim_extension(pdb), i + 2))
                if i == len(pdbs) - 1:
                    plt_file.write("\n")
                else:
                    plt_file.write(",")
            plt_file.write("reset\n")

    @staticmethod
    def print_canvas_script(pdbs, max_num):
        # ColorCoder.set_number(len(pdbs))
        hex_color = "#ZZZZZZ"
        with open("canvas.plt", "w") as plt_file:
            plt_file.write("set terminal canvas solid butt size 400,350 fsize 10 lw 1.5 "
                           "fontscale 1 name \"jsoutput_1\" jsdir \".\"\n")
            plt_file.write("set output 'jsoutput.1.js'\n")
            plt_file.write("set xlabel 'q [Å^{-1}]';set ylabel 'log intensity'; "
                           "set format y '10^{%L}'; set logscale y\n"
                           "set xtics nomirror;set ytics nomirror; set border 3\n")
            plt_file.write("set style line 11 lc rgb '#808080' lt 1;unset key;\n")
            plt_file.write("set border 3 back ls 11;\n")
            plt_file.write("plot ")
            for i, pdb in enumerate(pdbs):
                if i < max_num:
                    ColorCoder.html_hex_color(hex_color, i)
                    profile_file_name = pdb + ".dat"
                    plt_file.write("'{}' u 1:2 "
                                   "w lines lw 2.5 lc rgb '#{}'".format(profile_file_name, hex_color))
                    if i == len(pdbs) - 1 or i == max_num - 1:
                        plt_file.write("\n")
                    else:
                        plt_file.write(",")
            plt_file.write("reset\n")

    @staticmethod
    def print_fit_script(fp):
        pdb_name = IMP.saxs.trim_extension(fp.get_pdb_file_name())
        profile_name = IMP.saxs.trim_extension(os.path.basename(fp.get_profile_file_name()))
        plt_file_name = pdb_name + "_" + profile_name + ".plt"
        png_file_name = pdb_name + "_" + profile_name + ".png"
        eps_file_name = pdb_name + "_" + profile_name + ".eps"
        fit_file_name = pdb_name + "_" + profile_name + ".fit"

        with open(plt_file_name, "w") as plt_file:
            plt_file.precision(3)
            plt_file.write("set terminal png enhanced; set output \"{}\"\n".format(png_file_name))
            plt_file.write("set lmargin 7; set rmargin 2;set multiplot\n")

            # Lower residuals plot
            plt_file.write("set origin 0,0;set size 1,0.3; set tmargin 0; set bmargin 3;"
                           "set ylabel 'Residual';set format y '';set xtics nomirror;"
                           "set xlabel 'q [Å^{-1}]';"
                           "set ytics nomirror; set border 3\n"
                           "set style line 11 lc rgb '#808080' lt 1;"
                           "set border 3 back ls 11;f(x)=1\n")
            plt_file.write("plot f(x) notitle lc rgb '#333333'"
                           ", '{}' u 1:(($2-$4)/$3) notitle w lines lw 2.5 lc rgb '#e26261'\n".format(fit_file_name))

            # Upper fit plot
            plt_file.write("set origin 0,0.3;set size 1,0.69; set bmargin 0; set tmargin 1;"
                           "set xlabel ''; set format x ''; "
                           "set ylabel 'I(q) log-scale';"
                           "set format y \"10^{%L}\"; set logscale y\n")
            plt_file.write("plot '{}' u 1:2 t 'Experimental' lc rgb '#333333' pt 6 ps 0.8"
                           ", '{}' u 1:4 t 'FoXS χ^2 = {}' w lines lw 2.5 lc rgb '#e26261'\n".format(fit_file_name, fit_file_name, fp.get_chi_square()))

            plt_file.write("unset multiplot\n")
            plt_file.write("reset\n")

            # Combined eps plot for paper
            plt_file.write("set terminal postscript eps size 3.5,2.62 color enhanced solid "
                           "linewidth 2.5 font 'Helvetica,22'; set output \"{}\"\n".format(eps_file_name))
            plt_file.write("set lmargin 2; set rmargin 2;set multiplot\n")

            # Lower residuals plot
            plt_file.write("set origin 0,0;set size 1,0.3; set tmargin 0; set bmargin 3;"
                           "set ylabel '';set format y '';\n")
            plt_file.write("set xtics nomirror font 'Helvetica,18';\n")
            plt_file.write("set ytics nomirror; set border 3\n")
            plt_file.write("set style line 11 lc rgb '#808080' lt 1;\n")
            plt_file.write("set border 3 back ls 11;f(x)=1\n")
            plt_file.write("plot f(x) notitle lc rgb '#333333'" + ", '{}' u 1:(($2-$4)/$3) notitle w lines lw 2.5 lc rgb '#e26261'\n".format(fit_file_name))

            # Upper fit plot
            plt_file.write("set origin 0,0.3;set size 1,0.69; set bmargin 0; set tmargin 1;")
            plt_file.write("set xlabel ''; set format x ''; set ylabel 'I(q) log-scale';\n")
            plt_file.write("set format y \"10^{%L}\"; set logscale y\n")
            plt_file.write("plot '{}' u 1:2 t 'Experimental' lc rgb '#333333' pt 6 ps 0.8".format(fit_file_name))
            plt_file.write(", '{}' u 1:4 t 'FoXS {/Symbol c}^2 = {}' w lines lw 2.5 lc rgb '#e26261'\n".format(fit_file_name, fp.get_chi_square()))

            plt_file.write("unset multiplot\n")
            plt_file.write("reset\n")


    def print_fit_script(fps):
        # ColorCoder.set_number(len(fps))
        hex_color = "#ZZZZZZ"
        plt_file = open("fit.plt", "w")

        # Terminal type
        plt_file.write("set terminal png enhanced;set output \"fit.png\";\n")

        # Formatting
        plt_file.write("set lmargin 7; set rmargin 2;set multiplot\n")
        plt_file.write("set origin 0,0;set size 1,0.3; set tmargin 0; set bmargin 3;")
        plt_file.write("set ylabel 'Residual';set format y '';set xtics nomirror;")
        plt_file.write("set xlabel 'q [Å^{-1}]';")
        plt_file.write("set ytics nomirror; set border 3\n")
        plt_file.write("set style line 11 lc rgb '#808080' lt 1;")
        plt_file.write("set border 3 back ls 11\n")

        # Residuals
        plt_file.write("f(x)=1\n")
        plt_file.write("plot f(x) notitle lc rgb '#333333'")
        for i in range(len(fps)):
            hex_color = ColorCoder.html_hex_color(i)
            pdb_name = IMP.saxs.trim_extension(fps[i].get_pdb_file_name())
            profile_name = IMP.saxs.trim_extension(os.path.basename(fps[i].get_profile_file_name()))
            fit_file_name = f"{pdb_name}_{profile_name}.fit"
            plt_file.write(f", '{fit_file_name}' u 1:(($2-$4)/$3) notitle w lines lw 2.5 lc rgb '#{hex_color}'")
        plt_file.write("\n")

        # Actual plots
        plt_file.write("set origin 0,0.3;set size 1,0.69; set bmargin 0; set tmargin 1;")
        plt_file.write("set xlabel ''; set format x '';")
        plt_file.write("set ylabel 'I(q) log-scale';")
        plt_file.write("set format y \"10^{%L}\"; set logscale y\n")
        for i in range(len(fps)):
            hex_color = ColorCoder.html_hex_color(i)
            pdb_name = IMP.saxs.trim_extension(fps[i].get_pdb_file_name())
            profile_name = IMP.saxs.trim_extension(os.path.basename(fps[i].get_profile_file_name()))
            fit_file_name = f"{pdb_name}_{profile_name}.fit"
            if i == 0:
                plt_file.write(f"plot '{fit_file_name}' u 1:2 t 'Experimental' lc rgb '#333333' pt 6 ps 0.8 ")
            plt_file.write(f", '{fit_file_name}' u 1:4 t '{pdb_name} χ^2 = {fps[i].get_chi_square()}' w lines lw 2.5 lc rgb '#{hex_color}'")
        plt_file.write("\n")

        plt_file.write("unset multiplot;reset\n")
        plt_file.close()

    def print_canvas_script(fps, max_num):
        # ColorCoder.set_number(len(fps))
        hex_color = "#ZZZZZZ"
        plt_file = open("canvas.plt", "w")

        plt_file.write("set terminal canvas solid butt size 400,350 fsize 10 lw 1.5 "
                    "fontscale 1 name \"jsoutput_1\" jsdir \".\"\n")
        plt_file.write("set output 'jsoutput.1.js'\n")

        plt_file.write("set multiplot\n")
        plt_file.write("set lmargin 7\n")
        plt_file.write("set origin 0,0;set size 1,0.3; set tmargin 0;"
                    "set xlabel 'q [Å^{-1}]';set ylabel 'Residual';set format y '';"
                    "set xtics nomirror;set ytics nomirror;unset key;"
                    "set border 3; set style line 11 lc rgb '#808080' lt 1;"
                    "set border 3 back ls 11\n")
        # residuals
        plt_file.write("f(x)=1\n")
        plt_file.write("plot f(x) lc rgb '#333333'")
        for i in range(min(len(fps), max_num)):
            hex_color = ColorCoder.html_hex_color(i)
            pdb_name = IMP.saxs.trim_extension(fps[i].get_pdb_file_name())
            profile_name = IMP.saxs.trim_extension(os.path.basename(fps[i].get_profile_file_name()))
            fit_file_name = f"{pdb_name}_{profile_name}.fit"
            plt_file.write(f", '{fit_file_name}' u 1:(($2-$4)/$3) w lines lw 2.5 lc rgb '#{hex_color}'")
        plt_file.write("\n")

        # actual plots
        plt_file.write("set origin 0,0.3;set size 1,0.69; set bmargin 0; set tmargin 1; "
                    "set xlabel ''; set format x ''; set ylabel 'log intensity'; "
                    "set format y '10^{%L}'; set logscale y\n")
        for i in range(min(len(fps), max_num)):
            hex_color = ColorCoder.html_hex_color(i)
            pdb_name = IMP.saxs.trim_extension(fps[i].get_pdb_file_name())
            profile_name = IMP.saxs.trim_extension(os.path.basename(fps[i].get_profile_file_name()))
            fit_file_name = f"{pdb_name}_{profile_name}.fit"
            if i == 0:
                plt_file.write(f"plot '{fit_file_name}' u 1:2 lc rgb '#333333' pt 6 ps 0.8 ")
            plt_file.write(f", '{fit_file_name}' u 1:4 w lines lw 2.5 lc rgb '#{hex_color}'")
        plt_file.write("\n")

        plt_file.write("unset multiplot;reset\n")
        plt_file.close()
