"""
\file IMP/foxs/JmolWriter.cpp \brief outputs javascript for jmol display

Copyright 2007-2022 IMP Inventors. All rights reserved.
"""

from .ColorCoder import *
import IMP.saxs

class JmolWriter:
    display_selection_ = "frame 0#;restrict selection;\
        select selection and (protein, nucleic); ribbons only;\
        select selection and not (protein, nucleic); spacefill only;\
        if (!{*}.ribbons) { select selection and (protein, nucleic);spacefill only; };"

    MAX_DISPLAY_NUM_ = 30
    MAX_C2_ = 4.0

    @staticmethod
    def prepare_jmol_script(fps, particles_vec, filename):
        html_filename = filename + ".html"
        pdb_filename = filename + ".pdb"

        model_num = min(len(fps), JmolWriter.MAX_DISPLAY_NUM_)
        pdb_colors = JmolWriter.prepare_coloring_string(model_num)
        JmolWriter.prepare_PDB_file(fps, particles_vec, pdb_filename)

        with open(html_filename, "w") as outstream:
            init = "select all;" + JmolWriter.display_selection_
            selection_init = "define selection model =1"
            init += selection_init + "; " + pdb_colors + JmolWriter.display_selection_ + \
                "; background white; hide hydrogens;"

            # load applet with molecules
            outstream.write("<td width=350 height=350><div id=\"wrapper\" align=\"center\">\n")
            outstream.write("<script type=\"text/javascript\"> jmolApplet(350, 'load " +
                            pdb_filename + "; " + init + "');\n")
            outstream.write("</script> </div> </td> </tr> \n </table>\n")

            show_all_checkbox_str = JmolWriter.show_all_checkbox(model_num, True)

            outstream.precision(2)
            # output table
            showMolecule = True
            hex_color = "ZZZZZZ"
            outstream.write("<table align='center'>")
            outstream.write("<tr><th> PDB file </th> " +
                            "<th> " + show_all_checkbox_str + "</th>" +
                            "<th><center> &chi;<sup>2</sup> </th>" +
                            "<th><center><a href = \"https://modbase.compbio.ucsf.edu/foxs/help#c1c2\"> c<sub>1</sub> </a></th>" +
                            "<th><center><a href = \"https://modbase.compbio.ucsf.edu/foxs/help#c1c2\"> c<sub>2</sub> </a></th>" +
                            "<th><center>R<sub>g</sub></th>" +
                            "<th><center> # atoms </th> <th> fit file </th><th> png file </th></tr>\n")
            for i in range(len(fps)):
                hex_color = ColorCoder.html_hex_color(i)
                pdb_name = saxs.trim_extension(fps[i].get_pdb_file_name())
                profile_name = saxs.trim_extension(
                    basename(fps[i].get_profile_file_name().c_str()))
                fit_file_name = pdb_name + "_" + profile_name + ".fit"
                png_file_name = pdb_name + "_" + profile_name + ".png"
                rg = IMP.saxs.radius_of_gyration(particles_vec[fps[i].get_mol_index()])
                outstream.write("<tr><td>")
                if i < JmolWriter.MAX_DISPLAY_NUM_:
                    outstream.write("<font color=#" + hex_color + ">")
                outstream.write(pdb_name)
                if i < JmolWriter.MAX_DISPLAY_NUM_:
                    outstream.write("</font>")
                outstream.write("</td>\n<td><center>\n")
                if i > 0:
                    showMolecule = False
                if i < JmolWriter.MAX_DISPLAY_NUM_:
                    checkbox_string = JmolWriter.model_checkbox(i, showMolecule, True)
                    outstream.write(checkbox_string + "\n")
                outstream.write("</center></td><td><center> " + str(fps[i].get_chi_square()) +
                                "</center></td><td><center> " + str(fps[i].get_c1()))

                outstream.write("</center></td><td><center> ")
                if fps[i].get_c2() >= JmolWriter.MAX_C2_:
                    outstream.write("<a href = \"https://modbase.compbio.ucsf.edu/foxs/help#c1c2\"><font color=red>")
                    outstream.write(str(fps[i].get_c2()) + "!</font></a>")
                else:
                    outstream.write(str(fps[i].get_c2()))
                outstream.write("</center></td><td><center> " + str(rg) +
                                "</center></td><td><center> " +
                                str(len(particles_vec[fps[i].get_mol_index()])))
                outstream.write("</td><td><a href = \"dirname/" + fit_file_name + "\">" +
                                fit_file_name + "</a></td><td>" +
                                "<a href = \"dirname/" + png_file_name + "\">" +
                                png_file_name + "</a></td></tr>\n")

            outstream.write("</table>\n")
            outstream.write(JmolWriter.group_checkbox(model_num))

    def prepare_jmol_script(pdbs, particles_vec, filename):
        html_filename = filename + ".html"
        pdb_filename = filename + ".pdb"

        model_num = min(len(pdbs), JmolWriter.MAX_DISPLAY_NUM_)
        pdb_colors = prepare_coloring_string(model_num)

        prepare_PDB_file(particles_vec, pdb_filename)

        with open(html_filename, 'w') as outstream:
            outstream.write(prepare_gnuplot_init_selection_string(model_num, False))
            outstream.write(jsmol_script("/foxs/jsmol"))
            init = "select all;" + JmolWriter.display_selection_
            selection_init = "define selection model =1"

            init += selection_init + "; " + pdb_colors + JmolWriter.display_selection_ + \
                    "; background white; hide hydrogens;"

            outstream.write("<td width=350 height=350><div id=\"wrapper\" align=\"center\">\n")
            outstream.write("<script type=\"text/javascript\"> jmolApplet(350, 'load " +
                            pdb_filename + "; " + init + "');\n")
            outstream.write("</script> </div> </td> </tr> \n </table>\n")

            show_all_checkbox_str = show_all_checkbox(model_num, False)

            outstream.precision(2)
            outstream.setf(ios_fixed, ios_floatfield)

            showMolecule = True
            hex_color = "ZZZZZZ"
            outstream.write("<table align='center'>")
            outstream.write("<tr><th> PDB file </th>" +
                            "<th> " + show_all_checkbox_str + " </th>" +
                            "<th><center> R<sub>g</sub> </th>" +
                            "<th><center> # atoms </th> <th> Profile file</th></tr>\n")
            for i in range(len(pdbs)):
                hex_color = ColorCoder.html_hex_color(hex_color, i)
                pdb_name = saxs.trim_extension(pdbs[i])
                profile_name = pdbs[i] + ".dat"
                rg = IMP.saxs.radius_of_gyration(particles_vec[i])
                outstream.write("<tr><td>")
                if i < JmolWriter.MAX_DISPLAY_NUM_:
                    outstream.write("<font color=#" + hex_color + ">")
                outstream.write(pdb_name)
                if i < JmolWriter.MAX_DISPLAY_NUM_:
                    outstream.write("</font>")
                outstream.write("</td>\n<td><center>\n")
                if i > 0:
                    showMolecule = False
                if i < JmolWriter.MAX_DISPLAY_NUM_:
                    checkbox_string = model_checkbox(i, showMolecule, False)
                    outstream.write(checkbox_string + "\n")
                outstream.write("</center></td><td><center> " + str(rg) +
                                "</center></td><td><center> " + str(len(particles_vec[i])) +
                                "</td><td><a href = \"dirname/" + profile_name + "\">" +
                                profile_name + "</a></td></tr>\n")

            outstream.write("</table>\n")
            outstream.write(JmolWriter.group_checkbox(model_num))


    @staticmethod
    def prepare_coloring_string(model_num):
        coloring_string = ""
        for i in range(model_num):
            dec_color = ColorCoder.jmol_dec_color(i)
            coloring_string += "select model = " + str(i + 1) + "; color " + dec_color + ";"
        return coloring_string


    @staticmethod
    def prepare_PDB_file(fps, particles_vec, filename):
        with open(filename, 'w') as out_file:
            for i in range(len(fps)):
                if i >= JmolWriter.MAX_DISPLAY_NUM_:
                    break
                mol_index = fps[i].get_mol_index()

                # compute mean
                coordinates = [IMP.algebra.Vector3D(p.get_coordinates())
                            for p in particles_vec[mol_index]]
                m = sum(coordinates, IMP.algebra.Vector3D(0.0, 0.0, 0.0)) / len(coordinates)

                # output file
                out_file.write("MODEL    " + str(i + 1) + "\n")
                for j in range(len(particles_vec[mol_index])):
                    # centering
                    v = particles_vec[mol_index][j].get_coordinates() - m
                    ad = IMP.atom.Atom(particles_vec[mol_index][j])
                    rd = get_residue(ad)
                    c = get_chain(rd)
                    chain = c.get_id()[0]
                    out_file.write(IMP.atom.get_pdb_string(
                        v, ad.get_input_index(), ad.get_atom_type(),
                        rd.get_residue_type(), chain, rd.get_index(),
                        rd.get_insertion_code(), ad.get_occupancy(),
                        ad.get_temperature_factor(), ad.get_element()))

                out_file.write("ENDMDL\n")

    def prepare_PDB_file(particles_vec, filename):
        with open(filename, 'w') as out_file:
            for i in range(len(particles_vec)):
                if i >= JmolWriter.MAX_DISPLAY_NUM_:
                    break

                # compute mean
                coordinates = [IMP.algebra.Vector3D(p.get_coordinates())
                            for p in particles_vec[i]]
                m = sum(coordinates, IMP.algebra.Vector3D(0.0, 0.0, 0.0)) / len(coordinates)

                # output file
                out_file.write("MODEL    " + str(i + 1) + "\n")
                for j in range(len(particles_vec[i])):
                    # centering
                    v = particles_vec[i][j].get_coordinates() - m
                    ad = IMP.atom.Atom(particles_vec[i][j])
                    rd = get_residue(ad)
                    c = get_chain(rd)
                    chain = c.get_id()[0]
                    out_file.write(IMP.atom.get_pdb_string(
                        v, ad.get_input_index(), ad.get_atom_type(),
                        rd.get_residue_type(), chain, rd.get_index(),
                        rd.get_insertion_code(), ad.get_occupancy(),
                        ad.get_temperature_factor(), ad.get_element()))

                out_file.write("ENDMDL\n")

    @staticmethod
    def prepare_gnuplot_init_selection_string(model_num, exp):
        gnuplot_string = ""
        for i in range(1, model_num):
            if exp:
                gnuplot_string += "<script> gnuplot.hide_plot(\"jsoutput_1_plot_" + str(i + 2) + "\");</script>"
            else:
                gnuplot_string += "<script> gnuplot.hide_plot(\"jsoutput_1_plot_" + str(i + 1) + "\");</script>"
        return gnuplot_string

    def jmol_script(jmol_path):
        js_string = ""
        js_string += "<script src=\"" + jmol_path + "/Jmol.js\"></script>\n"
        js_string += "<script> jmolInitialize(\"" + jmol_path + "\", \"JmolAppletSigned.jar\"); </script>\n"
        return js_string

    @staticmethod
    def jsmol_script(jmol_path):
        js_string = ""
        js_string += "<script src=\"" + jmol_path + "/JSmol.min.js\"></script>\n"
        js_string += "<script src=\"" + jmol_path + "/Jmol2.js\"></script>\n"
        js_string += "<script type=\"text/javascript\">var Info = {}</script>\n"
        js_string += "<script> jmolInitialize(\"" + jmol_path + "\"); </script>\n"
        return js_string

    def model_checkbox(model_num, is_checked, exp):
        model_num_string = str(model_num + 1)
        model_num_string2 = model_num_string
        if exp:
            model_num_string2 = str(model_num + 2)
        checkbox_string = "<script>"
        checkbox_string += "\n jmolCheckbox("
        checkbox_string += "'javascript gnuplot.show_plot(\"jsoutput_1_plot_" + model_num_string2 + "\");"
        checkbox_string += "define selection selection, model=" + model_num_string + ";" + display_selection_ + "',"
        checkbox_string += "'javascript gnuplot.hide_plot(\"jsoutput_1_plot_" + model_num_string2 + "\");"
        checkbox_string += "define selection selection and not model=" + model_num_string + ";" + display_selection_ + "',''"
        if is_checked:
            checkbox_string += ",'isChecked'"
        checkbox_string += ") </script>\n"
        return checkbox_string

    def show_all_checkbox(model_num, exp):
        if model_num == 1:
            return "show/hide"
        show_string = "select all;define selection all;frame 0#;" + display_selection_
        hide_string = "define selection not all;frame 0#;restrict not all;"
        for i in range(model_num):
            plotnum = i + 1
            if exp:
                plotnum = i + 2
            plotnum_str = str(plotnum)
            show_string += "javascript gnuplot.show_plot(\"jsoutput_1_plot_" + plotnum_str + "\");"
            hide_string += "javascript gnuplot.hide_plot(\"jsoutput_1_plot_" + plotnum_str + "\");"
        ret = "<script>jmolCheckbox('" + show_string + "', '" + hide_string + "', \"\");</script> show all/hide all"
        return ret

    def group_checkbox(model_num):
        if model_num == 1:
            return ""
        ret = "<script> jmolSetCheckboxGroup(0"
        for i in range(model_num):
            ret += ", " + str(i + 1)
        ret += ");</script> \n"
        return ret
