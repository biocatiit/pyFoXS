import IMP.test
import IMP.foxs
import shutil
import sys
import os
import re


class Tests(IMP.test.ApplicationTestCase):

    def test_nup133(self):
        """Test the Nup133 example"""
        cmds = ["foxs -h",
                "foxs 3KFO.pdb 23922_merge.dat",
                "foxs -g 3KFO.pdb 3KFO-fill.B99990005.pdb 23922_merge.dat",
                "gnuplot fit.plt"]
        nup133 = IMP.foxs.get_example_path('nup133')
        with IMP.test.temporary_working_directory():
            shutil.copytree(nup133, 'nup133')
            os.chdir('nup133')
            # Skip running gnuplot on machines that don't have it
            have_gnuplot = os.path.exists('/usr/bin/gnuplot') \
                           and sys.platform != 'win32'
            for c in cmds:
                if have_gnuplot or 'gnuplot' not in c:
                    self.run_shell_command(c)
            with open('3KFO_23922_merge.dat') as fh:
                lines = fh.readlines()
            self.assertIn('Chi^2 = 8.76', lines[1])
            with open('3KFO-fill.B99990005_23922_merge.dat') as fh:
                lines = fh.readlines()
            self.assertIn('Chi^2 = 1.31', lines[1])
            expected = ['3KFO_23922_merge.dat',
                        '3KFO-fill.B99990005_23922_merge.dat', 'fit.plt',
                        'profiles.plt']
            if have_gnuplot:
                expected.append('fit.png')
            for e in expected:
                os.unlink(e)

if __name__ == '__main__':
    IMP.test.main()
