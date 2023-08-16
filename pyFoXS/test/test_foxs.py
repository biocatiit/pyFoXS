import os

import pytest
import numpy as np

pyfoxs_path = os.path.abspath(os.path.join('.', __file__, '..', '..'))
if pyfoxs_path not in os.sys.path:
    os.sys.path.append(pyfoxs_path)

import pyfoxs
from pyfoxs.pyfoxs_api import pyfoxs
from pyfoxs.utils.utils import read_files
from pyfoxs.utils.Profile import Profile

def test_no_exp_data_lys_pdb():
    lys_model = os.path.join('.', 'data', '6lyz.pdb')
    profiles, fit_profiles, fps, prs = pyfoxs([lys_model,], write_output=False,
        random_seed=42)

    comp_data = os.path.join('.', 'data', 'py_6lyz.dat')

    comp_prof = Profile(file_name=comp_data, fit_file=False, constructor=1)

    calc_prof = profiles[0]

    assert np.allclose(comp_prof.q_, calc_prof.q_)
    assert np.allclose(comp_prof.intensity_, calc_prof.intensity_)
    assert np.allclose(comp_prof.error_, calc_prof.error_)

def test_exp_data_lys_pdb():
    lys_model = os.path.join('.', 'data', '6lyz.pdb')
    lys_exp = os.path.join('.', 'data', 'lyzexp.dat')
    profiles, fit_profiles, fps, prs = pyfoxs([lys_model, lys_exp],
        write_output=False, random_seed=42)

    comp_data = os.path.join('.', 'data', 'py_6lyz_lyzexp.fit')

    comp_prof = Profile(file_name=comp_data, fit_file=True, constructor=1)

    calc_prof = fit_profiles[0]
    fit_params = fps[0]

    assert np.allclose(comp_prof.q_, calc_prof.q_)
    assert np.allclose(comp_prof.intensity_, calc_prof.intensity_)
    assert np.isclose(0.20220989389435595, fit_params.chi_square)
    assert np.isclose(1.0113119999999998, fit_params.c1)
    assert np.isclose(0.5871999999999999, fit_params.c2)

