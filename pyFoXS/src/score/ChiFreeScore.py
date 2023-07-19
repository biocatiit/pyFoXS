"""
\file ChiFreeScore.h   \brief Chi free SAXS score

Copyright 2007-2022 IMP Inventors. All rights reserved.
"""

import numpy as np

from ..utils.Profile import Profile

class ChiFreeScore:
    def __init__(self, K, ns):
        self.K_ = K
        self.ns_ = ns
        self.last_scale_updated_ = False
        self.last_scale_ = None

    def comp_function(self, a, b):
        return a[0] < b[0]

    def compute_score(self, exp_profile, model_profile, use_offset):
        if model_profile.size() != exp_profile.size():
            raise ValueError("ChiFreeScore.compute_score is supported only for profiles with the same q values!")

        self.last_scale_updated_ = False
        chis = np.zeros(self.K_, dtype=(float, float))
        bin_size = np.floor(exp_profile.size() / self.ns_)

        qs = np.zeros(self.ns_)
        errors = np.zeros(self.ns_)
        exp_intensities = np.zeros(self.ns_)
        model_intensities = np.zeros(self.ns_)

        for k in range(self.K_):
            for i in range(self.ns_):
                prob = np.random.uniform(0, 1)
                profile_index = int(i * bin_size + prob * bin_size)
                if profile_index < exp_profile.size():
                    qs[i] = exp_profile.q_[profile_index]
                    errors[i] = exp_profile.error_[profile_index]
                    exp_intensities[i] = exp_profile.intensity_[profile_index]
                    model_intensities[i] = model_profile.intensity_[profile_index]

            exp_profile_selection = Profile(qs, exp_intensities, errors)
            model_profile_selection = Profile(qs, model_intensities)

            offset = 0.0
            if use_offset:
                offset = self.compute_offset(exp_profile_selection, model_profile_selection)
            c = self.compute_scale_factor(exp_profile_selection, model_profile_selection, offset)

            chi_square = 0.0
            profile_size = min(model_profile_selection.size(), exp_profile_selection.size())
            for i in range(profile_size):
                square_error = exp_profile_selection.error_[i]**2
                weight_tilda = model_profile_selection.weight_[i] / square_error
                delta = exp_profile_selection.intensity_[i] - offset - c * model_profile_selection.intensity_[i]
                if abs(delta / exp_profile_selection.intensity_[i]) >= 1.0e-15:
                    chi_square += weight_tilda * delta**2
            chi_square /= profile_size
            chis[k] = (chi_square, c)

        n = self.K_ // 2
        chis_sorted = sorted(chis, key=lambda x: x[0])
        self.last_scale_updated_ = True
        self.last_scale_ = chis_sorted[n][1]
        return np.sqrt(chis_sorted[n][0])

    def compute_scale_factor(self, exp_profile, model_profile, offset):
        if self.last_scale_updated_:
            return self.last_scale_

        sum1 = 0.0
        sum2 = 0.0
        profile_size = min(model_profile.size(), exp_profile.size())
        for k in range(profile_size):
            square_error = exp_profile.error_[k]**2
            weight_tilda = model_profile.weight_[k] / square_error
            sum1 += weight_tilda * model_profile.intensity_[k] * (exp_profile.intensity_[k] - offset)
            sum2 += weight_tilda * model_profile.intensity_[k]**2

        return sum1 / sum2

    def compute_offset(self, exp_profile, model_profile):
        sum_iexp_imod = 0.0
        sum_imod = 0.0
        sum_iexp = 0.0
        sum_imod2 = 0.0
        sum_weight = 0.0
        profile_size = min(model_profile.size(), exp_profile.size())
        for k in range(profile_size):
            square_error = exp_profile.error_[k]**2
            weight_tilda = model_profile.weight_[k] / square_error
            sum_iexp_imod += weight_tilda * model_profile.intensity_[k] * exp_profile.intensity_[k]
            sum_imod += weight_tilda * model_profile.intensity_[k]
            sum_iexp += weight_tilda * exp_profile.intensity_[k]
            sum_imod2 += weight_tilda * model_profile.intensity_[k]**2
            sum_weight += weight_tilda

        offset = (sum_iexp_imod / sum_imod2) * sum_imod - sum_iexp
        offset /= (sum_weight - sum_imod * sum_imod / sum_imod2)
        return offset
