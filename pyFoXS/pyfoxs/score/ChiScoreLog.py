"""
\file ChiScoreLog \brief scoring with log intensity

Copyright 2007-2022 IMP Inventors. All rights reserved.
"""

import numpy as np

class ChiScoreLog:
    def compute_scale_factor(self, exp_profile, model_profile, offset=None):
        sum1 = 0.0
        sum2 = 0.0
        profile_size = min(model_profile.size(), exp_profile.size())
        for k in range(profile_size):
            square_error = (exp_profile.error_[k] / exp_profile.intensity_[k])**2
            weight_tilda = 1 / square_error # model_profile.weight_[k] / square_error
            sum1 += weight_tilda * np.log(exp_profile.intensity_[k] / model_profile.intensity_[k])
            sum2 += weight_tilda
        return np.exp(sum1 / sum2)


    def compute_score(self, exp_profile, model_profile, offset):
        c = self.compute_scale_factor(exp_profile, model_profile, offset)
        chi_square = 0.0
        profile_size = min(model_profile.size(), exp_profile.size())
        for k in range(profile_size):
            square_error = np.log(exp_profile.error_[k])**2
            weight_tilda = 1 / square_error # model_profile.weight_[k] / square_error
            delta = np.log(exp_profile.intensity_[k]) - np.log(c * model_profile.intensity_[k])
            if abs(delta / np.log(exp_profile.intensity_[k])) >= 1.0e-15:
                chi_square += weight_tilda * delta**2
        chi_square /= profile_size
        return np.sqrt(chi_square)


    def compute_score_with_q_range(self, exp_profile, model_profile, min_q, max_q):
        c = self.compute_scale_factor(exp_profile, model_profile)
        chi_square = 0.0
        profile_size = min(model_profile.size(), exp_profile.size())
        interval_size = 0
        for k in range(profile_size):
            if exp_profile.q_[k] > max_q:
                break
            if exp_profile.q_[k] >= min_q:
                square_error = np.log(exp_profile.error_[k])**2
                weight_tilda = 1 / square_error # model_profile.weight_[k] / square_error
                delta = np.log(exp_profile.intensity_[k]) - np.log(c * model_profile.intensity_[k])
                if abs(delta / np.log(exp_profile.intensity_[k])) >= 1.0e-15:
                    chi_square += weight_tilda * delta**2
                    interval_size += 1
        if interval_size > 0:
            chi_square /= interval_size
        return np.sqrt(chi_square)
