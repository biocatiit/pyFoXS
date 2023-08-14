"""
\file RatioVolatilityScore.h   \brief Chi free SAXS score

Copyright 2007-2022 IMP Inventors. All rights reserved.
"""

import numpy as np

from .ChiScore import ChiScore

class RatioVolatilityScore:
    def __init__(self, dmax=400):
        self.dmax_ = dmax

    def compute_score(self, exp_profile, model_profile, use_offset):
        if model_profile.size() != exp_profile.size():
            raise ValueError("RatioVolatilityScore::compute_score is supported only for profiles with the same q values!")

        bin_size = np.pi / self.dmax_
        number_of_bins = int(exp_profile.max_q_ / bin_size)
        number_of_points_in_bin = exp_profile.size() / number_of_bins

        ratio = np.zeros(number_of_bins)
        for i in range(number_of_bins):
            index1 = int(i * number_of_points_in_bin)
            index2 = int((i + 1) * number_of_points_in_bin)
            intensity1 = 0.0
            intensity2 = 0.0
            for j in range(index1, index2):
                intensity1 += exp_profile.intensity_[j]
                intensity2 += model_profile.intensity_[j]
            intensity1 /= (index2 - index1)
            intensity2 /= (index2 - index1)
            ratio[i] = intensity1 / intensity2

        vr = 0.0
        for i in range(number_of_bins - 1):
            vr += 2 * abs(ratio[i] - ratio[i + 1]) / (ratio[i] + ratio[i + 1])

        return 100 * vr / number_of_bins

    def compute_scale_factor(self, exp_profile, model_profile, offset):
        cs = ChiScore()
        return cs.compute_scale_factor(exp_profile, model_profile, offset)
