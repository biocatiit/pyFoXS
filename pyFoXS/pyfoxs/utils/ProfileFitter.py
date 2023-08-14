"""
\file IMP/saxs/ProfileFitter.h \brief a class for fitting two profiles

Copyright 2007-2022 IMP Inventors. All rights reserved.
"""

from .FitParameters import FitParameters
from .Profile import Profile
from ..score.ChiScore import ChiScore

class ProfileFitter:
    def __init__(self, exp_profile, scoring_function=ChiScore()):
        self.exp_profile_ = exp_profile
        self.scoring_function_ = scoring_function

    def compute_scale_factor(self, model_profile, offset=0.0):
        return self.scoring_function_.compute_scale_factor(self.exp_profile_, model_profile, offset)

    def compute_offset(self, model_profile):
        return self.scoring_function_.compute_offset(self.exp_profile_, model_profile)

    def get_profile(self):
        return self.exp_profile_

    def resample(self, model_profile, resampled_profile):
        model_profile.resample(self.exp_profile_, resampled_profile)

    def search_fit_parameters(self, partial_profile, min_c1, max_c1, min_c2, max_c2, use_offset, old_chi):
        c1_cells = 10
        c2_cells = 10
        if old_chi < float('inf') - 1:  # second iteration
            c1_cells = 5
            c2_cells = 5

        delta_c1 = (max_c1 - min_c1) / c1_cells
        delta_c2 = (max_c2 - min_c2) / c2_cells

        last_c1 = False
        last_c2 = False
        if delta_c1 < 0.0001:
            c1_cells = 1
            delta_c1 = max_c1 - min_c1
            last_c1 = True
        if delta_c2 < 0.001:
            c2_cells = 1
            delta_c2 = max_c2 - min_c2
            last_c2 = True

        best_c1 = 1.0
        best_c2 = 0.0
        best_chi = old_chi
        best_set = False

        c1 = min_c1
        for _ in range(c1_cells + 1):
            c2 = min_c2
            for _ in range(c2_cells + 1):
                partial_profile.sum_partial_profiles(c1, c2)
                curr_chi, fit_profile = self.compute_score(partial_profile, use_offset)
                if not best_set or curr_chi < best_chi:
                    best_set = True
                    best_chi = curr_chi
                    best_c1 = c1
                    best_c2 = c2
                c2 += delta_c2
            c1 += delta_c1

        if abs(best_chi - old_chi) > 0.0001 and not (last_c1 and last_c2):
            min_c1 = max(best_c1 - delta_c1, min_c1)
            max_c1 = min(best_c1 + delta_c1, max_c1)
            min_c2 = max(best_c2 - delta_c2, min_c2)
            max_c2 = min(best_c2 + delta_c2, max_c2)
            return self.search_fit_parameters(partial_profile, min_c1, max_c1, min_c2, max_c2, use_offset, best_chi)
        return FitParameters(best_chi, best_c1, best_c2)

    def fit_profile(self, partial_profile, min_c1, max_c1, min_c2, max_c2,
        use_offset, fit_file_name):
        # Compute chi value for default c1/c2
        default_c1 = 1.0
        default_c2 = 0.0
        partial_profile.sum_partial_profiles(default_c1, default_c2)
        default_chi, fit_profile = self.compute_score(partial_profile, use_offset)

        fp = self.search_fit_parameters(partial_profile, min_c1, max_c1,
            min_c2, max_c2, use_offset, float('inf'))
        best_c1 = fp.c1
        best_c2 = fp.c2
        fp.default_chi_square = default_chi
        # Compute a profile for the best c1/c2 combination
        partial_profile.sum_partial_profiles(best_c1, best_c2)
        score, fit_profile = self.compute_score(partial_profile, use_offset, fit_file_name)
        return fp, fit_profile

    def compute_score(self, model_profile, use_offset, fit_file_name=""):
        resampled_profile = Profile(
            qmin=self.exp_profile_.min_q_,
            qmax=self.exp_profile_.max_q_,
            delta=self.exp_profile_.delta_q_,
            constructor=0
        )
        # model_profile and resampled_profile might be different than the C++ version
        model_profile.resample(self.exp_profile_, resampled_profile)
        score = self.scoring_function_.compute_score(self.exp_profile_,
            resampled_profile, use_offset)

        offset = 0.0
        if use_offset:
            offset = self.scoring_function_.compute_offset(self.exp_profile_,
                resampled_profile)
        c = self.scoring_function_.compute_scale_factor(self.exp_profile_,
            resampled_profile, offset)

        if len(fit_file_name) > 0:
            self.write_SAXS_fit_file(fit_file_name, resampled_profile, score,
                c, offset)

        for i in range(resampled_profile.size()):
            resampled_profile.intensity_[i] = resampled_profile.intensity_[i]* c - offset
        return score, resampled_profile

    def write_SAXS_fit_file(self, file_name, model_profile, score, c, offset):
        with open(file_name, "w") as out_file:
            profile_size = min(model_profile.size(), self.exp_profile_.size())
            # header line
            out_file.write("# SAXS profile: number of points = " + str(profile_size) +
                        ", q_min = " + str(self.exp_profile_.min_q_) +
                        ", q_max = " + str(self.exp_profile_.max_q_) +
                        ", delta_q = " + str(self.exp_profile_.delta_q_) + "\n")

            out_file.write("# offset = " + str(offset) +
                        ", scaling c = " + str(c) +
                        ", Chi^2 = " + str(score) + "\n")
            out_file.write("#  q       exp_intensity   model_intensity error\n")

            # Main data
            for i in range(profile_size):
                out_file.write("{:<10.8f} ".format(self.exp_profile_.q_[i]))
                out_file.write("{:<15.8f} ".format(self.exp_profile_.intensity_[i]))
                out_file.write("{:<15.8f} ".format(model_profile.intensity_[i] * c - offset))
                out_file.write("{:<10.8f}\n".format(self.exp_profile_.error_[i]))

        # Write to the second file with .fit extension
        file_name2 = file_name[:-4] + ".fit"
        with open(file_name2, "w") as out_file2:
            # header line
            out_file2.write("# SAXS profile: number of points = " + str(profile_size) +
                            ", q_min = " + str(self.exp_profile_.min_q_) +
                            ", q_max = " + str(self.exp_profile_.max_q_) +
                            ", delta_q = " + str(self.exp_profile_.delta_q_) + "\n")

            out_file2.write("# offset = " + str(offset) +
                            ", scaling c = " + str(c) +
                            ", Chi^2 = " + str(score) + "\n")
            out_file2.write("#  q       exp_intensity   error model_intensity\n")

            # Main data
            for i in range(profile_size):
                out_file2.write("{:<10.8f} ".format(self.exp_profile_.q_[i]))
                out_file2.write("{:<15.8f} ".format(self.exp_profile_.intensity_[i]))
                out_file2.write("{:<10.8f} ".format(self.exp_profile_.error_[i]))
                out_file2.write("{:<15.8f}\n".format(model_profile.intensity_[i] * c - offset))
