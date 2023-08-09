"""
\file Distribution.cpp \brief computes

distribution class implementation

Copyright 2007-2022 IMP Inventors. All rights reserved.
"""

import numpy as np
from numba import jit

class RadialDistributionFunction:
    def __init__(self, bin_size=0.5):
        self.bin_size = bin_size
        self.distribution = np.array([0.0])
        self.max_distance_ = 50.0

    def scale(self, c):
        self.distribution *= c # [value * c for value in self.distribution]

    def add(self, other_rd):
        for i in range(min(self.distribution.shape[0], other_rd.distribution.shape[0])):
            self.distribution[i] += other_rd.distribution[i]

    def normalize(self):
        # Calculate area
        area = np.sum(self.distribution)

        # Normalize
        self.distribution /= area # [value / area for value in self.distribution]

    def add_to_distribution(self, r, value):
        index = get_index_from_distance(self.bin_size, r)
        if index >= self.distribution.shape[0]:
            ext = np.zeros(index-self.size()+1) # [0.0] * (index-self.size()+1)
            self.distribution = np.concatenate((self.distribution, ext))
            self.max_distance_ = get_distance_from_index(self.bin_size, index + 1)
        self.distribution[index] += value

    def size(self):
        return self.distribution.shape[0]

@jit(target_backend='cuda', nopython=True)
def get_index_from_distance(bin_size, dist):
        return round(dist * 1/bin_size)

@jit(target_backend='cuda', nopython=True)
def get_distance_from_index(bin_size, index):
    return index * bin_size
