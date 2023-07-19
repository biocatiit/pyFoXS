"""
\file Distribution.cpp \brief computes

distribution class implementation

Copyright 2007-2022 IMP Inventors. All rights reserved.
"""

import math

class RadialDistributionFunction:
    def __init__(self, bin_size=0.5):
        self.bin_size = bin_size
        self.distribution = [0.0]
        self.max_distance_ = 50.0

    def scale(self, c):
        self.distribution = [value * c for value in self.distribution]

    def add(self, other_rd):
        for i in range(min(len(self.distribution), len(other_rd.distribution))):
            self.distribution[i] += other_rd.distribution[i]

    def normalize(self):
        # Calculate area
        area = sum(self.distribution)

        # Normalize
        self.distribution = [value / area for value in self.distribution]

    def init(self, r):
        self.distribution = [0.0] * math.ceil(r / self.bin_size)

    def add_to_distribution(self, r, value):
        index = self.get_index_from_distance(r)
        if index >= len(self.distribution):
            ext = [0.0] * (index-self.size()+1)
            self.distribution.extend(ext)
        self.distribution[index] += value

    def get_distance_from_index(self, index):
        return index * self.bin_size

    def get_index_from_distance(self, dist):
        return round(dist * 1/self.bin_size)

    def size(self):
        return len(self.distribution)
