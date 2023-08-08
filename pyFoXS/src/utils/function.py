import math
import sys

class SincFunction:
    def __init__(self, max_value, bin_size):
        self.bin_size_ = bin_size
        self.one_over_bin_size_ = 1.0 / bin_size
        self.max_value_ = max_value

    def value2index(self, value):
        return round(value * self.one_over_bin_size_)

    def sinc(self, x):
        index = self.value2index(x)
        x = index * self.bin_size_
        # WARNING: this function might create a difference with the C++ code (precision difference)
        return sinc_pi(x)

def sinc_pi(x):
    if abs(x) >= 3.3 * sys.float_info.epsilon**(1/4):
        return math.sin(x)/x
    else:
        # |x| < (eps*120)^(1/4)
        return 1 - x * x / 6

class ExpFunction:
    def __init__(self, max_value, bin_size):
        self.bin_size_ = bin_size
        self.one_over_bin_size_ = 1.0 / bin_size
        self.max_value_ = max_value

    def value2index(self, value):
        return round(value * self.one_over_bin_size_)

    def exp(self, x):
        index = self.value2index(abs(x))
        y = index * self.bin_size_
        if x < 0.0:
            return 1.0 / math.exp(y)
        return math.exp(y)
