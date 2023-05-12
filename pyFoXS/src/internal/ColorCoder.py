"""
\file IMP/foxs/ColorCoder.cpp \brief

Copyright 2007-2022 IMP Inventors. All rights reserved.
"""

class ColorCoder:
    diff_ = 30

    @staticmethod
    def get_color_for_id(id):
        color_modules = [
            [1, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1], [1, 0, 1], [1, 1, 0]
        ]
        mod = id % 6
        i = id // 6
        mult = 255 - i * ColorCoder.diff_
        r = int(color_modules[mod][0] * mult)
        g = int(color_modules[mod][1] * mult)
        b = int(color_modules[mod][2] * mult)
        return r, g, b
