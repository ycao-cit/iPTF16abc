#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2017 Yi Cao <ycao16@uw.edu>
#
# Distributed under terms of the GNU General Public License 3.0 license.

"""
Examine spectral data
"""

import matplotlib.pyplot as plt
import numpy as np
import tables


def main():
    fp = tables.open_file("iPTF16abc.h5", mode="r")
    group = fp.root.spectroscopy
    for spec in group:
        plt.figure()
        w = np.array([row["wavelength"] for row in spec])
        f = np.array([row["flux_lambda"] for row in spec])
        plt.plot(w, f)
        plt.title(spec.attrs.TITLE)
    plt.show()

if __name__ == "__main__":
    main()
