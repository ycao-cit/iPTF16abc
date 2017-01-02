#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2016 Yi Cao <ycao16@uw.edu>
#
# Distributed under terms of the GNU General Public License 3.0 license.

"""
create a hdf5 file that contains filer curves
"""


import numpy as np
import tables


# define transmission curve columns
class TransmissionCurve(tables.IsDescription):
    wavelength = tables.Float32Col()
    transmission = tables.Float32Col()


class FilterData(object):

    def __init__(self, filename):
        self._filename = filename
        self._fp = tables.open_file(filename, mode="a",
                                    title="Filter Transmission Curves")

    def add_transmission_curve(self, telescope, instrument, filter_name,
                               wavelength, transmission, clobber=False):
        if len(wavelength) != len(transmission):
            raise InputError("len(wavelength) != len(transmission)")
        tel_group = self._maybe_create_group(self._fp.root, telescope)
        inst_group = self._maybe_create_group(tel_group, instrument)
        if filter_name in inst_group:
            if clobber:
                self._fp.remove_node(inst_group._f_get_child(filter_name))
            else:
                raise InputError("This filter %s/%s/%s already exists. "
                                 "Set clobber=True to overwrite it" %
                                 (telescope, instrument, filter_name))
        tbl = self._fp.create_table(inst_group, filter_name, TransmissionCurve)
        tbl.attrs.FIELD_0_UNIT = "Angstrom"
        tbl.attrs.FIELD_1_UNIT = "Photons per second per Angstrom"
        row = tbl.row
        for w, t in zip(wavelength, transmission):
            row["wavelength"] = w
            row["transmission"] = t
            row.append()

    def _maybe_create_group(self, root, child):
        if child not in root:
            self._fp.create_group(root, child)
        return root._f_get_child(child)


class InputError(Exception):

    def __init__(self, msg):
        self.msg = msg

    def __expr__(self):
        return self.msg


def main():
    file_object = FilterData("filters.h5")
    # P48
    for filter_name in ["g", "R"]:
        curve = np.genfromtxt("filters/P48/P48_%s.dat" % filter_name,
                              names=["wavelength", "transmission"])
        file_object.add_transmission_curve(
            "P48", "CFH12K", filter_name,
            curve["wavelength"], curve["transmission"])
    # P60
    for filter_name in ["g", "r", "i"]:
        curve = np.genfromtxt("filters/SEDm/%sband_eff.dat" % filter_name,
                              names=["wavelength", "transmission"])
        file_object.add_transmission_curve(
            "P60", "SEDm", filter_name,
            curve["wavelength"], curve["transmission"])
    return

if __name__ == "__main__":
    main()
