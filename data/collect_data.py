#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2016 Yi Cao <ycao16@uw.edu>
#
# Distributed under terms of the GNU General Public License 3.0 license.

"""
This program collects all data for the iPTF16abc paper into a single HDF5
file.
"""


import numpy as np
import os
import tables
from astropy.io import ascii
from astropy.time import Time

# global const
T_MAX = 57499.65


# define photometry columns
class Photometry(tables.IsDescription):
    time = tables.Float64Col()
    flux = tables.Float32Col()
    flux_err = tables.Float32Col()
    mag = tables.Float32Col()
    mag_err = tables.Float32Col()
    telescope = tables.StringCol(16)
    instrument = tables.StringCol(16)
    filter_name = tables.StringCol(16)


# define transmission curve columns
class TransmissionCurve(tables.IsDescription):
    wavelength = tables.Float32Col()
    transmission = tables.Float32Col()


# supernova data set
class SupernovaData(object):
    def __init__(self, filename, clobber=False, title=None):
        if os.path.isfile(filename) and not clobber:
            msg = ("%s: file exists. Change the filename "
                   "or set clobber=True") % filename
            raise FileError(msg)
        self._filename = filename
        self._fp = tables.open_file(filename, mode="a", title=title)

        self._phot_group = self._fp.create_group(
            self._fp.root, "photometry", "Photometric Data")
        self._lc_table = self._fp.create_table(
            self._phot_group, "lc", Photometry, "Light Curves")
        self._lc_table.attrs.FIELD_0_UNIT = "day"
        self._lc_table.attrs.FIELD_1_UNIT = "Jy"
        self._lc_table.attrs.FIELD_2_UNIT = "Jy"
        self._lc_table.attrs.FIELD_3_UNIT = "AB mag"
        self._lc_table.attrs.FIELD_4_UNIT = "AB mag"

        self._filters_group = self._fp.create_filters(
            self._phot_group, "filters", "Filter Transmission Curves")

        self._spec_group = self._fp.create_group(
            self._fp.root, "spectroscopy", "Spectroscopic Data")

    def load_photometry(self, time, filter_name, flux, flux_err,
                        mag, mag_err, telescope, instrument):
        row = self._lc_tables.row
        row["time"] = time
        row["filter_name"] = filter_name
        row["flux"] = flux
        row["flux_err"] = flux_err
        row["mag"] = mag
        row["mag_err"] = mag_err
        row["telescope"] = telescope
        row["instrument"] = instrument
        row.append()


class FileError(object):

    def __init__(self, msg):
        self.msg = msg

    def __expr__(self):
        return self.msg


def load_photometry_data(DataObject):
    # P48
    data_dict = {"g": "lc/forcepsffitdiff_d100151_f1_c11.out",
                 "R": "lc/forcepsffitdiff_d3381_f2_c10.out"}
    for filter_name, filename in data_dict.iteritems():
        data = ascii.read(filename, format="ipac")
        for item in data:
            if item["MJD"] < 57470:
                continue
            if item["flux"] >= 5. * item["sigflux"]:
                DataObject.load_photometry(
                    item["MJD"],
                    filter_name,
                    item["flux"] * 10**(-item["zpmag"] / 2.5) * 3631,
                    abs(item["flux"] * 10**(-item["zpmag"] / 2.5) * 3631) *
                    ((item["sigflux"] / item["flux"])**2 +
                     (item["zprms"] / 2.5)**2)**0.5,
                    -2.5 * np.log10(item["flux"]) + item["zpmag"],
                    ((2.5 * item["sigflux"] / item["flux"])**2 +
                     item["zprms"]**2)**0.5,
                    "P48",
                    "CFH12K")
            else:
                DataObject.load_photometry(
                    item["MJD"],
                    filter_name,
                    item["flux"] * 10**(-item["zpmag"] / 2.5) * 3631,
                    abs(item["flux"] * 10**(-item["zpmag"] / 2.5) * 3631) *
                    ((item["sigflux"] / item["flux"])**2 +
                     (item["zprms"] / 2.5)**2)**0.5,
                    -2.5 * np.log10(item["sigflux"] * 5.) + item["zpmag"],
                    99,
                    "P48",
                    "CFH12K")

    # P60
    with open("lc/Marshal_lc.txt", "r") as fp_txt:
        for item in fp_txt:
            items = item.split(",")
            if len(items) < 3:
                continue
            if not items[7].startswith("\"P60"):
                continue
            filter_name = items[2][1:-1]
            mag = float(items[4])
            mag_err = float(items[5])
            if mag > 90.:
                continue
            DataObject.load_photometry(
                float(items[1]) - 2400000.5,
                filter_name,
                10**(-mag/2.5) * 3631,
                0.921 * mag_err * 10**(-mag/2.5) * 3631,
                mag,
                mag_err,
                "P60",
                "SEDM")

    # LCOGT
    with open("lc/iPTF16abc_lcophot.txt", "r") as fp_txt:
        for item in fp_txt:
            items = item[:-1].split()
            MJD = float(items[1]) - 2400000.5
            mag = float(items[2])
            mag_err = float(item[3])
            filter_name = items[5]
            if filter_name == "B":
                flux = 10**(-mag / 2.5) * 4063  # Jy
            elif filter_name == "V":
                flux = 10**(-mag / 2.5) * 3636  # Jy
            else:
                flux = 10**(-mag / 2.5) * 3631
            flux_err = flux * 0.921 * mag_err
            DataObject.load_photometry(
                MJD,
                filter_name,
                flux,
                flux_err,
                mag,
                mag_err,
                "LCOGT-1m",
                "Sinistro")

    # Swift Data
    with open("swift_phot.txt", "r") as fp_txt:
        for item in fp_txt:
            if item.startswith("#"):
                continue
            items = item.split("\t")
            obs_time = Time(items[1], format="iso", scale="utc")
            DataObject.load_photometry(
                obs_time.mjd,
                items[8].strip(),
                float(items[3]) * 1e-3,
                float(items[4]) * 1e-3,
                -2.5 * np.log10(float(items[3]) / 3.631e6),
                1.0857 * float(items[4]) / float(item[3]),
                "Swift",
                "UVOT")

    filter_group = fp.create_group(phot_group, "filters",
                                   "Filter transmission curves")
    for filter_name in ["g", "R"]:
        table = fp.create_table(filter_group,
                                "P48_" + filter_name, TransmissionCurve,
                                "P48 %s filter transmission" % filter_name)
        table.attrs.FIELD_0_UNIT = "Angstrom"
        table.attrs.FIELD_1_UNIT = "photons per Angstrom"
        transmission = np.genfromtxt("../filters/P48/P48_%s.dat" % filter_name,
                                     names=["wavelength", "transmission"])
        row = table.row
        for item in transmission:
            row["wavelength"] = item["wavelength"]
            row["transmission"] = item["transmission"]
            row.append()
    for filter_name in ["g", "r", "i"]:
        table = fp.create_table(filter_group,
                                "P60_" + filter_name, TransmissionCurve,
                                "P60 %s filter transmission" % filter_name)
        table.attrs.FIELD_0_UNIT = "Angstrom"
        table.attrs.FIELD_1_UNIT = "photons per Angstrom"
        transmission = np.genfromtxt(
            "../filters/SEDm/%sband_eff.dat" % filter_name,
            names=["wavelength", "transmission"])
        row = table.row
        for item in transmission:
            row["wavelength"] = item["wavelength"]
            row["transmission"] = item["transmission"]
            row.append()

    fp.close()


def load_spectroscopy_data(data_filename):
    fp = tables.open_file(data_filename, mode="a")
    spec_group = fp.create_group(fp.root, "spectroscopy", "Spectral data")
