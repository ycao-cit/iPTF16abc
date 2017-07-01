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
import tables
from astropy.io import ascii
from astropy.time import Time
from makelc.ptfidelc import ptfide_light_curve as plc

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


class Spectroscopy(tables.IsDescription):
    wavelength = tables.Float32Col()
    flux_lambda = tables.Float32Col()


# supernova data set
class SupernovaData(object):
    def __init__(self, filename):
        self._filename = filename
        self._fp = tables.open_file(filename, mode="a")

        self._phot_table = self._fp.create_table(
            self._fp.root, "photometry", Photometry, "Photometric Data")
        self._phot_table.attrs.FIELD_0_UNIT = "day"
        self._phot_table.attrs.FIELD_1_UNIT = "Jy"
        self._phot_table.attrs.FIELD_2_UNIT = "Jy"
        self._phot_table.attrs.FIELD_3_UNIT = "AB mag"
        self._phot_table.attrs.FIELD_4_UNIT = "AB mag"

        self._spec_group = self._fp.create_group(
            self._fp.root, "spectroscopy", "Spectroscopic Data")

    def add_photometry(self, time, filter_name, flux, flux_err,
                       mag, mag_err, telescope, instrument):
        row = self._phot_table.row
        row["time"] = time
        row["filter_name"] = filter_name
        row["flux"] = flux
        row["flux_err"] = flux_err
        row["mag"] = mag
        row["mag_err"] = mag_err
        row["telescope"] = telescope
        row["instrument"] = instrument
        row.append()

    def add_spec(self, time, telescope, instrument,
                 wavelength, flux):
        if len(wavelength) != len(flux):
            raise InputError("len(wavelength) != "
                             "len(flux)")
        table_name = "%s_%i" % (instrument.replace("-", "_"), int(time))
        tbl = self._fp.create_table(self._spec_group,
                                    table_name,
                                    Spectroscopy,
                                    table_name)
        tbl.attrs.FIELD_0_UNIT = "Angstrom"
        tbl.attrs.FIELD_1_UNIT = "Arbitrary Unit"
        tbl.attrs.OBS_DATE = time
        tbl.attrs.TELESCOPE = telescope
        tbl.attrs.INSTRUMENT = instrument
        row = tbl.row
        for w, f in zip(wavelength, flux):
            row["wavelength"] = w
            row["flux_lambda"] = f
            row.append()


class InputError(Exception):

    def __init__(self, msg):
        self.msg = msg

    def __expr__(self):
        return self.msg


def main():
    file_object = SupernovaData("iPTF16abc.h5")

    # P48
    data_dict = {"g": "lc/forcepsffitdiff_d100151_f1_c11.out",
                 "R": "lc/forcepsffitdiff_d3381_f2_c10.out"}
    for filter_name, filename in data_dict.items():
        hjdDet, magDet, magUncDet, hjdLim, magLim, hjdFlux, flux, fluxUnc = plc(filename, 2457470)
        for hjd, f, f_unc in zip(hjdFlux, flux, fluxUnc):
            if hjd < 2457470:
                continue
            if hjd in hjdDet:
                file_object.add_photometry(
                    hjd - 2400000.5,
                    filter_name,
                    f,
                    f_unc,
                    magDet[hjdDet == hjd],
                    magUncDet[hjdDet == hjd],
                    "P48",
                    "CFH12K")
            elif hjd in hjdLim:
                file_object.add_photometry(
                    hjd - 2400000.5,
                    filter_name,
                    f,
                    f_unc,
                    magLim[hjdLim == hjd],
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
            file_object.add_photometry(
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
            mag_err = float(items[3])
            filter_name = items[5]
            if filter_name == "B":
                flux = 10**(-mag / 2.5) * 4063  # Jy
            elif filter_name == "V":
                flux = 10**(-mag / 2.5) * 3636  # Jy
            else:
                flux = 10**(-mag / 2.5) * 3631
            flux_err = flux * 0.921 * mag_err
            file_object.add_photometry(
                MJD,
                filter_name,
                flux,
                flux_err,
                mag,
                mag_err,
                "LCO-1m",
                "Sinistro")

    # Swift Data
    with open("lc/swift_phot.txt", "r") as fp_txt:
        for item in fp_txt:
            if item.startswith("#"):
                continue
            items = item.split("\t")
            obs_time = Time(items[1], format="iso", scale="utc")
            if float(items[3]) < 5. * float(items[4]):
                continue
            file_object.add_photometry(
                obs_time.mjd,
                items[8].strip(),
                float(items[3]) * 1e-3,
                float(items[4]) * 1e-3,
                -2.5 * np.log10(float(items[3]) / 3.631e6),
                1.0857 * float(items[4]) / float(items[3]),
                "Swift",
                "UVOT")
    # RATIR Data
    with open("lc/16abc_ratir_clean.dat", "r") as fp_txt:
        for item in fp_txt:
            if item.startswith("#"):
                continue
            items = item.split("  ")
            MJD = float(items[0])
            filter_name = items[1][-1]
            mag = 25 - 2.5*np.log10(float(items[2]))
            mag_err = 2.5/np.log(10)*float(items[3])/float(items[2])
            flux = 10**(-0.4*mag) * 3631 # Jy
            flux_err = flux*float(items[3])/float(items[2])
            file_object.add_photometry(
                MJD,
                filter_name,
                flux,
                flux_err,
                mag,
                mag_err,
                "SPM-1.5m",
                "RATIR")

    # Spectra
    filename = "spec/16abc_20160405_DCT_v1.ascii"
    curve = np.genfromtxt(filename,
                          names=["wavelength", "flux"])
    idx = np.logical_and(curve["wavelength"] > 3300,
                         curve["wavelength"] < 7500)
    curve = curve[idx]
    file_object.add_spec(57483.26, "DCT", "DeVeny",
                         curve["wavelength"], curve["flux"])

    filename = "spec/16abc_20160405_Gemini_N_v1.ascii"
    curve = np.genfromtxt(filename,
                          names=["wavelength", "flux"])
    file_object.add_spec(57483.88, "Gemini-North", "GMOS",
                         curve["wavelength"], curve["flux"])

    filename = "spec/16abc_20160406_Keck2_v1.ascii"
    curve = np.genfromtxt(filename,
                          names=["wavelength", "flux"])
    idx = np.logical_and(curve["wavelength"] > 5500,
                         curve["wavelength"] < 8100)
    curve = curve[idx]
    file_object.add_spec(57484.51, "Keck-II", "DEIMOS",
                         curve["wavelength"], curve["flux"])

    filename = "spec/16abc_20160408_Keck2_v2.ascii"
    curve = np.genfromtxt(filename,
                          names=["wavelength", "flux"])
    idx = np.logical_and(curve["wavelength"] > 5500,
                         curve["wavelength"] < 8100)
    curve = curve[idx]
    file_object.add_spec(57486.51, "Keck-II", "DEIMOS",
                         curve["wavelength"], curve["flux"])

    filename = "spec/16abc_20160410_Keck1_v1.ascii"
    curve = np.genfromtxt(filename,
                          names=["wavelength", "flux"])
    file_object.add_spec(57488.38, "Keck-I", "LRIS",
                         curve["wavelength"], curve["flux"])

    filename = "spec/16abc_20160411_FTN_v1.ascii"
    curve = np.genfromtxt(filename,
                          names=["wavelength", "flux"])
    idx = np.logical_and(curve["wavelength"] > 3300,
                         curve["wavelength"] < 9000)
    curve = curve[idx]
    file_object.add_spec(57489.506,
                         "LCOGT-2m", "FLOYDS",
                         curve["wavelength"], curve["flux"])

    filename = "spec/16abc_20160412_FTN_v1.ascii"
    curve = np.genfromtxt(filename,
                          names=["wavelength", "flux"])
    idx = np.logical_and(curve["wavelength"] > 3300,
                         curve["wavelength"] < 10000)
    curve = curve[idx]
    file_object.add_spec(57490.396,
                         "LCOGT-2m", "FLOYDS",
                         curve["wavelength"], curve["flux"])

    filename = "spec/16abc_20160413_FTS_v1.ascii"
    curve = np.genfromtxt(filename,
                          names=["wavelength", "flux"])
    idx = np.logical_and(curve["wavelength"] > 3300,
                         curve["wavelength"] < 10000)
    curve = curve[idx]
    file_object.add_spec(57491.551,
                         "LCOGT-2m", "FLOYDS",
                         curve["wavelength"], curve["flux"])

    filename = "spec/16abc_20160414_VLT_v1.ascii"
    curve = np.genfromtxt(filename,
                          names=["wavelength", "flux"])
    idx = np.logical_and(curve["wavelength"] > 3300,
                         curve["wavelength"] < 30000)
    curve = curve[idx]
    file_object.add_spec(57492.20,
                         "VLT", "X-shooter",
                         curve["wavelength"], curve["flux"])

    filename = "spec/16abc_20160416_VLT_v1.ascii"
    curve = np.genfromtxt(filename,
                          names=["wavelength", "flux"])
    file_object.add_spec(57494.00,
                         "VLT", "UVES",
                         curve["wavelength"], curve["flux"])

    filename = "spec/16abc_20160425_FTN_v1.ascii"
    curve = np.genfromtxt(filename,
                          names=["wavelength", "flux"])
    idx = np.logical_and(curve["wavelength"] > 3300,
                         curve["wavelength"] < 10000)
    curve = curve[idx]
    file_object.add_spec(57503.319,
                         "LCOGT-2m", "FLOYDS",
                         curve["wavelength"], curve["flux"])

    filename = "spec/16abc_20160428_NOT_v1.ascii"
    curve = np.genfromtxt(filename,
                          names=["wavelength", "flux"])
    idx = np.logical_and(curve["wavelength"] > 3600,
                         curve["wavelength"] < 8100)
    curve = curve[idx]
    file_object.add_spec(Time("2016-04-28", format="iso").mjd,
                         "NOT", "ALFOSC",
                         curve["wavelength"], curve["flux"])

    filename = "spec/16abc_20160430_FTN_v1.ascii"
    curve = np.genfromtxt(filename,
                          names=["wavelength", "flux"])
    idx = np.logical_and(curve["wavelength"] > 3300,
                         curve["wavelength"] < 10000)
    curve = curve[idx]
    file_object.add_spec(57508.272,
                         "LCOGT-2m", "FLOYDS",
                         curve["wavelength"], curve["flux"])

    filename = "spec/16abc_20160510_Keck1_v2.ascii"
    curve = np.genfromtxt(filename,
                          names=["wavelength", "flux"])
    file_object.add_spec(57518.42, "Keck-I", "LRIS",
                         curve["wavelength"], curve["flux"])

    filename = "spec/16abc_20160512_VLT_v1.ascii"
    curve = np.genfromtxt(filename,
                          names=["wavelength", "flux"])
    idx = np.logical_and(curve["wavelength"] > 3300,
                         curve["wavelength"] < 30000)
    curve = curve[idx]
    file_object.add_spec(57520.03,
                         "VLT", "X-shooter",
                         curve["wavelength"], curve["flux"])

    filename = "spec/16abc_20160521_FTN_v1.ascii"
    curve = np.genfromtxt(filename,
                          names=["wavelength", "flux"])
    idx = np.logical_and(curve["wavelength"] > 4000,
                         curve["wavelength"] < 9000)
    curve = curve[idx]
    file_object.add_spec(57529.405,
                         "LCOGT-2m", "FLOYDS",
                         curve["wavelength"], curve["flux"])

    filename = "spec/16abc_20160603_FTN_v1.ascii"
    curve = np.genfromtxt(filename,
                          names=["wavelength", "flux"])
    idx = np.logical_and(curve["wavelength"] > 4000,
                         curve["wavelength"] < 9000)
    curve = curve[idx]
    file_object.add_spec(57542.408,
                         "LCOGT-2m", "FLOYDS",
                         curve["wavelength"], curve["flux"])

    filename = "spec/16abc_20160611_FTN_v1.ascii"
    curve = np.genfromtxt(filename,
                          names=["wavelength", "flux"])
    idx = np.logical_and(curve["wavelength"] > 4000,
                         curve["wavelength"] < 9000)
    curve = curve[idx]
    file_object.add_spec(57550.402,
                         "LCOGT-2m", "FLOYDS",
                         curve["wavelength"], curve["flux"])

    filename = "spec/16abc_20160623_FTN_v1.ascii"
    curve = np.genfromtxt(filename,
                          names=["wavelength", "flux"])
    file_object.add_spec(57562.379,
                         "LCOGT-2m", "FLOYDS",
                         curve["wavelength"], curve["flux"])

    return

if __name__ == "__main__":
    main()
