# This file is part of SEALEVEL - a tool to estimates future sea-level rise
# constrained by past obervations and long-term sea-level commitment
# Copyright (C) 2016 Matthias Mengel working at PIK Potsdam

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# LICENSE.txt for more details.


""" matthias.mengel@pik
"""

import os
import glob
import sys
import numpy as np
# import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib import cm
import itertools
from scipy.io import loadmat
lib_path = os.path.abspath('../src')
sys.path.append(lib_path)
import dimarray as da
import optparse
from mpl_toolkits.axes_grid1 import make_axes_locatable
import timeit
import pickle as pickle

# make use of -n or -g flags with this script
parser = optparse.OptionParser()
parser.add_option('-n', '--realizations', type="int",
                  help='get number of realizations')
parser.add_option(
    '-g',
    help='get new data',
    dest='get_data',
    default=False,
    action='store_true')
(opts, args) = parser.parse_args()

contrib_ids = ["thermexp", "gic", "gis_sid", "gis_smb", "ant_sid", "ant_smb"]
rcpcoldict = {"RCP3PD": "#3b4ba8", "RCP45": "yellow", "RCP85": "#ee1d23"}
labels = {"gic": "Mountain glaciers", "thermexp": "Thermal expansion", "gis_sid": "Greenland SID",
          "gis_smb": "Greenland SMB", "ant_sid": "Antarctica SID", "ant_smb": "Antarctica SMB"}

if opts.get_data:
    realizations = 10000
    realizations = realizations if opts.realizations is None else opts.realizations
    print("use", realizations, "realizations")
    projection_data = pickle.load(
        open(
            "../data/projection/projected_slr_" +
            str(realizations) +
            "samples.pkl",
            "rb"))


def rd(array, percentile, years):

    if isinstance(years, int):
        numb = array[years]
    else:
        numb = array[years].mean(axis="time")

    perc = np.percentile(numb, percentile)
    # in mm
    return '%s' % float('%.3g' % perc)

## switch for the 2081-2100 mean
printyr = 2100
#printyr = np.arange(2081,2101,1)
print("## years", printyr)


print("RCP3PD", "RCP45", "RCP85")
for i, name in enumerate(contrib_ids):
    print(labels[name], ":", end=' ')
    for k, scen in enumerate(["RCP3PD", "RCP45", "RCP85", ]):

        contrib = projection_data[scen][name] * 1.e3
        # anomaly to 1986-2005
        contrib = contrib - contrib[1986:2005, :].mean(axis="time")

        if scen != "RCP85":
            print(rd(contrib, 50, printyr), "(", rd(contrib, 5, printyr), "to", rd(contrib, 95, printyr), ")", end=' ')
        else:
            print(rd(contrib, 50, printyr), "(", rd(contrib, 5, printyr), "to", rd(contrib, 95, printyr), ")")

#### total slr; sum up contributions first ####


def rdd(numb):
    # in mm
    return '%s' % float('%.4g' % numb)

print("total", ":", end=' ')
for k, scen in enumerate(["RCP3PD", "RCP45", "RCP85", ]):
    total_slr = da.zeros_like(projection_data[scen]["thermexp"])
    for i, name in enumerate(contrib_ids):
        # sum up all contributions
        single_contrib = projection_data[scen][name]
        # if nans occur, clear these contributions
        # single_contrib[np.isnan(single_contrib)] = 0.
        total_slr += single_contrib

    di_total_slr = da.DimArray(total_slr, dims=["time", "runnumber"])
    di_total_slr -= di_total_slr[1986:2005, :].mean(axis="time")
    if isinstance(printyr, int):
        mn = di_total_slr[printyr, :] * 1.e3
    else:
        mn = di_total_slr[printyr, :].mean(axis="time") * 1.e3
    # mn = di_total_slr[2081:2100,:].mean(axis="time")
    low = np.percentile(mn, 5)
    med = np.percentile(mn, 50)
    upp = np.percentile(mn, 95)
    print(rdd(med), "(", rdd(low), "to", rdd(upp), ")", end=' ')
