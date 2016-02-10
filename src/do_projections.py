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
import dimarray as da
import timeit
import cPickle as pickle
import optparse
import multiprocessing
import signal
import sealevel as sl
reload(sl)
import get_magicc_gmt_data as mag
reload(mag)
import contributor_functions as cf
reload(cf)


parser = optparse.OptionParser()
parser.add_option('-n', '--nrealizations', type="int",
                  help='get number of realizations')
(opts, args) = parser.parse_args()


contrib_ids = ["thermexp", "gic", "gis_smb", "gis_sid", "ant_smb", "ant_sid"]
project_these = ["thermexp", "gic", "gis_smb", "gis_sid", "ant_smb", "ant_sid"]

classname = {
    "gic": "glaciers_and_icecaps",
    "thermexp": "thermal_expansion",
    "gis": "surfacemassbalance_gis",
    "gis_sid": "solid_ice_discharge_gis",
    "gis_smb": "surfacemassbalance_gis",
    "ant_sid": "solid_ice_discharge_ais",
    "ant_smb": "surfacemassbalance_ais"}

proj_period = np.arange(1900, 2101, 1)

# the number of monte carlo samples, 10000 used in PNAS paper
nrealizations = 100
nrealizations = nrealizations if opts.nrealizations is None else opts.nrealizations
realizations = np.arange(nrealizations)

projection_data = {}

for scen in ["RCP3PD", "RCP45", "RCP85"]:

    print "scenario", scen

# def do_projection(scen):
    projection_data[scen] = {}
    gmt = mag.magicc_gmt[scen]
    #selected_numbers = np.random.random_integers(0,len(magicc_gmt.runnumber)-1,realizations)
    #magicc_selected  = magicc_gmt[:,selected_numbers]

    for i, contrib_name in enumerate(project_these):

        print "conribution", contrib_name
        calibdata = pickle.load(
            open(
                "../data/calibration/" +
                contrib_name +
                ".pkl",
                "rb"))

        # dd = p.map_async(project,selected_numbers)
        proj = np.zeros([len(proj_period), nrealizations])

        for n in realizations:
            proj[:, n] = sl.project(gmt, proj_period, calibdata, n)

        #proj = map(project,selected_numbers)
        pdata = da.DimArray(proj, axes=[proj_period, realizations],
                            dims=["time", "runnumber"])
        # pdata = pdata.dropna(axis=1)
        projection_data[scen][contrib_name] = pdata

fname = "../data/projection/projected_slr_" + \
    str(nrealizations) + "samples.pkl"
print "save to pickle."
pickle.dump(projection_data, open(fname, "wb"), protocol=2)
