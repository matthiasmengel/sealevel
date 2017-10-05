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

""" Minimal way to do sea level projectoons using the sealevel.projection.project
    function.
"""

import os
import numpy as np
import cPickle as pickle
import settings
import dimarray as da
import pandas as pd
reload(settings)
import sealevel.projection as pr
reload(pr)
import sealevel.get_magicc_gmt_data as mag
reload(mag)
import sealevel.contributor_functions as cf
reload(cf)
import sealevel.calib_settings as cs
reload(cs)

realizations = np.arange(settings.nrealizations)
projection_data = {}

for scen in settings.scenarios:

    print "scenario", scen

    projection_data[scen] = {}
    gmt = mag.magicc_gmt[scen]

    for i, contrib_name in enumerate(settings.project_these):

        print "conribution", contrib_name

        calibdata = pd.read_csv(
            os.path.join(settings.calibfolder, contrib_name+".csv"),
            index_col=[0])

        temp_anomaly_year = cs.temp_anomaly_year[contrib_name]
        sl_contributor = cf.contributor_functions[contrib_name]

        proj = np.zeros([len(settings.proj_period), settings.nrealizations])

        for n in realizations:
            slr, gmt_n, obs_choice, params = pr.project(
                gmt, settings.proj_period, calibdata, temp_anomaly_year,
                sl_contributor, n)
            proj[:, n] = slr

        pdata = da.DimArray(proj, axes=[settings.proj_period, realizations],
                            dims=["time", "runnumber"])
        projection_data[scen][contrib_name] = pdata

        fname = "projected_slr_"+scen+"_n"+str(settings.nrealizations)+".nc"
        da.Dataset(projection_data[scen]).write_nc(os.path.join(
            settings.projected_slr_folder,fname))
