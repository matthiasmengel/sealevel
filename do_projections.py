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


import os
import numpy as np
import cPickle as pickle
import settings
reload(settings)
import src.sealevel as sl
reload(sl)
import src.get_magicc_gmt_data as mag
reload(mag)
import src.contributor_functions as cf
reload(cf)

realizations = np.arange(settings.nrealizations)
projection_data = {}

for scen in settings.scenarios:

    print "scenario", scen

    projection_data[scen] = {}
    gmt = mag.magicc_gmt[scen]

    for i, contrib_name in enumerate(settings.project_these):

        print "conribution", contrib_name
        calibdata = pickle.load( open(
                 os.path.join(settings.calibfolder, contrib_name+".pkl"),
                    "rb" ) )

        proj = np.zeros([len(proj_period), nrealizations])

        for n in realizations:
            slr,driving_temp = sl.project(gmt, proj_period, calibdata, n)
            proj[:, n] = slr

        pdata = da.DimArray(proj, axes=[proj_period, realizations],
                            dims=["time", "runnumber"])
        projection_data[scen][contrib_name] = pdata

fname = os.path.join(settings.projected_slr_folder,
                     "projected_slr_"+str(nrealizations)+"samples.pkl"
print "save to pickle."
pickle.dump(projection_data, open(fname, "wb"), protocol=2)
