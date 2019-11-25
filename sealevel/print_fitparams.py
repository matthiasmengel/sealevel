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

import numpy as np
import pickle as pickle

contrib_ids = ["thermexp", "gic", "gis_sid", "gis_smb", "ant_sid", "ant_smb"]
rcpcoldict = {"RCP3PD": "#3b4ba8", "RCP45": "yellow", "RCP85": "#ee1d23"}
labels = {"gic": "Mountain glaciers", "thermexp": "Thermal expansion", "gis_sid": "Greenland SID",
          "gis_smb": "Greenland SMB", "ant_sid": "Antarctica SID", "ant_smb": "Antarctica SMB"}


def rd(numb):
    return '%s' % float('%.3g' % numb)


def minmax(array):
    return rd(np.min(array)) + " to " + rd(np.max(array))

for name in contrib_ids:
    print(labels[name], end=' ')
    calibdata = pickle.load(open("../data/calibration/" + name + ".pkl", "rb"))
    calibrated_param = np.array([])
    for obs in calibdata["params"]:
        # print obs,
        cal_par_obs = calibdata["params"][obs].calibrated_parameter
        calibrated_param = np.append(calibrated_param, cal_par_obs)

    print(minmax(calibdata["params"][obs].commitment_parameter), end=' ')
    print(minmax(calibrated_param))
