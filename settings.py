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

calibfolder = "data/calibration"

# The emulator of Deconto & Pollard ice sheet contribution is not
# calibrated here. See github.com/matthiasmengel/fast_ant_sid for
# the calibration of ant_dp16.
calibrate_these = [
    "thermexp",
    "gic",
    "gis_smb",
    "gis_sid",
    # "ant_smb",
    # "ant_sid",
    ]

projected_slr_folder = "data/projection"
scenarios = ["rcp26", "rcp45", "rcp85"]
contrib_ids = ["thermexp", "gic", "gis_smb", "gis_sid", "ant_smb", "ant_sid"]

# PNAS 2016 like Antarctic ice sheet contributon
project_these = ["thermexp", "gic", "gis_smb", "gis_sid", "ant_smb",
                    "ant_sid", "ant_dp16"]
# project_these = ["ant_smb",
#                     "ant_sid"]

# Nature Communications 2018 like Antarctic ice sheet contributon
# project_these = ["thermexp", "gic", "gis_smb", "gis_sid", "ant_dp16"]

proj_period = np.arange(1900, 2101, 1)

# the number of monte carlo samples, 10000 used in PNAS paper
nrealizations = 10

# only possible if you have an ensemble of global mean temperature
# projections
probablistic_climate = False

# No edits below this line needed.
# Convert relative to absolute paths, to make them accessible
# across the operating system
project_root = os.path.dirname(os.path.abspath(__file__))
calibfolder = os.path.join(project_root, calibfolder)
projected_slr_folder = os.path.join(project_root, projected_slr_folder)
