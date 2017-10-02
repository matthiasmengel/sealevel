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

calibfolder = "data/calibration"
projected_slr_folder = "data/projection"

scenarios = ["RCP3PD", "RCP45", "RCP85"]

contrib_ids = ["thermexp", "gic", "gis_smb", "gis_sid", "ant_smb", "ant_sid"]
project_these = ["thermexp", "gic", "gis_smb", "gis_sid", "ant_smb", "ant_sid"]

calibrate_these = [
    "thermexp",
    "gic",
    "gis_smb",
    "gis_sid",
    "ant_smb",
    "ant_sid"]

# TODO: finde a better name for classname
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

