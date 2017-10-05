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


project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
inputdatadir = os.path.join(project_dir, "data/input/")


######## IPCC mean sl contributions ########

# see Chapter 13, Fifth IPCC Report of WG1, Table 13.5

ipccdata = np.loadtxt(
    inputdatadir +
    "ipcc_ar5/slr_contributions_ch13.csv",
    skiprows=1,
    usecols=(
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15))

ipcc_columnames = [
    "RCP3PD_med",
    "RCP3PD_low",
    "RCP3PD_high",
    "RCP45_med",
    "RCP45_low",
    "RCP45_high",
    "RCP60_med",
    "RCP60_low",
    "RCP60_high",
    "RCP85_med",
    "RCP85_low",
    "RCP85_high"]

ipcc_rownames = [
    "thermexp", "gic", "gis_smb", "ant_smb", "gis_sid", "ant_sid",
    "LandWaterStorage", "mean_slr_2081_2100", "Greenland_Sum", "Antarctica_Sum", "Ice-sheet_rapid_dyn",
    "rate_slr_2081_2100", "mean_slr_2046_2065",
    "mean_slr_2100"]


def get_ipcc_range(rcp, contribution):

    lowind = ipcc_columnames.index(rcp + "_low")
    medind = ipcc_columnames.index(rcp + "_med")
    highind = ipcc_columnames.index(rcp + "_high")
    rowind = ipcc_rownames.index(contribution)
    contrib = ipccdata[rowind, :]
    return np.array([contrib[lowind], contrib[medind],
                     contrib[highind]]) * 1e3  # in mm

ipcc_contrib_estimates = {}

for contrib in ["thermexp", "gic", "gis_smb", "ant_smb", "gis_sid", "ant_sid"]:
    ipcc_contrib_estimates[contrib] = {}
    for rcp in ["RCP3PD", "RCP45", "RCP85"]:
        ipcc_contrib_estimates[contrib][rcp] = get_ipcc_range(rcp, contrib)

ipcc_contrib_estimates["gis"] = {}
for rcp in ["RCP3PD", "RCP45", "RCP85"]:
    ipcc_contrib_estimates["gis"][rcp] = (
        ipcc_contrib_estimates["gis_sid"][rcp] + ipcc_contrib_estimates["gis_smb"][rcp])
