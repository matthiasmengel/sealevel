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

""" get the MAGICC v 6.0 600 member global mean temperature projection ensemble,
    see Meinshausen et al., Nature (2009) for details.
"""

import os
import numpy as np
import dimarray as da
import get_gmt_data as ggd
reload(ggd)

project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
inputdatadir = os.path.join(project_dir, "data/input/")


# Filter out invalid runs.
invalid_runs = {
    "RCP3PD": [76, 148, 198, 199, 333, 369],
    "RCP45": [148, 198, 333],
    "RCP60": [146, 147, 148, 149, 197, 198, 199, 211, 282, 333],
    "RCP85": [3, 36, 54, 56, 61, 75, 76, 87, 101, 123, 136, 146, 147, 148, 149, 150, 151, 153, 165, 197, 198, 199, 204, 211, 212, 239, 258, 282, 289, 291, 302, 313, 321, 332, 333, 347, 357, 375, 401, 449, 462, 468, 483, 510, 514, 540, 543, 544, 551, 560, 565, 581, 588]
}

magicc_gmt = {}
magicc_scen = ["RCP3PD", "RCP45", "RCP60", "RCP85"]

for scen in magicc_scen:

    magicc_runs = np.loadtxt(
        inputdatadir +
        '/RcpHistoricallyConstrained/RCP.2500/' +
        scen +
        '_histConstrain4Katja_DAT_SURFACE_TEMP.txt')
    magicc_years = np.array(magicc_runs[:, 0], dtype="int")
    magicc_runs = magicc_runs[:, 1:]
    magicc_runs = np.delete(magicc_runs, invalid_runs[scen], axis=1)
    magicc_gmt_scen = da.DimArray(
        magicc_runs, dims=[
            "time", "runnumber"], axes=[
            magicc_years, np.arange(
                magicc_runs.shape[1])])
    # harmonize along 1950-1980 mean, make them equvalent to giss
    magicc_gmt[scen] = magicc_gmt_scen - magicc_gmt_scen[1951:1980,
                                                         :].mean(axis=0) + ggd.preind_to_1951_1980
