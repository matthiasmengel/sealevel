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

""" Code for Fig. S1 as in
    M. Mengel et al.
    Future sea-level rise constrained by observations and long-term commitment
    PNAS (2016)
    (C) Matthias Mengel working at Potsdam Institute for Climate Impact Research

    Note: You may not be able to fully plot this figure as it contains data that is not
          openly available. Request data from authors or comment out.
"""

import os
import glob
import sys
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib import cm
import itertools
import cPickle as pickle
from scipy.io import loadmat
lib_path = os.path.abspath('../src')
sys.path.append(lib_path)
import get_calibration_data as gd
reload(gd)
import calib_settings as cs
reload(cs)
import contributor_functions as cf
reload(cf)
import dimarray as da
import collections

plt.rcParams['xtick.major.pad'] = 10
plt.rcParams['font.size'] = 12
plt.rcParams['lines.markeredgewidth'] = 2
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['figure.figsize'] = 8, 8
plt.rcParams['figure.facecolor'] = "white"
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['pdf.fonttype'] = '42'

# contrib_upper700 = collections.OrderedDict([
# ("domingues08",gd.thermo_obs_domingues08),
# ("ishii09",gd.thermo_obs_ishii09),
# ("levitus12",gd.thermo_obs_levit12_700),
# # ("church11",gd.church_observed_700m)]
# ])

# contrib_700_2000m = collections.OrderedDict([
# ("levitus12",gd.thermo_obs_levit12_2000 - gd.thermo_obs_levit12_700),
# ("church11",gd.church_observed_700m_2000m)])

# contrib_below2000m = collections.OrderedDict([
# ("zero",zero),
# ("purkey10",gd.purkey10_below2000m)])


fig = plt.figure(6)
plt.subplots_adjust(bottom=0.1, top=0.95)
plt.clf()
ax = plt.subplot(111)


observations = cs.contrib_upper700
# for obs in cs.contrib_upper700:
cols = [cm.Accent(np.float(k) / 6) for k in np.arange(6)]


# for contrib in plot_these:
#     # def get_obs_and_calibrated(contrib):

#     calibdata = pickle.load(open("../calibrationdata/"+contrib+".pkl","rb"))
#     observations = cs.__dict__[contrib+"_observations"]
#     cols = [cm.Accent(np.float(k)/len(observations)) for k in np.arange(len(observations))]

#     fig = plt.figure(5)
#     plt.subplots_adjust(bottom=0.06,top=0.95)
#     plt.clf()
#     ax = plt.subplot(111)

def rl(ts):
    return ts - ts[1986:2005].mean()

of = 100.
i = 0
offset = i / of
ax.plot(gd.purkey10_below2000m.time, rl(gd.purkey10_below2000m) + offset,
        label="purkey10 below 2000m", lw=2,
        color=cols[i])
i += 1
for obs in cs.contrib_700_2000m:
    offset = i / of
    ax.plot(cs.contrib_700_2000m[obs].time, rl(cs.contrib_700_2000m[obs]) + offset,
            label=obs + " 700m-2000m", lw=2,
            color=cols[i])
    i += 1
for obs in cs.contrib_upper700:
    offset = i / of
    ax.plot(cs.contrib_upper700[obs].time, rl(cs.contrib_upper700[obs]) + offset,
            label=obs + " 0m-700m", lw=2,
            color=cols[i])
    i += 1


# ncol = 2 if contrib == "thermexp" else 1
l2 = ax.legend(loc="upper left", ncol=2)
l2.draw_frame(0)

# if contrib == "thermexp":
ax.set_ylim(-0.008, 0.062)

ax.set_xlabel("Time in years")
ax.set_ylabel("Sea-level contribution in m")

# ax.set_xlim(1880,2015)

plt.draw()
plt.show()

runname = __file__.split("/")[-1][0:-3]
plt.savefig("../figures/" + runname + ".pdf")
