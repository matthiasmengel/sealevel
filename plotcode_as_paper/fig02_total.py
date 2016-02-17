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


""" Code for Fig. 2 as in
    M. Mengel et al.
    Future sea-level rise constrained by observations and long-term commitment
    PNAS (2016)
    (C) Matthias Mengel working at Potsdam Institute for Climate Impact Research

    Note: You may not be able to fully plot this figure as it contains data that is not
    openly available. Request data from authors or comment out.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

lib_path = os.path.abspath('../src')
sys.path.append(lib_path)
import get_calibration_data as gd
reload(gd)
import get_gmt_data as ggd
reload(ggd)
import sealevel as sl
reload(sl)
import calib_settings as cs
reload(cs)
import dimarray as da
import cPickle as pickle
import matplotlib.font_manager as font_manager

path = '/Users/mengel/Library/Fonts/OpenSans-Bold.ttf'
prop = font_manager.FontProperties(fname=path)
plt.rcParams['font.family'] = prop.get_name()
plt.rcParams['font.sans-serif'] = prop.get_name()

plt.rcParams['xtick.major.pad'] = 10
plt.rcParams['font.size'] = 12
plt.rcParams['lines.markeredgewidth'] = 2
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['figure.figsize'] = 8, 8
plt.rcParams['figure.facecolor'] = "white"
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['pdf.fonttype'] = '42'


contrib_ids = ["thermexp", "gic", "gis_smb", "gis_sid", "ant_smb", "ant_sid"]

nrealizations = 10000
realizations = np.arange(nrealizations)

all_contributions = {}

# Monte Carlo Sampling from all the different observational datasets and
# fitting/calibrated parametes.
for i, contrib_name in enumerate(contrib_ids):

    print "conribution", contrib_name
    calibdata = pickle.load(
        open(
            "../data/calibration/" +
            contrib_name +
            ".pkl",
            "rb"))

    proj = np.zeros([nrealizations,len(ggd.giss_temp.time)])

    for n in realizations:
        proj[n,:] = sl.project(ggd.giss_temp, ggd.giss_temp.time, calibdata, n)

    pdata = da.DimArray(proj, axes=[realizations, ggd.giss_temp.time],
                        dims=["runnumber","time"])
    all_contributions[contrib_name] = pdata

all_contributions = da.DimArray(all_contributions, dims=[
                                "contribution", "runnumber", "time"])

obs_period = np.arange(1900, 2009)

fig = plt.figure(2)
plt.clf()
ax = plt.subplot(111)

natural_gic_contrib = gd.marzeion_gic_nat_up.mean(axis="model") / 1.e3
# gd.
past_sl_anom = gd.church_past_sl[obs_period] - natural_gic_contrib[obs_period]
past_sl_anom -= past_sl_anom[1986:2005].mean()
lower_bound = past_sl_anom - gd.church_past_sl_std[obs_period]
upper_bound = past_sl_anom + gd.church_past_sl_std[obs_period]

#### observed ####
ax.plot(obs_period, past_sl_anom * 1e3, lw=2, color="black", alpha=1.,
        label="observed\nChurch et al.", marker="|", markevery=10, markersize=10)

ax.fill_between(obs_period, lower_bound * 1e3, upper_bound *
                1e3, lw=.5, color="black", alpha=.3)

hay15_past_sl = gd.hay15_past_sl - natural_gic_contrib[obs_period]
hay15_past_sl -= hay15_past_sl[1986:2005].mean()
ax.plot(hay15_past_sl.time, hay15_past_sl * 1e3, lw=2, color="blue", alpha=1.,
        label="observed\nHay et al.", marker="x", markevery=10, markeredgewidth=2, markersize=5)

hay15_lower = gd.hay15_lower - \
    natural_gic_contrib[obs_period] - gd.hay15_past_sl[1986:2005].mean()
hay15_upper = gd.hay15_upper - \
    natural_gic_contrib[obs_period] - gd.hay15_past_sl[1986:2005].mean()
ax.fill_between(hay15_lower.time, hay15_lower * 1e3,
                hay15_upper * 1e3, lw=.5, color="blue", alpha=.2)

#### calculated from GMT ####
tslr = all_contributions.sum(axis="contribution")
tslr -= tslr[:, 1986:2005].mean(axis="time")

upper_perc = da.DimArray(np.percentile(tslr, 95, axis=0),
                         axes=tslr.time, dims="time")
lower_perc = da.DimArray(np.percentile(tslr, 5, axis=0),
                         axes=tslr.time, dims="time")
median = da.DimArray(np.percentile(tslr, 50, axis=0),
                     axes=tslr.time, dims="time")


ax.plot(tslr.time, median * 1e3, lw=2, color="#c98345",
        label="reconstructed", alpha=1., )
ax.fill_between(lower_perc.time, lower_perc * 1e3,
                upper_perc * 1e3, lw=.5, color="#c98345", alpha=.3)

ax.set_xlim(obs_period[0], obs_period[-1] + 1)
ax.set_xlabel("Time in years")
ax.set_ylabel("Sea level in mm")
# axb.set_ylabel("gmt")
l1 = ax.legend(ncol=1, loc="upper left")
l1.draw_frame(0)
for l in l1.get_lines():
    l.set_alpha(1)

plt.draw()
plt.show()

runname = __file__.split("/")[-1][0:-3]
figname = "../figures/" + runname + ".pdf"

plt.savefig(figname)
