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

""" Code for Fig. 3 as in
    M. Mengel et al.
    Future sea-level rise constrained by observations and long-term commitment
    PNAS (2016)
    (C) Matthias Mengel working at Potsdam Institute for Climate Impact Research

    Note: Please run src/do_projections first, then (in ipython)
          run fig03_rcp_total.py -g -n numberofsamples
"""

import os
import glob
import sys
import numpy as np
# import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib import cm
import itertools
from scipy.io import loadmat
lib_path = os.path.abspath('../src')
sys.path.append(lib_path)
import contributor_functions as cf
reload(cf)
import get_ipcc_data as ipcc
reload(ipcc)
#import get_calibration_data as ggd; reload(gd)
import dimarray as da
import optparse
from mpl_toolkits.axes_grid1 import make_axes_locatable
import timeit
import cPickle as pickle
import matplotlib.font_manager as font_manager

path = '/Users/mengel/Library/Fonts/OpenSans-Bold.ttf'
prop = font_manager.FontProperties(fname=path)
plt.rcParams['font.family'] = prop.get_name()
plt.rcParams['font.sans-serif'] = prop.get_name()

parser = optparse.OptionParser()
parser.add_option('-n', '--realizations',  type="int",
                  help='get number of realizations')
parser.add_option('-g', help='get new data', dest='get_data',
                  default=False, action='store_true')
(opts, args) = parser.parse_args()


plt.rcParams['xtick.major.pad'] = 10
plt.rcParams['font.size'] = 12
plt.rcParams['lines.markeredgewidth'] = 2
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['figure.figsize'] = 12, 10
plt.rcParams['figure.facecolor'] = "white"
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['pdf.fonttype'] = '42'

contrib_ids = ["thermexp", "gic", "gis_smb", "gis_sid", "ant_smb", "ant_sid"]

# rcpcoldict =
# {"RCP3PD":"#3b4ba8","RCP45":"#8cbbe7","RCP60":"#f46f25","RCP85":"#ee1d23"}
rcpcoldict = {"RCP3PD": "#2256A6", "RCP45": "#73B2E1", "RCP85": "#EE322D"}
rcpnamedict = {"RCP3PD": "RCP26", "RCP45": "RCP45", "RCP85": "RCP85"}

if opts.get_data:
    realizations = 10000
    realizations = realizations if opts.realizations == None else opts.realizations
    print "use", realizations, "realizations"
    projection_data = pickle.load(open("../data/projection/projected_slr_" +
                                       str(realizations) + "samples.pkl", "rb"))

plot_period = np.arange(2000, 2101, 1)

plt.figure(4, figsize=(8, 8))
plt.clf()
plt.subplots_adjust(left=0.1, bottom=0.1, right=1.0, top=0.97,
                    wspace=0.1, hspace=None)
ax6 = plt.subplot(111)
# ax6b = ax6.twinx()

# ipcc ranges
divider = make_axes_locatable(ax6)
axy = divider.append_axes("right", size=0.8, pad=0.0, sharey=ax6)
axy.axis("off")
axy.axvspan(0, 10, facecolor='0.5', alpha=0.2, lw=0)
# contributions_file.write("\n## total contribution ##\n")
total_contribs = {}
xloc = 10
oloc = 4

for k, scen in enumerate(["RCP85", "RCP45", "RCP3PD"]):
    total_slr = da.zeros_like(projection_data[scen]["thermexp"])

    for i, name in enumerate(contrib_ids):
        # sum up all contributions
        single_contrib = projection_data[scen][name]
        # if nans occur, clear these contributions
        # single_contrib[np.isnan(single_contrib)] = 0.
        total_slr += single_contrib
        # contrib = projection_data[scen][name].sum(axis=1)

    contrib = total_slr - total_slr[1986:2005,
                                    :].mean(axis=0)  # [:,np.newaxis]
    upper_perc = da.DimArray(np.percentile(contrib, 95, axis=1),
                             axes=contrib.time, dims="time")  # [obs_period.searchsorted(plot_period)]
    lower_perc = da.DimArray(np.percentile(contrib, 5, axis=1),
                             axes=contrib.time, dims="time")  # [obs_period.searchsorted(plot_period)]
    median = da.DimArray(np.percentile(contrib, 50, axis=1),
                         axes=contrib.time, dims="time")  # [obs_period.searchsorted(plot_period)]

    h = ax6.fill_between(plot_period, lower_perc[plot_period] * 1e3, upper_perc[
                         plot_period] * 1e3, color=rcpcoldict[scen], alpha=.4, lw=0.5)
    ax6.plot(plot_period, median[plot_period] * 1e3, lw=3, color=rcpcoldict[scen],
             alpha=1., label=rcpnamedict[scen])
    # ax6.plot(plot_period,sem_median*1e3,lw=2,color=rcpcoldict[scen],alpha=1.,label=scen+" VR09",marker="|",markeredgecolor="black",markevery=20)

    total_contribs[scen] = np.array(
        [median[2100], lower_perc[2100], upper_perc[2100]])
    # contributions_file.write(scen+" "+rd(median[-1])+" "+rd(lower_perc[-1])+" "+rd(upper_perc[-1])+"\n")

    oloc -= 1.  # plot_offset[name]
    xloc -= 1.  # plot_offset[name]
    low = lower_perc[2081:2100].mean() * 1.e3
    med = median[2081:2100].mean() * 1.e3
    upp = upper_perc[2081:2100].mean() * 1.e3
    axy.plot([oloc - 1, oloc + 1], [med, med], lw=3, color=rcpcoldict[scen])
    # axy.plot([oloc,oloc],[low,upp],lw=3,color=rcpcoldict[scen])
    axy.fill_between([oloc - 1, oloc + 1], [low, low],
                     [upp, upp], color=rcpcoldict[scen], alpha=.4, lw=0.)
    # axy.plot([oloc,oloc],[low,upp],lw=3,color=rcpcoldict[scen])
    low, med, upp = ipcc.get_ipcc_range(scen, "mean_slr_2081_2100")
    # low,med,upp = gd.ipcc_contrib_estimates[name][scen] #gd.get_ipcc_range(scen,name)
    # axy.plot([xloc,xloc],[low,upp],color=rcpcoldict[scen],lw=3,alpha=0.5)
    axy.fill_between([xloc - 1, xloc + 1], [low, low],
                     [upp, upp], color=rcpcoldict[scen], alpha=.4, lw=0.)
    axy.plot([xloc - 1, xloc + 1], [med, med],
             color=rcpcoldict[scen], lw=3, alpha=1.)
    axy.set_xlim(-1, 11)

    # low,med,upp = gd.get_ipcc_range(scen,"mean_slr_2100")
    # xloc -= 1.#plot_offset[name]
    # axy.plot([xloc,xloc],[low,upp],color=rcpcoldict[scen],lw=3)
    # axy.plot([xloc-1,xloc+1],[med,med],color=rcpcoldict[scen],lw=3)
    # axy.set_xlim(-1,6)

axy.text(1.5, 1200, "M15", rotation="vertical", horizontalalignment='center',
         verticalalignment='center')
axy.text(7.5, 1200, "IPCC", rotation="vertical", horizontalalignment='center',
         verticalalignment='center')

# past_sl_anom  = gd.church_past_sl-gd.anthrop_gic_contrib
# ax6.plot(past_sl_anom.time,past_sl_anom*1e3,lw=2,color="black",alpha=1.,label="observed",marker="|",markevery=10)
# # ax6b.plot(obs_period,gd.giss_temp[obs_period],lw=2,color="red",alpha=1.)
# ax6.text(0.05,0.8,"total", transform=ax6.transAxes,
#     fontdict={'family':'sans-serif','weight':'bold'})
ax6.set_xlim(plot_period[0], plot_period[-1])
ax6.set_xlabel("Time in years")
ax6.set_ylabel("Sea level in mm")


l1 = ax6.legend(ncol=1, loc="center left")
l1.draw_frame(0)
for l in l1.get_lines():
    l.set_alpha(1)

plt.draw()
plt.show()

runname = __file__.split("/")[-1][0:-3]
figname = "../figures/" + runname + ".pdf"
plt.savefig(figname, dpi=150)
