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
          run fig03_rcp_components.py -g -n numberofsamples
"""


import os, glob, sys
import numpy as np
# import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib import cm
import itertools
from scipy.io import loadmat
lib_path = os.path.abspath('../src')
sys.path.append(lib_path)
import contributor_functions as cf; reload(cf)
#import get_calibration_data as gd; reload(gd)
import get_ipcc_data as ipcc; reload(ipcc)
import dimarray as da
import optparse
from mpl_toolkits.axes_grid1 import make_axes_locatable
import timeit
import cPickle as pickle
import pandas as pd
import matplotlib.font_manager as font_manager

# make use of -n or -g flags with this script
parser = optparse.OptionParser()
parser.add_option('-n', '--realizations',  type="int",  help='get number of realizations')
parser.add_option('-g', help='get new data', dest='get_data', default=False, action='store_true')
(opts, args) = parser.parse_args()

path = '/Users/mengel/Library/Fonts/OpenSans-Bold.ttf'
prop = font_manager.FontProperties(fname=path)
plt.rcParams['font.family'] = prop.get_name()
plt.rcParams['font.sans-serif'] = prop.get_name()

# plot settings
plt.rcParams['xtick.major.pad']  = 10
plt.rcParams['font.size']= 12
plt.rcParams['lines.markeredgewidth']=2
plt.rcParams['legend.fontsize']=12
plt.rcParams['figure.figsize'] = 10,10
plt.rcParams['figure.facecolor'] = "white"
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['pdf.fonttype'] = '42'

contrib_ids = ["thermexp","gic","gis_sid","gis_smb", "ant_sid","ant_smb", ]
rcpnamedict = {"RCP3PD":"RCP26","RCP45":"RCP45","RCP85":"RCP85"}

rcpcoldict = {"RCP3PD":"#2256A6","RCP45":"#73B2E1","RCP85":"#EE322D"}
labels = {"gic":"Mountain glaciers","thermexp":"Thermal expansion","gis_sid":"Greenland SID",
"gis_smb":"Greenland SMB","ant_sid":"Antarctica SID","ant_smb":"Antarctica SMB"}

if opts.get_data:
  realizations = 10000
  realizations = realizations  if opts.realizations==None else opts.realizations
  print "use",realizations,"realizations"
  projection_data = pickle.load(open("../data/projection/projected_slr_"+
                                str(realizations)+"samples.pkl","rb"))

# plt.close("all")
plt.figure(3,figsize=(8,10))
plt.clf()
plt.subplots_adjust(left=0.1, bottom=0.06, right=1.0, top=0.97,
                  wspace=0.1, hspace=None)
# save numbers in single_contribs to print later
# single_contribs = {}

plot_period = np.arange(2000,2101,1)
axs = []
plotlb = "abcdef"
for i,name in enumerate(contrib_ids[0:6]):
  ax = plt.subplot(3,2,i+1)
  # ax = plt.subplot(111)
  divider = make_axes_locatable(ax)
  axy = divider.append_axes("right", size=0.6, pad=0.0, sharey=ax,axisbg="grey")
  axy.axis("off")
  axy.axvspan(0, 10, facecolor='0.5', alpha=0.2,lw=0)
  # print name
  xloc=10 # ipcc
  oloc=4 # our estimates
  # single_contribs[name] = {}
  for k,scen in enumerate(["RCP85","RCP45","RCP3PD"]):

    contrib    = projection_data[scen][name]
    # plot_period_indices = contrib.time.searchsorted(plot_period)
    # anomaly to 1986-2005
    contrib    = contrib - contrib[1986:2005,:].mean(axis="time")
    # pandas handles NaNs within percentile calculus
    contrib     = contrib.T.to_pandas()
    percentiles = da.from_pandas(contrib.quantile([0.05,0.5,0.95]))[:,plot_period]
    ## plot in mm SLR
    ax.fill_between(plot_period,percentiles[0.05]*1e3,percentiles[0.95]*1e3,color=rcpcoldict[scen],alpha=.4,lw=.5)
    ax.plot(plot_period,percentiles[0.5]*1e3,lw=3,color=rcpcoldict[scen],alpha=1.,label=rcpnamedict[scen])#,label=name)

    print "##",name,scen,": ", percentiles[0.05][2100]*1e3,percentiles[0.5][2100]*1e3,percentiles[0.95][2100]*1e3
    # single_contribs[name][scen] = np.array([median[-1],lower_perc[-1],upper_perc[-1]])

    # if name != "gis":
    oloc -= 1.#plot_offset[name]
    xloc -= 1.#plot_offset[name]
    med = percentiles[0.5][np.arange(2081,2100,1)].mean()*1e3
    low = percentiles[0.05][np.arange(2081,2100,1)].mean()*1e3
    upp = percentiles[0.95][np.arange(2081,2100,1)].mean()*1e3
    axy.plot([oloc-1,oloc+1],[med,med],lw=3,color=rcpcoldict[scen])
    # axy.plot([oloc,oloc],[low,upp],lw=3,color=rcpcoldict[scen])
    axy.fill_between([oloc-1,oloc+1],[low,low],[upp,upp],color=rcpcoldict[scen],alpha=.4,lw=0.)
    # axy.plot([oloc,oloc],[low,upp],lw=3,color=rcpcoldict[scen])
    low,med,upp = ipcc.ipcc_contrib_estimates[name][scen] #gd.get_ipcc_range(scen,name)
    # axy.plot([xloc,xloc],[low,upp],color=rcpcoldict[scen],lw=3,alpha=0.5)
    axy.fill_between([xloc-1,xloc+1],[low,low],[upp,upp],color=rcpcoldict[scen],alpha=.4,lw=0.)
    axy.plot([xloc-1,xloc+1],[med,med],color=rcpcoldict[scen],lw=3,alpha=1.)
    axy.set_xlim(-1,11)

  # contributions[name] = contrib_param
  # ax.text(0.05,0.8,name, transform=ax.transAxes,fontweight='bold',)
  ax.text(0.05,0.9,plotlb[i], transform=ax.transAxes,
      fontdict={'family':'sans-serif','weight':'bold', "size":16})
  ax.text(0.05,0.8,labels[name], transform=ax.transAxes,
      fontdict={'family':'sans-serif','weight':'bold'})
  if i==0:
    axy.text(1.5,440,"M16",rotation="vertical",horizontalalignment='center',
      verticalalignment='center')
    axy.text(7.5,440,"IPCC",rotation="vertical",horizontalalignment='center',
      verticalalignment='center')

  l1 = ax.legend(ncol=1,loc="center left")
  l1.draw_frame(0)
  for l in l1.get_lines(): l.set_alpha(1)
  # ax.set_ylabel("sea level in mm")
  # ax.set_xlabel("time in yr")
  # if name != "Antarctica SMB":
  #   ax.set_ylim(bottom=0.0)
  # ax.set_xticklabels([])
  axs.append(ax)
  ax.set_xlim(plot_period[0],plot_period[-1])

for ax in axs[-2:]:
  ax.set_xlabel("Time in years")
for ax in [axs[0],axs[2],axs[4]]:
  ax.set_ylabel("Sea level in mm")

axs[3].set_ylim(0,500)
# for ax in axs[0:2]:
#   ax.set_ylim(0,500)
#   for ax in axs[3:6]:
#     ax.set_ylim(0,200)

plt.draw()
plt.show()

runname = __file__.split("/")[-1][0:-3]
figname = "../figures/"+runname+".pdf"
plt.savefig(figname,dpi=150)
