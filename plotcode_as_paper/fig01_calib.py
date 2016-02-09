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


""" Code for Fig. 1 as in
    M. Mengel et al.
    Future sea-level rise constrained by observations and long-term commitment
    PNAS (2016)
    (C) Matthias Mengel working at Potsdam Institute for Climate Impact Research

    Note: You may not be able to fully plot this figure as it contains data that is not
    openly available. Request data from authors or comment out.
"""

import os, glob, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import itertools
import cPickle as pickle
from scipy.io import loadmat
lib_path = os.path.abspath('../src')
sys.path.append(lib_path)
import get_calibration_data as gd; reload(gd)
import calib_settings as cs; reload(cs)
import contributor_functions as cf; reload(cf)
import get_gmt_data as ggd; reload(ggd)
import sealevel as sl; reload(sl)
import dimarray as da
import collections
import matplotlib.font_manager as font_manager

path = '/Users/mengel/Library/Fonts/OpenSans-Bold.ttf'
prop = font_manager.FontProperties(fname=path)
plt.rcParams['font.family'] = prop.get_name()
plt.rcParams['font.sans-serif'] = prop.get_name()

plt.rcParams['xtick.major.pad']  = 10
plt.rcParams['font.size']= 12
plt.rcParams['lines.markeredgewidth']=1
plt.rcParams['legend.fontsize']=12
plt.rcParams['figure.figsize'] = 8,10
plt.rcParams['figure.facecolor'] = "white"
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['pdf.fonttype'] = '42'
# plt.rcParams['font.family'] = 'Open Sans'


contrib_ids = ["thermexp","gic","gis_sid", "gis_smb", "ant_sid", "ant_smb"]

labels = {"gic":"Mountain glaciers","thermexp":"Thermal expansion","gis_sid":"Greenland SID",
"gis_smb":"Greenland SMB","ant_sid":"Antarctica SID","ant_smb":"Antarctica SMB"}
legend_loc = {"gic":"lower right","thermexp":0,"gis_sid":"lower right",
"gis_smb":(0,0.45),"ant_sid":"center left","ant_smb":"lower right"}


obslbl = {"leclerqu11":"Leclercq 11","cogley09":"Cogley 09","marzeion12":"Marzeion 12",
    "sasgen12":"Sasgen 12","church11":"Church 11","angelen14":"Angelen 14",
    "box_colgan13":"Box 13","mouginot_rignot14":"Mouginot 14","harig_simons15":"Harig 15"}

plot_anomaly_time = np.arange(1986,2006,1)

fig = plt.figure(1)
plt.clf()
plt.subplots_adjust(left=None, bottom=0.06, right=None, top=0.97,
                  wspace=None, hspace=None)

plotlb = "abcdef"
axs = []
for p,contrib in enumerate(contrib_ids):
    # def get_obs_and_calibrated(contrib):
    print p,contrib
    ax = plt.subplot(3,2,p+1)
    axs.append(ax)


    calibdata = pickle.load(open("../data/calibration/"+contrib+".pkl","rb"))
    observations = cs.__dict__[contrib+"_observations"]
    cols = [cm.Accent_r(np.float(k)/len(observations)) for k in np.arange(len(observations))]


    for i,obs in enumerate(observations):

        if contrib == "thermexp" and i % 2 == 0:
            continue
        print obs
        params = calibdata["params"][obs]

        offset = i*1./400.
        if contrib == "thermexp":
            offset = i*1./75.
        if contrib == "gis_sid":
            offset = i*1./400. + 1./400
        if contrib == "gic":
            offset = i*1./100.
        if contrib == "ant_smb":
            offset = -0.005

        offset *= 1.e3

        if contrib=="thermexp":
            ax.plot(gd.ipcc_thermal_expansion.time, gd.ipcc_thermal_expansion*1.e3 + offset,
                lw=1,color="grey")


        gmt_anom = ggd.giss_temp
        temp_anomaly_year = params.temp_anomaly_year
        # print temp_anomaly_year
        # apply an offset to temperature levels for box & colgan
        if obs == "box_colgan13":
            gmt_anom = ggd.giss_temp + sl.gis_colgan_temperature_offset

        sl_calculated = np.zeros([gmt_anom.shape[0],len(params.commitment_parameter)])
        for j,alpha in enumerate(params.commitment_parameter):
            tau = params.calibrated_parameter[j]

            sl_contributor = params.sl_contributor(alpha,temp_anomaly_year)
            sl_calculated[:,j] = sl_contributor.calc_contribution(gmt_anom,tau)
            sl_calculated[:,j] -= sl_calculated[
                np.searchsorted(gmt_anom.time,plot_anomaly_time),j].mean()
            ax.plot(gmt_anom.time,sl_calculated[:,j]*1.e3+offset,color=cols[i], lw=2, alpha=0.2)

        # plot_anomaly_time = np.maximum(params.cal_period[0],obsv.time[0])
        # calculated_first_year_of_obs = sl_calculated[
        #     np.searchsorted(gmt_anom.time,plot_anomaly_time),:].mean()
        # # print calculated_first_year_of_obs
        # obs_aligned = (obsv - obsv[plot_anomaly_time] +
        #     calculated_first_year_of_obs)
        if contrib == "ant_smb":
            continue

        obsv = observations[obs]
        if contrib == "gic":
            obsv = (obsv.diff()*gd.anthropogenic_frac).dropna().cumsum()
            #print obsv.values

        if obs == "harig_simons15":
            ## only starts in 2003, so 1986-2005 mean not available,
            ## use that of observation before (dirty!)
            obs_aligned = (obsv - obsv[2003] +
                sl_calculated[np.searchsorted(gmt_anom.time,2003),:].mean())
        else:
            obs_aligned = obsv - obsv[plot_anomaly_time].mean()
        #     calculated_first_year_of_obs)

        lbl = obslbl[obs] if contrib != "thermexp" else ""
        ax.plot(obsv.time, obs_aligned*1.e3 + offset,label=lbl,lw=1,
            color=cols[i],marker="|",markevery=5,markersize=10)

    if contrib not in ["ant_smb","thermexp"]:
        ncol = 2 if contrib == "thermexp" else 1
        l2 = ax.legend(loc=legend_loc[contrib],ncol=ncol)
        l2.draw_frame(0)


    ax.text(0.05,0.9,plotlb[p], transform=ax.transAxes,
        fontdict={'family':'sans-serif','weight':'bold', "size":16})
    ax.text(0.05,0.8,labels[contrib], transform=ax.transAxes,
        fontdict={'family':'sans-serif','weight':'bold'})
    # if contrib == "thermexp":
    #     ax.set_ylim(top=0.45)
    if "ant" not in contrib:
        ax.set_ylim(ax.get_ylim()[0]*0.8,ax.get_ylim()[1]*0.8)
    if "ant_smb" == contrib:
        ax.set_ylim(-15,5)
    if "gic" == contrib:
        ax.set_ylim(-40,25)
    if "thermexp" == contrib:
        ax.set_ylim(-30,180)

    ax.set_xlim(1880,2015)


for ax in axs[4:6]:
    ax.set_xlabel("Time in years")
for ax in [axs[0],axs[2],axs[4]]:
    ax.set_ylabel("Sea level contribution in mm")

plt.draw()
plt.show()

runname = __file__.split("/")[-1][0:-3]
# plt.savefig("../figures/"+runname+"2.svg")
# plt.savefig("../figures/"+runname+"2.png")
plt.savefig("../figures/"+runname+"5.pdf")
