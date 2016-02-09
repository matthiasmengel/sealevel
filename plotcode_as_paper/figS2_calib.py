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

""" Code for Fig. S2 as in
    M. Mengel et al.
    Future sea-level rise constrained by observations and long-term commitment
    PNAS (2016)
    (C) Matthias Mengel working at Potsdam Institute for Climate Impact Research

    Note: You may not be able to fully plot this figure as it contains data that is not
          openly available. Request data from authors or comment out.
"""

import os, glob, sys
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib import cm
import itertools
import cPickle as pickle
from scipy.io import loadmat
lib_path = os.path.abspath('../src')
sys.path.append(lib_path)
#import get_calibration_data as gd; reload(gd)
import get_gmt_data as ggd; reload(ggd)
import sealevel as sl; reload(sl)
import calib_settings as cs; reload(cs)
import contributor_functions as cf; reload(cf)
import dimarray as da
import collections

plt.rcParams['xtick.major.pad']  = 10
plt.rcParams['font.size']= 12
plt.rcParams['lines.markeredgewidth']=2
plt.rcParams['legend.fontsize']=12
plt.rcParams['figure.figsize'] = 10,10
plt.rcParams['figure.facecolor'] = "white"
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['pdf.fonttype'] = '42'

contrib_ids = ["thermexp","gic","gis_smb", "gis_sid", "ant_smb", "ant_sid"]

## Note: "gic" has specific plot script.
plot_these = ["thermexp"]#,"thermexp","gic","gis_smb","gis_sid"]



for contrib in plot_these:
    # def get_obs_and_calibrated(contrib):

    calibdata = pickle.load(open("../data/calibration/"+contrib+".pkl","rb"))
    observations = cs.__dict__[contrib+"_observations"]
    cols = [cm.Accent(np.float(k)/len(observations)) for k in np.arange(len(observations))]

    fig = plt.figure(5)
    plt.subplots_adjust(bottom=0.06,top=0.95)
    plt.clf()
    ax = plt.subplot(111)

    for i,obs in enumerate(observations):

        print obs
        params = calibdata["params"][obs]

        offset = i*1./400.
        if contrib == "thermexp":
            offset = i*1./75.
        if contrib == "gis_sid":
            offset = i*1./400. + 1./400

        gmt_anom = ggd.giss_temp
        temp_anomaly_year = params.temp_anomaly_year
        print temp_anomaly_year
        # apply an offset to temperature levels for box & colgan
        if obs == "box_colgan13":
            gmt_anom = ggd.giss_temp + sl.gis_colgan_temperature_offset

        sl_calculated = np.zeros([gmt_anom.shape[0],len(params.commitment_parameter)])
        for j,alpha in enumerate(params.commitment_parameter):
            tau = params.calibrated_parameter[j]

            sl_contributor = params.sl_contributor(alpha,temp_anomaly_year)
            sl_calculated[:,j] = sl_contributor.calc_contribution(gmt_anom,tau)
            ax.plot(gmt_anom.time,sl_calculated[:,j]+offset,color=cols[i])


        obsv = observations[obs]
        plot_anomaly_time = np.maximum(params.cal_period[0],obsv.time[0])
        calculated_first_year_of_obs = sl_calculated[
            np.searchsorted(gmt_anom.time,plot_anomaly_time),:].mean()

        # print calculated_first_year_of_obs
        obs_aligned = (obsv - obsv[plot_anomaly_time] +
            calculated_first_year_of_obs)

        ax.plot(obsv.time, obs_aligned +offset,label=obs,lw=2,
            color=cols[i],marker="|",markevery=10,markersize=10)


    ncol = 2 if contrib == "thermexp" else 1
    l2 = ax.legend(loc="upper left",ncol=ncol)
    l2.draw_frame(0)

    if contrib == "thermexp":
        ax.set_ylim(-0.01,0.22)

    ax.set_xlabel("Time in years")
    ax.set_ylabel("Sea-level contribution in m")

    ax.set_xlim(1880,2015)

    plt.draw()
    plt.show()

    runname = __file__.split("/")[-1][0:-3]
    plt.savefig("../figures/"+runname+"_"+contrib+"2.pdf")
