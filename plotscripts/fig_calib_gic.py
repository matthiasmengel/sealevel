
""" matthias.mengel@pik
"""

import os, glob, sys
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib import cm
import itertools
import cPickle as pickle
# from scipy.io import loadmats
lib_path = os.path.abspath('../src')
sys.path.append(lib_path)
import get_data as gd; reload(gd)
import calib_settings as cs; reload(cs)
import contributor_functions as cf; reload(cf)
import dimarray as da
import collections


plt.rcParams['xtick.major.pad']  = 10
plt.rcParams['font.size']= 12
plt.rcParams['lines.markeredgewidth']=2
plt.rcParams['legend.fontsize']=12
plt.rcParams['figure.figsize'] = 8,10
plt.rcParams['figure.facecolor'] = "white"
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['pdf.fonttype'] = '42'
runname = __file__.split("/")[-1][0:-3]

calibdata = pickle.load(open("../calibrationdata/gic.pkl","rb"))
# temp_anomaly_year = gd.temp_anomaly_year["gic"]
# all timeseries plotted as anomalies to offsettime
offsettime=2000

# plt.close("all")
plt.figure(7)
plt.clf()
ax = plt.subplot(211)

cols = ["red","green","blue","black"]
namelbl = {"leclerqu11":"Leclercq et al. 2011","marzeion12":"Marzeion et al. 2012",
"cogley09":"Cogley 2009"}

for i,obs in enumerate(cs.gic_observations):

    temp_anomaly_year = cs.temp_anomaly_year["gic"][obs]
    obs_period = cs.observation_period["gic"][obs]

    # gic_par = gicdata["gic_params"][obs]
    # for modelno in gic_par[0]:
    #     tau = gic_par[1][modelno]
    #     gic_cont = cf.glaciers_and_icecaps(modelno)
    #     gic_calculated = da.DimArray(gic_cont.calc_contribution(gd.giss_temp,tau),axes=gd.giss_temp.time,dims="time")
    #     gic_calculated -= gic_calculated[2000]
    #     ax.plot(gic_calculated.time,gic_calculated+i/100.,alpha=0.5,color=cols[i])
    dat = cs.gic_observations[obs]#[obs_period]
    dat -= dat[offsettime]
    ax.plot(dat.time,dat+i/100.,label=namelbl[obs],color=cols[i],lw=2)

ax.set_ylim(-0.11,0.045)
# plt.grid()
l1 = plt.legend(loc=0,ncol=2)
l1.draw_frame(0)

# gmt_anom = gd.giss_temp
# gmt_anom = da.DimArray(np.linspace(0,0,gd.giss_temp.shape[0]),axes = gd.giss_temp.time,dims="time")-0.8
gmt = gd.giss_temp #- gd.giss_temp[1880:1900].mean()
ax2 = plt.subplot(212)
for i,obs in enumerate(cs.gic_observations):


    temp_anomaly_year = cs.temp_anomaly_year["gic"][obs]
    obs_period = cs.observation_period["gic"][obs]
    # print obs
    # if obs != "marzeion12":
    #     continue


    # dat = gd.gic_observations[obs]#*gd.anthropogenic_fraction
    dat = (cs.gic_observations[obs].diff()*gd.anthropogenic_frac).dropna().cumsum()
    # dat2 =  cs.gic_observations[obs]

    # offsettime = dat.time[0]
    # if offsettime < 1940:
    #     offsettime = 1940


    # dat -= dat.ix[0]
    dat = dat - dat[offsettime]
    lbl = "antropogenic fraction\nof observation" if i==2 else ""
    ax2.plot(dat.time,dat + i/100.,color=cols[i],lw=2,ls="--", label=lbl)
    # ax2.plot(dat2.time,dat2 + i/100.,color=cols[i],lw=2,ls="-")

    params = calibdata["params"][obs]
    for modelno in params.commitment_parameter:
        tau = params.calibrated_parameter[modelno]
        # print tau
        gic_cont = params.sl_contributor(modelno,temp_anomaly_year)

        lbl = "reconstructed\nfrom GMT" if modelno==0 and i==2 else ""

        gic_calculated = da.DimArray(gic_cont.calc_contribution(gmt,tau),axes=gd.giss_temp.time,dims="time")
        gic_calculated -= gic_calculated[offsettime] #- dat[offsettime]

        # if modelno ==0 :
        #     print gic_calculated.values

        ax2.plot(gic_calculated.time,gic_calculated+i/100.,alpha=0.5,color=cols[i],label=lbl)

# plt.grid()
# ax2.set_ylim(ax.get_ylim())

l2 = ax2.legend(loc=0,ncol=2)
l2.draw_frame(0)

ax2.set_xlabel("Time in years")
ax.set_ylabel("Sea-level contribution in m")
ax2.set_ylabel("Sea-level contribution in m")
# ax.set_xticklabels([])

ax.set_xlim(1870,2015)
ax2.set_xlim(1870,2015)

ax.text(0.05,0.75,"observations", transform=ax.transAxes,fontweight='bold',fontsize=16)
ax2.text(0.05,0.75,"anthropogenic fraction", transform=ax2.transAxes,fontweight='bold',fontsize=16)


# plt.clf()
# for i,obs_gic in enumerate(gd.gic_observations):
#     print "####",obs_gic
#     sl_observation = gd.gic_observations[obs_gic]
#     sl_observation -= sl_observation[1961]
#     plt.plot(sl_observation.time,sl_observation-sl_observation[1961],color=cols[i])

#     anthro = (sl_observation.diff()*gd.anthropogenic_fraction).dropna().cumsum()
#     # anthro -= anthro.ix[0]
#     plt.plot(anthro.time,anthro,"--",color=cols[i])


# plt.grid()
# plt.legend(loc=0,ncol=3)
plt.savefig("../figures/"+runname+".png")
plt.savefig("../figures/"+runname+".pdf")

plt.draw()
plt.show()

