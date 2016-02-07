
""" matthias.mengel@pik
    currently excludes the Surface Mass Balance of Antarctica
"""

import os, glob, sys
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib import cm
import itertools
from scipy.io import loadmat
lib_path = os.path.abspath('../code_python')
sys.path.append(lib_path)
import single_contributions as sc; reload(sc)
import get_data as gd; reload(gd)
import dimarray as da

plt.rcParams['xtick.major.pad']  = 10
plt.rcParams['font.size']= 10
plt.rcParams['lines.markeredgewidth']=2
plt.rcParams['legend.fontsize']=10
plt.rcParams['figure.figsize'] = 8,12
plt.rcParams['figure.facecolor'] = "white"


markers = ['+', '*', 'o', 'x', '.', '|',"v","h","^","p","_",'D']
processdict = {}
for m,lbl in enumerate(["ThermalExpansion","Mountainglaciers","Softening","Smb_Dynamic","Smb_Melting","Smb_Accumulation","Oceanwarming","Sum"]):
  processdict[lbl] = markers[m]
## colors for components
componentdict  = {"Greenland":"#84ce28","Antarctica":"#556eb9","Mountain glaciers (anth)":"#6dcfd7","Thermal expansion":"#ee252b","Greenland SID":"#84ce28","Greenland SMB":"#84ce28","Antarctica SID":"#556eb9","Antarctica SMB":"#556eb9",}

classname = {"Mountain glaciers (anth)":"glaciers_and_icecaps","Thermal expansion":"thermal_expansion","Greenland SID":"solid_ice_discharge_gis","Greenland SMB":"surfacemassbalance_gis","Antarctica SID":"solid_ice_discharge_ais","Antarctica SMB":"surfacemassbalance_ais"}

contrib_names = ["Thermal expansion","Mountain glaciers (anth)","Greenland SMB", "Greenland SID", "Antarctica SID", "Antarctica SMB"]


obs_period          = np.arange(1900,2009)
sm_iceage_offset = "0p0"

# plt.close("all")
fig = plt.figure(2,figsize=(8,12))
plt.clf()
plt.subplots_adjust(left=None, bottom=0.05, right=None, top=0.95,
                  wspace=0.25, hspace=None)

contributions = {}
axs = []
for i,name in enumerate(contrib_names):
  ax = plt.subplot(6,1,i+1)
  print name

  # obs_component = gd.church_observed[name]
  col = componentdict[name]

  # obs_period     = gd.observation_period[name]
  # giss_temp_anom = gd.giss_temp[obs_period] - gd.giss_temp[1880:1900].mean()

  # if name == "Thermal expansion":
  rclass = sc.__dict__[classname[name]]

  ## load data
  calib_data = np.load("../data_python/"+classname[name]+sm_iceage_offset+".npz")
  ind_param  = calib_data["independent_param"]
  dep_param  = calib_data["dependent_param"]

  contrib_param = np.zeros([len(ind_param),len(obs_period)])
  for j,indp in enumerate(ind_param):
    if "Greenland SID" == name:
      tau     = np.random.randint(41)+10 # free sampled time lag between 10 and 50
      cl      = rclass(indp,tau)
    else:
      cl      = rclass(indp)

    contrib            = cl.calc_contribution(gd.giss_temp[obs_period],dep_param[j])
    contrib_param[j,:] = contrib
    contrib            = da.DimArray(contrib,dims="time",axes=obs_period)
    contrib           -= contrib[1986:2005].mean()
    # if "Greenland SID" == name:
    #   print contrib.values
    # add offset from observations
    # contrib = contrib + obs_component[obs_period[0]]
    lbl = "fitted" if j==0 else ""

    ax.plot(obs_period,contrib*1e3,lw=2,color=col,alpha=0.5,label=lbl)#,label=name)

  contributions[name] = contrib_param

  if name == "Mountain glaciers (anth)":
    marzeion_obs = gd.anth_gic_contrib[obs_period] - gd.anth_gic_contrib[1986:2005].mean()
    ax.plot(obs_period,marzeion_obs*1e3,lw=1,color="grey",alpha=.5,label=lbl)
  if name == "Thermal expansion":
    te_ipcc = gd.ipcc_thermal_expansion - gd.ipcc_thermal_expansion[1986:2005].mean()
    ax.plot(te_ipcc.time,te_ipcc*1e3,lw=1,color="grey",alpha=.5,label=lbl)
    te_lower = gd.ipcc_te_lower - gd.ipcc_te_lower.ix[-1] + contrib[1990]*1e3
    te_upper = gd.ipcc_te_upper - gd.ipcc_te_upper.ix[-1] + contrib[1990]*1e3
    ax.fill_between(te_lower.time,te_lower,te_upper,color="grey",alpha=.2,lw=0.5)

  ax.text(0.05,0.8,name, transform=ax.transAxes,fontweight='bold',)
  # l1 = plt.legend(ncol=1,loc="center left")
  # l1.draw_frame(0)
  # for l in l1.get_lines(): l.set_alpha(1)
  ax.set_ylabel("sea level in mm")
  ax.set_xlim(obs_period[0],obs_period[-1])
  if i != 5:
    ax.set_xticklabels([])
  else:
    ax.set_xlabel("time in years")
 # ax.set_xticklabels([])
  axs.append(ax)

# ax6 = plt.subplot(7,1,7)
# # ax6b = ax6.twinx()

# # total_sl = contributions.sum(axis=0)
# realizations = 100
# for n in np.arange(realizations):

#   total_sl = np.zeros([6,len(obs_period)])
#   # for j,indp in enumerate(ind_param):
#   for i,name in enumerate(contrib_names):
#     contrib_param = contributions[name]
#     chosen        = np.random.randint(contrib_param.shape[0])
#     total_sl[i,:] = contrib_param[chosen,:]

#     # if "SID" in name:
#     #   # can occur or not
#     #   total_sl[i,:] = total_sl[i,:] *np.random.randint(2)
#   total_sl_an = da.DimArray(total_sl.sum(axis=0),dims="time",axes=obs_period)
#   total_sl_an = total_sl_an - total_sl_an[1986:2005].mean()
#   lbl = "fitted" if n == 0 else ""
#   ax6.plot(obs_period,total_sl_an*1e3,lw=2,color="#c98345",alpha=.1,label=lbl)


# past_sl_anom  = gd.church_past_sl[obs_period]-gd.nat_gic_contrib[obs_period]#gd.church_past_sl[1950]#gd.church_past_sl[1880:1890].mean()
# past_sl_anom -= past_sl_anom[1986:2005].mean()
# lower_bound   = past_sl_anom-gd.church_past_sl_std[obs_period]
# upper_bound   = past_sl_anom+gd.church_past_sl_std[obs_period]

# # past_sl_anom -= past_sl_anom[1990:2000].mean() #- total_sl.sum(axis=0)[-10:].mean()
# ax6.plot(obs_period,past_sl_anom*1e3,lw=2,color="black",alpha=1.,label="observed\nanthropogenic",marker="|",markevery=10)
# # ax6.plot(obs_period,lower_bound*1e3,lw=.5,color="black",alpha=1.)
# # ax6.plot(obs_period,upper_bound*1e3,lw=.5,color="black",alpha=1.)
# ax6.fill_between(obs_period,lower_bound*1e3,upper_bound*1e3,lw=.5,color="black",alpha=.3)
# # ax6.plot(obs_period,gd.anthrop_gic_contrib[obs_period]*1e3,lw=2,color="green",alpha=1.,label="gic non anthropogenic")
# # ax6b.plot(obs_period,gd.giss_temp[obs_period],lw=2,color="red",alpha=1.)
# ax6.text(0.05,0.8,"total", transform=ax6.transAxes,fontweight='bold',)
# ax6.set_xlim(obs_period[0],obs_period[-1]+1)
# ax6.set_xlabel("time in years")
# ax6.set_ylabel("sea level in mm")
# # ax6b.set_ylabel("gmt")
# l1 = ax6.legend(ncol=2,loc="upper center",bbox_to_anchor=(0.35,1.0))
# l1.draw_frame(0)
# for l in l1.get_lines(): l.set_alpha(1)

# # for ax in axs:
# #   ax.set_ylim(ax6.get_ylim())

runname = __file__.split("/")[-1][0:-3]
figname = "../figures_paper/"+runname+"_"+sm_iceage_offset+".png"

plt.savefig(figname)
