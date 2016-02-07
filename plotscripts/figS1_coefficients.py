
""" matthias.mengel@pik
    currently excludes the Surface Mass Balance of Antarctica
"""

import os, glob, sys
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import ticker
import itertools
from scipy.io import loadmat
lib_path = os.path.abspath('../code_python')
sys.path.append(lib_path)
import single_contributions as sc; reload(sc)
import get_data as gd; reload(gd)
import dimarray as da
import optparse
from mpl_toolkits.axes_grid1 import make_axes_locatable

parser = optparse.OptionParser()
parser.add_option('-g', help='get new data', dest='get_data', default=False, action='store_true')
(opts, args) = parser.parse_args()


runname = __file__.split("/")[-1][0:-3]
figname = "../figures_paper/"+runname

plt.rcParams['xtick.major.pad']  = 10
plt.rcParams['font.size']= 10
plt.rcParams['lines.markeredgewidth']=2
plt.rcParams['legend.fontsize']=10
plt.rcParams['figure.figsize'] = 12,10
plt.rcParams['figure.facecolor'] = "white"


markers = ['+', '*', 'o', 'x', '.', '|',"v","h","^","p","_",'D']
processdict = {}
for m,lbl in enumerate(["ThermalExpansion","Mountainglaciers","Softening","Smb_Dynamic","Smb_Melting","Smb_Accumulation","Oceanwarming","Sum"]):
  processdict[lbl] = markers[m]
## colors for components
componentdict  = {"Greenland":"#84ce28","Antarctica":"#556eb9","Mountain glaciers (anth)":"#6dcfd7","Thermal expansion":"#ee252b","Greenland SID":"#84ce28","Greenland SMB":"#84ce28","Antarctica SID":"#556eb9",}

classname = {"Mountain glaciers (anth)":"glaciers_and_icecaps","Thermal expansion":"thermal_expansion","Greenland SID":"solid_ice_discharge_gis","Greenland SMB":"surfacemassbalance_gis","Antarctica SID":"solid_ice_discharge_ais","Antarctica SMB":"surfacemassbalance_ais"}

contrib_names = ["Thermal expansion","Mountain glaciers (anth)","Greenland SMB", "Greenland SID", "Antarctica SID"]

indep_names =  {
"Thermal expansion":"alpha te",
"Mountain glaciers (anth)":"model number",
"Greenland SMB": "smb coefficient",
 "Greenland SID": "alpha sid",
 "Antarctica SID": "alpha sid"}

dep_names   = {
"Thermal expansion":"tau te",
"Mountain glaciers (anth)":"tau gic",
"Greenland SMB":"tau smb",
 "Greenland SID":"prefactor",
 "Antarctica SID":"tau sid"}




# rcpcoldict = {"RCP3PD":"#3b4ba8","RCP45":"#8cbbe7","RCP60":"#f46f25","RCP85":"#ee1d23"}
rcpcoldict = {"RCP3PD":"#3b4ba8","RCP45":"yellow","RCP85":"#ee1d23"}

# plot_offset = {"ThermalExpansion_Schewe":1,"Mountainglaciers_Marzeion":2,"Antarctica_Sum":3,"Greenland_Sum":4}

sm_iceage_offset="0p0"


# plt.close("all")
plt.figure(5,figsize=(12,12))
plt.clf()
plt.subplots_adjust(wspace=0.4)
y_formatter = ticker.ScalarFormatter(useOffset=True)
axs = []
for i,name in enumerate(contrib_names):
  ax = plt.subplot(3,2,i+1)
  ax.yaxis.set_major_formatter(y_formatter)
  print name
  calib_data  = np.load("../data_python/"+classname[name]+sm_iceage_offset+".npz")
  ind_params  = calib_data["independent_param"]
  dep_params  = calib_data["dependent_param"]

  # let model number start at one
  if name == "Mountain glaciers (anth)":
    ind_params += 1.

  # col = componentdict[name]

  # for k,scen in enumerate(["RCP85","RCP45","RCP3PD"]):

  #   upper_perc = np.percentile(contributions[scen][:,i,:],95,axis=0)
  #   lower_perc = np.percentile(contributions[scen][:,i,:],05,axis=0)
  #   median     = np.percentile(contributions[scen][:,i,:],50,axis=0)
  #   ax.fill_between(obs_period,lower_perc*1e3,upper_perc*1e3,color=rcpcoldict[scen],alpha=.5)
  ax.plot(ind_params,dep_params,".",label=name,markersize=10)#,label=name)
    # for n in np.arange(realizations):

    #   contrib = contributions[scen][n,i,:]
    #   # lbl = "fitted" if k==0 else ""
    #   ax.plot(obs_period,contrib*1e3,lw=3,color=rcpcoldict[scen],alpha=.5,label=lbl)#,label=name)

  # ytikcs
  if name == "Greenland SID":
    locs,labels = plt.yticks()
    plt.yticks(locs, map(lambda x: "%.1e" % x, locs))
    ax.set_xlim([ax.get_xlim()[0]*1.05,ax.get_xlim()[1]*0.95])
  else:
    # ylabel('microseconds (1E-9)')
    ax.set_xlim([ax.get_xlim()[0]*.95,ax.get_xlim()[1]*1.05])
  # contributions[name] = contrib_param

  if name == "Mountain glaciers (anth)":
    ax.set_xlim(left=.8)
    ax.set_xticks([1,2,3,4])

  ax.set_ylabel(dep_names[name])
  ax.set_xlabel(indep_names[name])
  ax.text(0.05,0.8,name, transform=ax.transAxes,fontweight='bold',)
  # l1 = plt.legend(ncol=1,loc="center left")
  # l1.draw_frame(0)
  # for l in l1.get_lines(): l.set_alpha(1)
  # ax.set_ylabel("sea level in mm")
  # if name != "Antarctica SMB":
  #   ax.set_ylim(bottom=0.0)
  # ax.set_xticklabels([])
  # ax.get_yaxis().get_major_formatter().set_scientific(True)
  axs.append(ax)


figname = "../figures_paper/"+runname+"_"+sm_iceage_offset+".png"
plt.savefig(figname,dpi=150)

