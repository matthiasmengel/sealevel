
""" matthias.mengel@pik
    currently excludes the Surface Mass Balance of Antarctica
"""

# import os, glob, sys
import numpy as np
# import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib import cm
# import itertools
# from scipy.io import loadmat
import get_data as gd; reload(gd)
import calib_settings as cs; reload(cs)
import dimarray as da
import cPickle as pickle
import matplotlib.font_manager as font_manager

path = '/Users/mengel/Library/Fonts/OpenSans-Bold.ttf'
prop = font_manager.FontProperties(fname=path)
plt.rcParams['font.family'] = prop.get_name()
plt.rcParams['font.sans-serif'] = prop.get_name()

plt.rcParams['xtick.major.pad']  = 10
plt.rcParams['font.size']= 12
plt.rcParams['lines.markeredgewidth']=2
plt.rcParams['legend.fontsize']=12
plt.rcParams['figure.figsize'] = 8,8
plt.rcParams['figure.facecolor'] = "white"
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['pdf.fonttype'] = '42'


contrib_ids = ["thermexp","gic","gis_smb", "gis_sid", "ant_smb", "ant_sid"]

def project(dummynumber):

    """ Same as in create_projections2.py, only use giss_gmt as driving temperature.
      for one random choice of observations obs_choice.
      and one random tuple of independent and dependent parameter.
    """

    driving_temperature = gd.giss_temp

    # print contrib_name, temp_anomaly_year
    # use one of the observational dataset
    obs_choice = np.random.choice(calibdata["params"].keys())
    params = calibdata["params"][obs_choice]
    temp_anomaly_year = params.temp_anomaly_year

    if obs_choice == "box_colgan13":
        driving_temperature += gd.gis_colgan_temperature_offset

    # print np.max(driving_temperature)

    # choose a random parameter set
    tuple_choice = np.random.randint(len(params.commitment_parameter))
    independent_param = params.commitment_parameter[tuple_choice]
    dependent_param   = params.calibrated_parameter[tuple_choice]
    # print "tuple",independent_param,dependent_param

    contributor = params.sl_contributor(independent_param, temp_anomaly_year)
    contrib = contributor.calc_contribution(driving_temperature,dependent_param)
    # print contrib
    return contrib

### Monte Carlo Sampling from all the different observational datasets and
### fitting/calibrated parametes.
all_contributions = {}
realizations = 10000
for i,contrib_name in enumerate(contrib_ids):

    print "conribution", contrib_name
    calibdata = pickle.load(open("../calibrationdata/"+contrib_name+".pkl","rb"))

    # dd = p.map_async(project,selected_numbers)
    proj = map(project,np.arange(realizations))
    all_contributions[contrib_name] = da.DimArray(proj,
        axes=[np.arange(realizations),gd.giss_temp.time],dims=["realization","time"])

all_contributions = da.DimArray(all_contributions,dims=["contribution","realization","time"])


obs_period          = np.arange(1900,2009)
# plt.close("all")
fig = plt.figure(2)
plt.clf()
ax = plt.subplot(111)

natural_gic_contrib = gd.marzeion_gic_nat_up.mean(axis="model")/1.e3
past_sl_anom  = gd.church_past_sl[obs_period]-natural_gic_contrib[obs_period]#gd.
past_sl_anom -= past_sl_anom[1986:2005].mean()
lower_bound   = past_sl_anom-gd.church_past_sl_std[obs_period]
upper_bound   = past_sl_anom+gd.church_past_sl_std[obs_period]

#### observed ####
# past_sl_anom -= past_sl_anom[1990:2000].mean() #- total_sl.sum(axis=0)[-10:].mean()
ax.plot(obs_period,past_sl_anom*1e3,lw=2,color="black",alpha=1.,
    label="observed\nChurch et al.",marker="|",markevery=10,markersize=10)

ax.fill_between(obs_period,lower_bound*1e3,upper_bound*1e3,lw=.5,color="black",alpha=.3)

hay15_past_sl = gd.hay15_past_sl - natural_gic_contrib[obs_period]
hay15_past_sl -= hay15_past_sl[1986:2005].mean()
ax.plot(hay15_past_sl.time,hay15_past_sl*1e3,lw=2,color="blue",alpha=1.,
    label="observed\nHay et al.",marker="x",markevery=10,markeredgewidth=2,markersize=5)

hay15_lower = gd.hay15_lower- natural_gic_contrib[obs_period] - gd.hay15_past_sl[1986:2005].mean()
hay15_upper = gd.hay15_upper- natural_gic_contrib[obs_period] - gd.hay15_past_sl[1986:2005].mean()
# hay15_past_sl -= hay15_lower[1986:2005].mean()
ax.fill_between(hay15_lower.time,hay15_lower*1e3,hay15_upper*1e3,lw=.5,color="blue",alpha=.2)

# ax.plot(hay15_past_sl.time,hay15_lower*1e3,lw=2,color="blue",alpha=0.6,
#     label="observed\nHay et al.",marker="x",markevery=10,markeredgewidth=2,markersize=5)

#### calculated from GMT ####

tslr = all_contributions.sum(axis="contribution")
tslr -= tslr[:,1986:2005].mean(axis="time")

upper_perc = da.DimArray(np.percentile(tslr,95,axis=0),
    axes=tslr.time,dims="time")#[obs_period.searchsorted(plot_period)]
lower_perc = da.DimArray(np.percentile(tslr,5,axis=0),
    axes=tslr.time,dims="time")#[obs_period.searchsorted(plot_period)]
median     = da.DimArray(np.percentile(tslr,50,axis=0),
    axes=tslr.time,dims="time")#[obs_period.searchsorted(plot_period)]

# upper_perc -= median[1986:2005].mean(axis="time")
# lower_perc -= median[1986:2005].mean(axis="time")
# median -= median[1986:2005].mean(axis="time")

# total_slr_calculated -= tslr - tslr[:,1986:2005].mean(axis="time")

# for i,slr_calculated in enumerate(total_slr_calculated):
#     lbl = "reconstructed" if i==0 else ""
ax.plot(tslr.time,median*1e3,lw=2,color="#c98345",
        label="reconstructed",alpha=1., )
ax.fill_between(lower_perc.time,lower_perc*1e3,upper_perc*1e3,lw=.5,color="#c98345",alpha=.3)

# axb.plot(obs_period,gd.giss_temp[obs_period],lw=2,color="red",alpha=1.)
# ax.text(0.05,0.8,"total", transform=ax.transAxes,fontweight='bold',)



ax.set_xlim(obs_period[0],obs_period[-1]+1)
ax.set_xlabel("Time in years")
ax.set_ylabel("Sea level in mm")
# axb.set_ylabel("gmt")
l1 = ax.legend(ncol=1,loc="upper left")#,bbox_to_anchor=(0.35,1.0))
l1.draw_frame(0)
for l in l1.get_lines(): l.set_alpha(1)

# for ax in axs:
#   ax.set_ylim(ax.get_ylim())
plt.draw()
plt.show()

runname = __file__.split("/")[-1][0:-3]
figname = "../figures/"+runname+"2.png"
figname = "../figures/"+runname+"2.pdf"

plt.savefig(figname)
