
""" matthias.mengel@pik
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
import get_data as gd; reload(gd)
import contributor_functions as cf; reload(cf)
import dimarray as da
import collections


plt.rcParams['xtick.major.pad']  = 10
plt.rcParams['font.size']= 12
plt.rcParams['lines.markeredgewidth']=2
plt.rcParams['legend.fontsize']=12
plt.rcParams['figure.figsize'] = 8,10
plt.rcParams['figure.facecolor'] = "white"
runname = __file__.split("/")[-1][0:-3]


cols = ["red","green","blue","black"]

nat = -gd.marzeion_gic_nat_mb.mean(axis="model")
full = -gd.marzeion_gic_full_mb.mean(axis="model")
anth = full - nat


plt.clf()
ax = plt.subplot(211)


ax.plot(full.time,full,label="FULL")
ax.plot(nat.time,nat,label="NAT")
ax.plot(anth.time,anth,label="ANTH")

# baseyear = 1870

# nat = gd.marzeion_gic_nat_up.mean(axis="model")
# full = gd.marzeion_gic_full_up.mean(axis="model")
# full -= full[baseyear].mean()
# nat -= nat[baseyear].mean()
# ax.plot(full.time,full,label="FULL")
# ax.plot(nat.time,nat,label="NAT")
# ax.plot(nat.time,full-nat,label="ANTH")
# ax.legend(loc=0)
anth_fraction = (full-nat)/full

anth_frac_nonan = anth_fraction.dropna()[1880:]
fn = np.polyfit(anth_frac_nonan.time, anth_frac_nonan, 1)
lin_fit = np.poly1d(fn)
ax2 = plt.subplot(212)

# full_old = gd.full_gic_contrib - gd.full_gic_contrib[baseyear]
# anth_old = gd.anth_gic_contrib - gd.anth_gic_contrib[baseyear]
# anth_frac_old = anth_old/full_old
ax2.plot(full.time,anth_fraction,label="calculated anth fraction\nfrom updated data")
# ax2.plot(anth_frac_old.time,anth_frac_old,label="calculated anth fraction from\nMarzeion Science 2014, Fig 1d")
ax2.plot(gd.anthropogenic_fraction.time,gd.anthropogenic_fraction,label="linear fit, used so far.")
ax2.plot(anth_frac_nonan.time,lin_fit(anth_frac_nonan.time),label="new linear fit")
ax2.plot(gd.frac_gic_contrib.time,gd.frac_gic_contrib,label="from Marzeion Science 2014, Fig 1c")
# ax2.legend(loc=0)

ax2.set_xlim(ax.get_xlim())
ax2.set_ylim(-0.1,1.1)
ax2.set_xlabel("time in years")
ax2.set_ylabel("anthropogenic fraction")
# ax.set_ylabel("SLE contribution in mm")

# ax.set_title("baseyear "+str(baseyear),fontsize=20)

# for i,obs_gic in enumerate(gd.gic_observations):
#     print "####",obs_gic
#     sl_observation = gd.gic_observations[obs_gic]
#     sl_observation -= sl_observation[2000]+0.1
#     plt.plot(sl_observation.time,sl_observation-sl_observation[1961],color=cols[i])
#     anthro = sl_observation*gd.anthropogenic_fraction
#     anthro -= anthro[1961]
#     plt.plot(sl_observation.time,anthro[sl_observation.time],"--",color=cols[i])


# plt.grid()
# plt.legend(loc=0,ncol=3)
# plt.savefig("../figures/"+runname+"_3.png")

plt.draw()
plt.show()
