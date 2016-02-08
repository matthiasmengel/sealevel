
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
lib_path = os.path.abspath('../src')
sys.path.append(lib_path)
import get_calibration_data as gd; reload(gd)
import dimarray as da
import collections


plt.rcParams['xtick.major.pad']  = 10
plt.rcParams['font.size']= 12
plt.rcParams['lines.markeredgewidth']=2
plt.rcParams['legend.fontsize']=6
plt.rcParams['figure.figsize'] = 8,8
plt.rcParams['figure.facecolor'] = "white"

church_observed_700m = gd.church_observed["thermexp"] - gd.church_observed["te_700_3000"] - gd.church_observed["te_below_3000"]
church_observed_700m_2000m = gd.church_observed["te_700_3000"] + gd.church_observed["te_below_3000"] - gd.purkey10_below2000m

zero = church_observed_700m_2000m.copy()
zero[:] = 0.


# plt.close("all")
contrib_upper700 = collections.OrderedDict([
("domingues08",gd.thermo_obs_domingues08),
("ishii09",gd.thermo_obs_ishii09),
("levitus12",gd.thermo_obs_levit12_700),
("church11",church_observed_700m)])


contrib_700_2000m = collections.OrderedDict([
("levitus12",gd.thermo_obs_levit12_2000 - gd.thermo_obs_levit12_700),
("church11",church_observed_700m_2000m)])

contrib_below2000m = collections.OrderedDict([
("zero",zero),
("purkey10",gd.purkey10_below2000m)])


def remix_observations(above700m,between700_2000m,below2000m):

    obs = collections.OrderedDict()
    for ab700 in above700m:
        for bet700_2000 in between700_2000m:
            for bel2000 in below2000m:
                lbl = ab700+"_"+bet700_2000+"_"+bel2000
                obs[lbl] = above700m[ab700] + between700_2000m[bet700_2000] + below2000m[bel2000]
                obs[lbl] = obs[lbl].dropna()
                # print lbl,obs[lbl].values
    return obs

cols    = [cm.Accent(np.float(k)/24) for k in np.arange(24)]

plt.clf()
ax  = plt.subplot(111)

obs_rmx = remix_observations(contrib_upper700,contrib_700_2000m,contrib_below2000m)

for i,lb in enumerate(obs_rmx):

    cn = obs_rmx[lb] - obs_rmx[lb][1986:2005].mean()
    ax.plot(cn.time,cn,label=lb,lw=2,color=cols[i])


plt.grid()
plt.legend(loc=0,ncol=3)
plt.draw()
plt.show()



runname = __file__.split("/")[-1][0:-3]
plt.savefig("../figures/"+runname+".png")
