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


""" matthias.mengel@pik
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
import get_data as gd; reload(gd)
import dimarray as da
import optparse
from mpl_toolkits.axes_grid1 import make_axes_locatable
import timeit
import cPickle as pickle

# make use of -n or -g flags with this script
parser = optparse.OptionParser()
parser.add_option('-n', '--realizations',  type="int",  help='get number of realizations')
parser.add_option('-g', help='get new data', dest='get_data', default=False, action='store_true')
(opts, args) = parser.parse_args()

contrib_ids = ["thermexp","gic","gis_sid","gis_smb",  "ant_sid", "ant_smb"]
rcpcoldict = {"RCP3PD":"#3b4ba8","RCP45":"yellow","RCP85":"#ee1d23"}
labels = {"gic":"Mountain glaciers","thermexp":"Thermal expansion","gis_sid":"Greenland SID",
"gis_smb":"Greenland SMB","ant_sid":"Antarctica SID","ant_smb":"Antarctica SMB"}

if opts.get_data:
  realizations = 10000
  realizations = realizations  if opts.realizations==None else opts.realizations
  print "use",realizations,"realizations"
  projection_data = pickle.load(open("../projectiondata/projected_contributons_"+str(realizations)+".pkl","rb"))


def rd(array,percentile,years):

    if type(years) == int:
        numb = array[years]
    else:
        numb = array[years].mean(axis="time")

    perc = np.percentile(numb,percentile)
    # in mm
    return '%s' % float('%.3g' % perc)

printyr = 2100
#printyr = np.arange(2081,2101,1)
print "## years",printyr
# printdata = {}

# print "#### Projections ####"
print "RCP3PD","RCP45","RCP85"
for i,name in enumerate(contrib_ids):
  # printdata[name] = {}
  print labels[name],":",
  for k,scen in enumerate(["RCP3PD","RCP45","RCP85",]):

    contrib    = projection_data[scen][name]*1.e3
    # plot_period_indices = contrib.time.searchsorted(plot_period)
    # anomaly to 1986-2005
    contrib    = contrib - contrib[1986:2005,:].mean(axis="time")
    # pandas handles NaNs within percentile calculus
    # contrib     = contrib.T.to_pandas()
    # percentiles = da.from_pandas(contrib.quantile([0.05,0.5,0.95])*1.e3,
    #   dims=["percentile","time"])

    if scen != "RCP85":
      print rd(contrib,50,printyr),"(",rd(contrib,5,printyr),"to",rd(contrib,95,printyr),")",
    else:
      print rd(contrib,50,printyr),"(",rd(contrib,5,printyr),"to",rd(contrib,95,printyr),")"

    # printdata[name][scen]=percentiles
    # single_contribs[name][scen] = np.array([median[-1],lower_perc[-1],upper_perc[-1]])



#### total slr; sum up contributions first ####


def rdd(numb):
  # in mm
  return '%s' % float('%.4g' % numb)

print "total",":",
for k,scen in enumerate(["RCP3PD","RCP45","RCP85",]):
    total_slr = da.zeros_like(projection_data[scen]["thermexp"])
    for i,name in enumerate(contrib_ids):
        # sum up all contributions
        single_contrib = projection_data[scen][name]
        # if nans occur, clear these contributions
        # single_contrib[np.isnan(single_contrib)] = 0.
        total_slr += single_contrib

    di_total_slr = da.DimArray(total_slr,dims=["time","runnumber"])
    di_total_slr -= di_total_slr[1986:2005,:].mean(axis="time")
    if type(printyr) == int:
      mn = di_total_slr[printyr,:]*1.e3
    else:
      mn = di_total_slr[printyr,:].mean(axis="time")*1.e3
    # mn = di_total_slr[2081:2100,:].mean(axis="time")
    low = np.percentile(mn,5)
    med = np.percentile(mn,50)
    upp = np.percentile(mn,95)
    print rdd(med),"(",rdd(low),"to",rdd(upp),")",


# pdata = da.DimArray(printdata,dims=["contribution","scenario","percentile","time"])



# for scen in ["RCP3PD","RCP45","RCP85"]:
#     total = printdata[:,scen,:].sum(axis="contribution")
#     print rdd(total['median']),"(",rdd(total['0.05']),"to",rdd(total['0.95']),")",

