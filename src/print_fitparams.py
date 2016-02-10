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

import numpy as np
# import get_data as gd; reload(gd)
import cPickle as pickle

contrib_ids = ["thermexp", "gic", "gis_sid", "gis_smb", "ant_sid", "ant_smb"]
rcpcoldict = {"RCP3PD": "#3b4ba8", "RCP45": "yellow", "RCP85": "#ee1d23"}
labels = {"gic": "Mountain glaciers", "thermexp": "Thermal expansion", "gis_sid": "Greenland SID",
          "gis_smb": "Greenland SMB", "ant_sid": "Antarctica SID", "ant_smb": "Antarctica SMB"}


def rd(numb):
    return '%s' % float('%.3g' % numb)


def minmax(array):
    return rd(np.min(array)) + " to " + rd(np.max(array))

for name in contrib_ids:
    print labels[name],
    calibdata = pickle.load(open("../calibrationdata/" + name + ".pkl", "rb"))
    calibrated_param = np.array([])
    for obs in calibdata["params"]:
        # print obs,
        cal_par_obs = calibdata["params"][obs].calibrated_parameter
        calibrated_param = np.append(calibrated_param, cal_par_obs)

    print minmax(calibdata["params"][obs].commitment_parameter),
    print minmax(calibrated_param)


# contrib_names = ["Thermal expansion","Mountain glaciers (anth)","Greenland SMB", "Greenland SID", "Antarctica SMB", "Antarctica SID"]

# classname = {
# "Mountain glaciers (anth)":"glaciers_and_icecaps",
# "Thermal expansion":"thermal_expansion",
# "Greenland SID":"solid_ice_discharge_gis",
# "Greenland SMB":"surfacemassbalance_gis",
# "Antarctica SID":"solid_ice_discharge_ais",
# "Antarctica SMB":"surfacemassbalance_ais"}

# indpnames =  {
# "Mountain glaciers (anth)":"model number",
# "Thermal expansion":"alpha te",
# "Greenland SID":"alpha sid gis",
# "Greenland SMB":"gis smb coefficient",
# "Antarctica SID":"alpha sid ais",
# "Antarctica SMB":"ais precipitation factor"}

# depnames =  {
# "Mountain glaciers (anth)":"tau gic",
# "Thermal expansion":"tau te",
# "Greenland SID":"prefactor sid gis",
# "Greenland SMB":"tau smb gis",
# "Antarctica SID":"tau sid ais",
# "Antarctica SMB":"dummy, no calibration here"}


# def rd(numb): return '%s' % float('%.3g' % numb)

# # write contributions to file
# # width = 10
# fitparams_file = open('../data_python/fitparams.txt', 'w')
# fitparams_file.write('{0:30} {1:20} {2:20}\n'.format("contribution","commitment factor","calibrated parameter"))
# # for scen in ["RCP3PD","RCP45","RCP85"]
# for name in contrib_names:
#   calib_data  = np.load("../data_python/"+classname[name]+gd.small_iceage_offset_str+".npz")
#   ind_params  = calib_data["independent_param"]
#   dep_params  = calib_data["dependent_param"]

#   # ct = single_contribs[name]
#   l  = [name,rd(ind_params[0]),rd(ind_params[-1]),rd(np.min(dep_params)),rd(np.max(dep_params))]
#   fitparams_file.write('{0:30} [{1:5} to {2:6}]    [{3:8} to {4:7}]\n'.format(l[0],l[1],l[2],l[3],l[4]))

# # ct = total_contribs
# # l  = ["total",rd(ct["RCP3PD"][0]),rd(ct["RCP3PD"][1]),rd(ct["RCP3PD"][2]),rd(ct["RCP45"][0]),rd(ct["RCP45"][1]),rd(ct["RCP45"][2]),rd(ct["RCP85"][0]),rd(ct["RCP85"][1]),rd(ct["RCP85"][2])]
# # fitparams_file.write('{0:30} {1:7} ({2:7}-{3:7}) {4:7} ({5:7}-{6:7}) {7:7} ({8:7}-{9:7})\n'.format(l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7],l[8],l[9]))
# # fitparams_file.write("median, lower perc, upper perc\n")
# # fitparams_file.write("## single contributions ##\n")
# fitparams_file.close()


# for i,name in enumerate(contrib_names):
#   print "#####",name, "#####"
#   print "independent parameter:",indpnames[name]
#   print "dependent parameter:",depnames[name]
#   print indpnames[name],"|",depnames[name]
#   ## load data
#   calib_data = np.load("../data_python/"+classname[name]+gd.small_iceage_offset_str+".npz")
#   ind_params  = calib_data["independent_param"]
#   dep_params  = calib_data["dependent_param"]

#   for j,indp in enumerate(ind_params):
#     print indp,"|",dep_params[j]
