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

""" Calibrate the sea level contriubtions as in
    M. Mengel et al.
    Future sea-level rise constrained by observations and long-term commitment
    PNAS (2016)
    (C) Matthias Mengel working at Potsdam Institute for Climate Impact Research

"""

import calibration; reload(calibration)
import contributor_functions as cf; reload(cf)
import get_calibration_data as gcd; reload(gcd)
import calib_settings as cs; reload(cs)
import get_gmt_data as ggd; reload(ggd)
import sealevel as sl; reload(sl)
import numpy as np
import collections
import dimarray as da
import cPickle as pickle

contrib_ids = ["thermexp","gic","gis_smb", "gis_sid", "ant_smb", "ant_sid"]

calibrate_these = ["thermexp","gic","gis_smb", "gis_sid", "ant_smb", "ant_sid"]

gmt = ggd.giss_temp

calibdatadir = "../data/calibration/"

##### Thermal expansion #####

if "thermexp" in calibrate_these:

    # sl_observation = gcd.church_observed["Thermal expansion"]

    # This is the basis of Fig 10.34 of AR4, and in (c) they have pretty much levelled off by year 3000
    # Bern2D-CC 0.488; CLIMBER-2 0.458; CLIMBER-3a 0.200; MIT 0.214; MoBidiC 0.626; UCL 0.386
    alpha_te = np.array([0.488, 0.458, 0.200, 0.214, 0.626, 0.386])
    sl_contributor = cf.thermal_expansion

    te_params = {}
    # te_condensed = np.array([])

    for obs_te in cs.thermexp_observations:
        print obs_te

        observation_period = cs.observation_period["thermexp"][obs_te]
        temp_anomaly_year = cs.temp_anomaly_year["thermexp"][obs_te]
        sl_observation = cs.thermexp_observations[obs_te]
        calib = calibration.Calibration(gmt, "thermexp",
                    sl_observation, alpha_te, sl_contributor, observation_period, temp_anomaly_year)
        calib.calibrate()
        te_params[obs_te] = calib
        # te_condensed = np.append(te_condensed,te_params[obs_te][1])

    # first axis along alpha_te
    # te_condensed = te_condensed.reshape(len(alpha_te),-1)

    outfile = calibdatadir+"thermexp.pkl"
    pickle.dump({"params":te_params}, open(outfile,"wb"),protocol=2)

##### Glaciers and ice caps #####

if "gic" in calibrate_these:

    sl_contributor = cf.glaciers_and_icecaps
    gic_modelno    = np.arange(len(sl.gic_equi_functions))

    gic_anth_params = {}

    for i,obs_gic in enumerate(cs.gic_observations):
        print "####",obs_gic
        observation_period = cs.observation_period["gic"][obs_gic]
        temp_anomaly_year = cs.temp_anomaly_year["gic"][obs_gic]
        sl_observation = cs.gic_observations[obs_gic]

        anthro_gic = (sl_observation.diff()*gcd.anthropogenic_frac).dropna().cumsum()
        calib = calibration.Calibration(gmt, "gic", anthro_gic, gic_modelno,
            sl_contributor,observation_period, temp_anomaly_year)
        calib.calibrate()
        gic_anth_params[obs_gic] = calib

    outfile = calibdatadir+"gic.pkl"
    pickle.dump({"params":gic_anth_params}
        ,open(outfile,"wb"),protocol=2)


##### Greenland SMB #####

if "gis_smb" in calibrate_these:

    sl_contributor = cf.surfacemassbalance_gis
    # levermann 13 coefficients.
    gis_smb_coeff = np.arange(0.05,0.22,0.02)
    gis_smb_params = {}

    for i,obs_gis in enumerate(cs.gis_smb_observations):

        gmt_gis = ggd.giss_temp
        observation_period = cs.observation_period["gis_smb"][obs_gis]
        temp_anomaly_year = cs.temp_anomaly_year["gis_smb"][obs_gis]
        ## apply an offset to temperature levels for box & colgan
        if obs_gis == "box_colgan13":
            gmt_gis = ggd.giss_temp + sl.gis_colgan_temperature_offset
            # observation_period = np.arange(1880,2013)
            # temp_anomaly_year  = None

        print "####",obs_gis
        # print observation_period
        sl_observation = cs.gis_smb_observations[obs_gis]
        calib = calibration.Calibration(gmt_gis, "gis_smb", sl_observation,
            gis_smb_coeff, sl_contributor, observation_period, temp_anomaly_year)
        calib.calibrate()
        gis_smb_params[obs_gis] = calib

    outfile = calibdatadir+"gis_smb.pkl"
    # only params will be used for projections
    pickle.dump({"params":gis_smb_params}
        ,open(outfile,"wb"),protocol=2)


##### Greenland solid ice discharge #####

if "gis_sid" in calibrate_these:

    sl_contributor = cf.solid_ice_discharge_gis
    # levermann 13 coefficients.
    gis_sid_coeff = np.arange(-0.9,-0.45,0.05) # alpha
    gis_sid_params = {}

    for i,obs_gis in enumerate(cs.gis_sid_observations):

        gmt_gis = ggd.giss_temp
        observation_period = cs.observation_period["gis_sid"][obs_gis]
        temp_anomaly_year = cs.temp_anomaly_year["gis_sid"][obs_gis]
        ## apply an offset to temperature levels for box & colgan
        if obs_gis == "box_colgan13":
            gmt_gis = ggd.giss_temp + sl.gis_colgan_temperature_offset

        print "####",obs_gis
        # print observation_period
        sl_observation = cs.gis_sid_observations[obs_gis]
        calib = calibration.Calibration(gmt_gis, "gis_sid", sl_observation,
            gis_sid_coeff, sl_contributor, observation_period, temp_anomaly_year)
        calib.calibrate()
        gis_sid_params[obs_gis] = calib

    outfile = calibdatadir+"gis_sid.pkl"
    # only params will be used for projections
    pickle.dump({"params":gis_sid_params}
        ,open(outfile,"wb"),protocol=2)


##### Antarctica solid ice discharge #####

if "ant_sid" in calibrate_these:

    sl_contributor = cf.solid_ice_discharge_ais
    # levermann 13 coefficients.
    ant_sid_coeff = np.arange(1.0,1.55,0.05)
    ant_sid_params = {}

    for i,obs_ant_sid in enumerate(cs.ant_sid_observations):

        observation_period = cs.observation_period["ant_sid"][obs_ant_sid]
        temp_anomaly_year = cs.temp_anomaly_year["ant_sid"][obs_ant_sid]

        print "####",obs_ant_sid
        # print observation_period
        sl_observation = cs.ant_sid_observations[obs_ant_sid]
        calib = calibration.Calibration(gmt, obs_ant_sid, sl_observation,
            ant_sid_coeff, sl_contributor, observation_period, temp_anomaly_year)
        calib.calibrate()
        ant_sid_params[obs_ant_sid] = calib

    outfile = calibdatadir+"ant_sid.pkl"
    # only params will be used for projections
    pickle.dump({"params":ant_sid_params}
        ,open(outfile,"wb"),protocol=2)


##### Antarctica SMB #####

""" Antarcita SMB, we only pass the precipitation scaling parameters
 for projections, no real calibration. """

if "ant_smb" in calibrate_these:

    sl_contributor = cf.surfacemassbalance_ais
    # inferred from Ligtenberg 13, in m/yr/K
    ais_prec_scaling=np.arange(2.,6.5,.5)*1e-3
    ais_smb_params = {}

    observation_period = "dummy"
    temp_anomaly_year = "dummy"
    sl_observation = "dummy"
    calib = calibration.Calibration(gmt, "ant_smb", sl_observation,
        ais_prec_scaling, sl_contributor, observation_period, temp_anomaly_year)

    ais_smb_params = {"ligtenberg13":calib}

    outfile = calibdatadir+"ant_smb.pkl"
    pickle.dump({"params":ais_smb_params}
        ,open(outfile,"wb"),protocol=2)
