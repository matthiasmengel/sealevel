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

import os
import numpy as np
import cPickle as pickle
import pandas as pd

import settings
reload(settings)
import src.calibration as calibration
reload(calibration)
import src.contributor_functions as cf
reload(cf)
import src.get_calibration_data as gcd
reload(gcd)
import src.calib_settings as cs
reload(cs)
import src.get_gmt_data as ggd
reload(ggd)
import src.sealevel as sl
reload(sl)


gmt = ggd.giss_temp

##### Thermal expansion #####

if "thermexp" in settings.calibrate_these:

    # This is the basis of Fig 10.34 of AR4, and in (c) they have pretty much levelled off by year 3000
    # Bern2D-CC 0.488; CLIMBER-2 0.458; CLIMBER-3a 0.200; MIT 0.214; MoBidiC
    # 0.626; UCL 0.386
    alpha_te = np.array([0.488, 0.458, 0.200, 0.214, 0.626, 0.386])
    sl_contributor = cf.thermal_expansion

    te_params = pd.DataFrame(index=pd.MultiIndex.from_product([cs.thermexp_observations.keys(),alpha_te],
        names=["observation","independent_param"]), columns=["dependent_param"])

    for obs_te in cs.thermexp_observations:
        print obs_te

        observation_period = cs.observation_period["thermexp"][obs_te]
        temp_anomaly_year = cs.temp_anomaly_year["thermexp"][obs_te]
        sl_observation = cs.thermexp_observations[obs_te]
        calib = calibration.Calibration(gmt, "thermexp",
                                        sl_observation, alpha_te, sl_contributor, observation_period, temp_anomaly_year)
        ip,dp = calib.calibrate()
        # ip is same as alpha_te. TODO: remove
        te_params.loc[obs_te,:]= dp

    te_params.to_csv(os.path.join(settings.calibfolder, "thermexp.csv"))

##### Glaciers and ice caps #####

if "gic" in settings.calibrate_these:

    sl_contributor = cf.glaciers_and_icecaps
    gic_modelno = np.arange(len(sl.gic_equi_functions))

    # gic_anth_params = {}
    gic_anth_params = pd.DataFrame(index=pd.MultiIndex.from_product([cs.gic_observations.keys(),gic_modelno],
        names=["observation","independent_param"]), columns=["dependent_param"])

    for i, obs_gic in enumerate(cs.gic_observations):
        print "####", obs_gic
        observation_period = cs.observation_period["gic"][obs_gic]
        temp_anomaly_year = cs.temp_anomaly_year["gic"][obs_gic]
        sl_observation = cs.gic_observations[obs_gic]

        anthro_gic = (sl_observation.diff() *
                      gcd.anthropogenic_frac).dropna().cumsum()
        calib = calibration.Calibration(gmt, "gic", anthro_gic, gic_modelno,
                                        sl_contributor, observation_period, temp_anomaly_year)
        ip,dp = calib.calibrate()
        gic_anth_params.loc[obs_gic,:]= dp

    gic_anth_params.to_csv(os.path.join(settings.calibfolder, "gic.csv"))

##### Greenland SMB #####

if "gis_smb" in settings.calibrate_these:

    sl_contributor = cf.surfacemassbalance_gis
    # levermann 13 coefficients.
    gis_smb_coeff = np.arange(0.05, 0.22, 0.02)

    gis_smb_params = pd.DataFrame(index=pd.MultiIndex.from_product([cs.gis_smb_observations.keys(),gis_smb_coeff],
        names=["observation","independent_param"]), columns=["dependent_param"])

    for i, obs_gis in enumerate(cs.gis_smb_observations):

        gmt_gis = ggd.giss_temp
        observation_period = cs.observation_period["gis_smb"][obs_gis]
        temp_anomaly_year = cs.temp_anomaly_year["gis_smb"][obs_gis]
        # apply an offset to temperature levels for box & colgan
        if obs_gis == "box_colgan13":
            gmt_gis = ggd.giss_temp + sl.gis_colgan_temperature_offset

        print "####", obs_gis
        # print observation_period
        sl_observation = cs.gis_smb_observations[obs_gis]
        calib = calibration.Calibration(gmt_gis, "gis_smb", sl_observation,
                                        gis_smb_coeff, sl_contributor, observation_period, temp_anomaly_year)
        ip,dp = calib.calibrate()
        gis_smb_params.loc[obs_gis,:]= dp

    gis_smb_params.to_csv(os.path.join(settings.calibfolder, "gis_smb.csv"))


##### Greenland solid ice discharge #####

if "gis_sid" in settings.calibrate_these:

    sl_contributor = cf.solid_ice_discharge_gis
    # levermann 13 coefficients.
    gis_sid_coeff = np.arange(-0.9, -0.45, 0.05)  # alpha
    # gis_sid_params = {}

    gis_sid_params = pd.DataFrame(index=pd.MultiIndex.from_product([cs.gis_sid_observations.keys(),gis_sid_coeff],
        names=["observation","independent_param"]), columns=["dependent_param"])

    for i, obs_gis in enumerate(cs.gis_sid_observations):

        gmt_gis = ggd.giss_temp
        observation_period = cs.observation_period["gis_sid"][obs_gis]
        temp_anomaly_year = cs.temp_anomaly_year["gis_sid"][obs_gis]
        # apply an offset to temperature levels for box & colgan
        if obs_gis == "box_colgan13":
            gmt_gis = ggd.giss_temp + sl.gis_colgan_temperature_offset

        print "####", obs_gis
        # print observation_period
        sl_observation = cs.gis_sid_observations[obs_gis]
        calib = calibration.Calibration(gmt_gis, "gis_sid", sl_observation,
                                        gis_sid_coeff, sl_contributor, observation_period, temp_anomaly_year)
        ip,dp = calib.calibrate()
        gis_sid_params.loc[obs_gis,:]= dp

    gis_sid_params.to_csv(os.path.join(settings.calibfolder, "gis_sid.csv"))


##### Antarctica solid ice discharge #####

if "ant_sid" in settings.calibrate_these:

    sl_contributor = cf.solid_ice_discharge_ais
    # Levermann et al. PNAS (2013) coefficients.
    ant_sid_coeff = np.arange(1.0, 1.55, 0.05)

    ant_sid_params = pd.DataFrame(index=pd.MultiIndex.from_product([cs.ant_sid_observations.keys(),ant_sid_coeff],
        names=["observation","independent_param"]), columns=["dependent_param"])

    for i, obs_ant_sid in enumerate(cs.ant_sid_observations):

        observation_period = cs.observation_period["ant_sid"][obs_ant_sid]
        temp_anomaly_year = cs.temp_anomaly_year["ant_sid"][obs_ant_sid]

        print "####", obs_ant_sid
        # print observation_period
        sl_observation = cs.ant_sid_observations[obs_ant_sid]
        calib = calibration.Calibration(gmt, obs_ant_sid, sl_observation,
                                        ant_sid_coeff, sl_contributor, observation_period, temp_anomaly_year)
        ip,dp = calib.calibrate()
        ant_sid_params.loc[obs_ant_sid,:]= dp

    ant_sid_params.to_csv(os.path.join(settings.calibfolder, "ant_sid.csv"))


##### Antarctica SMB #####

""" Antarcita SMB, we only pass the precipitation scaling parameters
 for projections, no real calibration. """

if "ant_smb" in settings.calibrate_these:

    # TODO: update later or remove. We do not really do a calibration here.

    sl_contributor = cf.surfacemassbalance_ais
    # inferred from Ligtenberg 13, in m/yr/K
    ais_prec_scaling = np.arange(2., 6.5, .5) * 1e-3

    # gis_sid_params = pd.DataFrame(index=pd.MultiIndex.from_product([cs.ant_sid_observations.keys(),ais_prec_scaling],
    #     names=["observation","independent_param"]), columns=["dependent_param"])

    # observation_period = "dummy"
    # temp_anomaly_year = "dummy"
    # sl_observation = "dummy"
    # calib = calibration.Calibration(gmt, "ant_smb", sl_observation,
    #                                 ais_prec_scaling, sl_contributor, observation_period, temp_anomaly_year)

    # ais_smb_params = {"ligtenberg13": calib}

    # outfile = os.path.join(settings.calibfolder, "ant_smb.pkl")
    # pickle.dump({"params": ais_smb_params}, open(outfile, "wb"), protocol=2)
