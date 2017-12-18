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
import pandas as pd
import dimarray as da
import sealevel.contributor_functions as cf
reload(cf)


project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
inputdatadir = os.path.join(project_dir, "data/")

######## parameters that need to be known for calibrations and projection

# add temperature offset to be used with the box & colgan 2013 data to fit past observations
# similar offset is used by e.g. Rahmstorf 2007, Science.
gis_colgan_temperature_offset = 0.5


######## sea level projection using for Monte-Carlo sampling ########

def project(gmt, proj_period, calibdata, temp_anomaly_year, sl_contributor,
            sample_number, contrib_name):
    """
    Monte Carlo sampling for slr contribution
    for a single global mean temperature (gmt) timeseries or
    an ensemble of gmt timeseries,
    for one random choice of observations obs_choice,
    and one random tuple of independent and dependent parameter.
    the contributor function (i.e. thermal expansion) is chosen through
    contrib_name.

    Parameters
    ----------
    gmt : single or ensemble of gmt timeseries
    proj_period : time period for which slr projection is done
    calibdata : calibration data for the several observations per component
    temp_anomaly_year: year in which global mean temperature passes zero,
                       depending on observation.
    sl_contributor: function to calculate transient sea level rise.
    sample_number : number for seed to be created to make sampling reproducible

    Returns
    -------
    contrib : timeseries of sea level contribution with length proj_period
    """

    np.random.seed(sample_number)

    try:
        gmt_ensemble_size = gmt.shape[1]
        gmt_choice = np.random.randint(gmt_ensemble_size)
        # print gmt_choice
        driving_temperature = gmt[:, gmt_choice]
    except IndexError:
        # this is the case if single gmt is supplied
        gmt_choice = 0
        driving_temperature = gmt

    # print contrib_name, temp_anomaly_year

    # convert index to str to avoid floating point issues
    calibdata.index = [str(i) for i in calibdata.index]
    # use one of the observational dataset
    obs_choice = np.random.choice(calibdata.index.unique())
    params_of_obs = calibdata.loc[obs_choice]
    # print params_of_obs
    # temp_anomaly_year = params.temp_anomaly_year

    if obs_choice == "box_colgan13":
        driving_temperature += gis_colgan_temperature_offset

        # for dp16, the different ensemble members are interpreted
        # as different observations, so selection already happened
        # above through obs_choice
    if contrib_name == "ant_dp16":
        params = params_of_obs
    else:
        # choose a random parameter set
        paramset_choice = np.random.randint(len(params_of_obs.index))
        # can be variable number of parameters per each observation
        params = params_of_obs.iloc[paramset_choice,:]
        # print "pp",params
    contributor = sl_contributor(params, temp_anomaly_year.loc[obs_choice][0])
    contrib = contributor.calc_contribution(
        driving_temperature,proj_period)
    # print contrib
    return [contrib, gmt_choice, obs_choice, params]



def project_slr(scen, gmt, settings):

    projection_data = {}

    temp_anomaly_years = pd.read_csv(os.path.join(
        settings.calibfolder, "temp_anomaly_years.csv"),index_col=[0,1])
    temp_anomaly_years = temp_anomaly_years.where(
        pd.notnull(temp_anomaly_years), None)

    for i, contrib_name in enumerate(settings.project_these):

        print "conribution", contrib_name

        realizations = np.arange(settings.nrealizations)
        calibdata = pd.read_csv(
            os.path.join(settings.calibfolder, contrib_name+".csv"),
            index_col=[0])

        temp_anomaly_year = temp_anomaly_years.loc[contrib_name]
        sl_contributor = cf.contributor_functions[contrib_name]

        proj = np.zeros([len(settings.proj_period), settings.nrealizations])

        for n in realizations:
            slr, gmt_n, obs_choice, params = project(
                gmt, settings.proj_period, calibdata, temp_anomaly_year,
                sl_contributor, n, contrib_name)
            proj[:, n] = slr

        pdata = da.DimArray(proj, axes=[settings.proj_period, realizations],
                            dims=["time", "runnumber"])

        projection_data[contrib_name] = pdata

    if not os.path.exists(settings.projected_slr_folder):
        os.makedirs(settings.projected_slr_folder)

    fname = "projected_slr_"+scen+"_n"+str(settings.nrealizations)+".nc"
    da.Dataset(projection_data).write_nc(os.path.join(
        settings.projected_slr_folder,fname))
    print "Sea level projection data written to"
    print settings.projected_slr_folder