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

project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
inputdatadir = os.path.join(project_dir, "data/")

######## parameters that need to be known for calibrations and projection

# add temperature offset to be used with the box & colgan 2013 data to fit past observations
# similar offset is used by e.g. Rahmstorf 2007, Science.
gis_colgan_temperature_offset = 0.5


######## sea level projection using for Monte-Carlo sampling ########

def project(gmt, proj_period, calibdata, temp_anomaly_year, sl_contributor,
            sample_number):
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
        driving_temperature = gmt[proj_period, gmt_choice]
    except IndexError:
        # this is the case if single gmt is supplied
        gmt_choice = 0
        driving_temperature = gmt[proj_period]

    # print contrib_name, temp_anomaly_year
    # use one of the observational dataset
    obs_choice = np.random.choice(calibdata.index.unique())
    params = calibdata.loc[obs_choice]
    # temp_anomaly_year = params.temp_anomaly_year

    if obs_choice == "box_colgan13":
        driving_temperature += gis_colgan_temperature_offset

    # choose a random parameter set
    tuple_choice = np.random.randint(len(params))
    independent_param = params["independent_param"][tuple_choice]
    dependent_param = params["dependent_param"][tuple_choice]
    # print "tuple",independent_param,dependent_param

    contributor = sl_contributor(independent_param, temp_anomaly_year[obs_choice])
    contrib = contributor.calc_contribution(
        driving_temperature, dependent_param)
    # print contrib
    return [contrib, gmt_choice, obs_choice, independent_param, dependent_param]
