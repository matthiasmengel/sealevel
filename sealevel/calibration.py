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
from scipy.optimize import curve_fit, leastsq
import dimarray as da


class Calibration(object):

    def __init__(self, global_mean_temperature, label,
                 sl_observation, commitment_parameter, sl_contributor,
                 observation_period=None, temp_anomaly_year=None):

        self.gmt = global_mean_temperature
        self.label = label
        # self.obs_period = observation_period
        self.sl_observation = sl_observation

        self.commitment_parameter = commitment_parameter
        self.sl_contributor = sl_contributor
        self.temp_anomaly_year = temp_anomaly_year

        if label == "ant_smb":
            # workaround for Antarctic SMB, we do not want to calibrate here.
            self.calibrated_parameter = np.zeros_like(
                self.commitment_parameter)
            pass
        else:
            # restrict calibration period to available observations.
            if observation_period is None:
                print "observation period none"
                cal_period_0 = np.max(
                    [self.gmt.time[0], sl_observation.time[0]])
                cal_period_1 = np.min(
                    [self.gmt.time[-1], sl_observation.time[-1]])

            else:
                # also restrict to provided observation period
                cal_period_0 = np.max([observation_period[0], self.gmt.time[0],
                                       sl_observation.time[0]])
                cal_period_1 = np.min([observation_period[-1], self.gmt.time[-1],
                                       sl_observation.time[-1]])

            self.cal_period = np.arange(cal_period_0, cal_period_1 + 1, 1)


    def calibrate(self):

        ## first guess for tau
        tau0 = 100.
        optimal_tau = np.zeros_like(self.commitment_parameter, dtype="float")
        for i, cparam in enumerate(self.commitment_parameter):

            cl = self.sl_contributor([cparam,"dummy"], self.temp_anomaly_year)

            def calc_residuals(param_to_calibrate):
                """ Only residuals in calibration period play a role for fit.
                    Residuals start always with zero in the first year of cal_period. """

                proj_period = self.gmt.time
                contrib = da.DimArray(
                    cl.calc_contribution(self.gmt, proj_period, param_to_calibrate),
                    axes=proj_period,
                    dims="time")
                residuals = da.DimArray(
                    np.zeros_like(
                        contrib.values),
                    axes=proj_period,
                    dims="time")

                residuals[self.cal_period] = self.sl_observation[
                    self.cal_period] - contrib[self.cal_period]
                residuals[self.cal_period] -= residuals[self.cal_period[0]]

                return residuals.values

            optimal_tau[i], pcov = leastsq(calc_residuals, tau0)

            print "cparam=", cparam, " tau=", optimal_tau[i]

        self.calibrated_parameter = optimal_tau
        return [self.commitment_parameter, optimal_tau]
