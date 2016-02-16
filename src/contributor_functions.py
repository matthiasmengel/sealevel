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


""" Sea level equilibrium responses and respective transient sea level functions.
    Equilibrium responses as in equ. 2,3,5 and 6 in M. Mengel et al., PNAS (2016).
    The transient response is constructed via equ. 1.
"""

import os
import numpy as np
import dimarray as da
import sealevel as sl
reload(sl)

project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
inputdatadir = os.path.join(project_dir, "inputdata/")

dtime = 1
area_ocean = 3.61132e14  # % m^2
area_antarctica = 0.14e14  # m^2


################################################

class thermal_expansion(object):

    """ Thermal expansion equilibrium sea level response and transient response.
        See M. Mengel et al., PNAS (2016), Materials and Methods and equ. 2. """

    def __init__(self, alpha, temp_anomaly_year=None):

        self.dtime = dtime
        self.alpha = alpha

    def te_equilibrium_sl(self, delta_temp):
        """ equilibrium response
            alpha is the equilibrium sl sensitivity in m/K """

        return self.alpha * delta_temp

    def calc_contribution(self, delta_gmt, tau):
        """ transient response, see equ. 1 in
            M. Mengel et al., PNAS (2016) """

        delta_gmt = delta_gmt.values
        sl_contrib = np.zeros_like(delta_gmt, dtype="float")
        sl_rate = 0.0

        for t in np.arange(1, len(delta_gmt), 1):

            sl_rate = (self.te_equilibrium_sl(
                delta_gmt[t - 1]) - sl_contrib[t - 1]) / tau
            # sl_rate = np.maximum(sl_rate,0.)
            sl_contrib[t] = sl_contrib[t - 1] + self.dtime * sl_rate

        return sl_contrib


################################################

class glaciers_and_icecaps(object):

    """ Glaciers and ice caps equilibrium sea level response and transient response.
        See M. Mengel et al., PNAS (2016), Materials and Methods and supplementary
        material. """

    def __init__(self, modelno, temp_anomaly_year=None):

        self.dtime = dtime
        self.modelno = modelno
        self.temp_anomaly_year = temp_anomaly_year

    def gic_equilibrium_sl(self, delta_temp):
        """ glacier equilbrium response via custom functions, see
            M. Mengel et al., PNAS (2016), equation in supplementary material. """

        return sl.gic_equi_functions[self.modelno](delta_temp)

    def calc_contribution(self, delta_gmt, tau):
        """ transient response, see equ. 1 in
            M. Mengel et al., PNAS (2016) """

        if self.temp_anomaly_year is not None:
            # do not let temperature drive ice loss before first year of SL
            # observation
            delta_gmt -= delta_gmt[self.temp_anomaly_year]
            delta_gmt[delta_gmt.time < self.temp_anomaly_year] = 0.

        delta_gmt = delta_gmt.values
        sl_contrib = np.zeros_like(delta_gmt, dtype="float")
        sl_rate = 0.0

        for t in np.arange(1, len(delta_gmt), 1):
            sl_rate = (self.gic_equilibrium_sl(
                delta_gmt[t]) - sl_contrib[t - 1]) / tau
            # sl_rate = np.maximum(sl_rate,0.)
            sl_contrib[t] = sl_contrib[t - 1] + self.dtime * sl_rate

        return sl_contrib


################################################

class surfacemassbalance_gis(object):

    """ Greenland ice sheet surface mass balaance
        equilibrium sea level response and transient response.
        See M. Mengel et al., PNAS (2016), Materials and Methods
    """

    def __init__(self, smb_coeff, temp_anomaly_year=None):

        self.dtime = dtime
        self.smb_coeff = smb_coeff
        self.temp_anomaly_year = temp_anomaly_year

    def equilibrium_sl(self, delta_temp):
        """ equilibrium response as in equ.3.
            alpha is the equilibrium sl sensitivity in m/K^2 """

        return self.smb_coeff * np.sign(delta_temp) * delta_temp**2

    def calc_contribution(self, delta_gmt, tau):
        """ transient response, see equ. 1 in
            M. Mengel et al., PNAS (2016) """

        if self.temp_anomaly_year is not None:
            # do not let temperature drive ice loss before first year of SL
            # observation
            delta_gmt -= delta_gmt[self.temp_anomaly_year]
            delta_gmt[delta_gmt.time < self.temp_anomaly_year] = 0.

        delta_gmt = delta_gmt.values
        sl_contrib = np.zeros_like(delta_gmt, dtype="float")
        sl_rate = 0.0

        for t in np.arange(1, len(delta_gmt), 1):

            sl_rate = (self.equilibrium_sl(
                delta_gmt[t]) - sl_contrib[t - 1]) / tau
            sl_contrib[t] = sl_contrib[t - 1] + self.dtime * sl_rate

        return sl_contrib


################################################

class surfacemassbalance_ais(object):

    """ Antarctic ice sheet surface mass balance. This is not calibrated and
        we do not apply the pursuit curve method. We use scaling through the
        Clausius Clapeyron temperature dependence of water carrying capacity of air."""

    def __init__(self, smb_coeff, temp_anomaly_year=None):

        self.dtime = dtime
        self.smb_coeff = smb_coeff
        # 25 per cent will be lost due to increased discharge.
        self.winkelmann_factor = 0.75

    def calc_contribution(self, delta_gmt, dummy):
        """ Only scaling, no pursuit curve approach. See equ.5 in
            Materials and Methods."""

        delta_gmt = delta_gmt.values
        sl_contrib = np.zeros_like(delta_gmt, dtype="float")
        sl_rate = 0.0

        for t in np.arange(1, len(delta_gmt), 1):

            sl_rate = self.winkelmann_factor * \
                self.smb_coeff * delta_gmt[t - 1]
            sl_contrib[t] = sl_contrib[t - 1] + self.dtime * sl_rate

        # to m slr
        return -sl_contrib * area_antarctica / area_ocean


################################################

class solid_ice_discharge_gis(object):

    """ Solid ice discharge from the Greenland ice sheet. There is no equilibrium
        estimate available for solid ice discharge. We therefore use a response
        function approach, see equ. 4 in Materials and Methods. """

    def __init__(self, alpha, temp_anomaly_year=None):

        self.alpha = alpha
        self.dtime = dtime
        # tau in years; does only play a role in saturation, so not until 2100
        self.tau = 40.
        # self.calibrate = calibrate
        self.temp_anomaly_year = temp_anomaly_year

    def calc_contribution(self, temperature, prefactor):
        """ No pursuit curve, but response function approach. """

        oceantemp = temperature

        if self.temp_anomaly_year is not None:
            # do not let temperature drive ice loss before first year of SL
            # observation
            oceantemp -= oceantemp[self.temp_anomaly_year]
            oceantemp[oceantemp.time < self.temp_anomaly_year] = 0.

        otemp = oceantemp.values
        # print otemp
        discharge = np.zeros_like(otemp, dtype="float")

        for t in np.arange(1, len(otemp), 1):
            discharge_rate = 0
            # integrate over t to yield slr, see winkelmann response function
            # paper p. 2581
            tp = np.arange(t)
            discharge_rate = prefactor * \
                np.trapz(otemp[tp] * ((t - tp))**self.alpha, dx=self.dtime)
            discharge_rate = np.maximum(discharge_rate, 0.0)
            discharge[t] = discharge[t - 1] + self.dtime * discharge_rate

        return discharge


################################################

class solid_ice_discharge_ais(object):

    """ Antarctic ice sheet solid ice discharge,
        equilibrium sea level response and transient response.
        See M. Mengel et al., PNAS (2016), Materials and Methods
    """

    def __init__(self, alpha, temp_anomaly_year=None):

        self.alpha = alpha
        self.dtime = dtime
        self.temp_anomaly_year = temp_anomaly_year

    def equilibrium_sl(self, delta_temp):
        """ see equ. 6 in Materials and Methods. """
        return self.alpha * delta_temp

    def calc_contribution(self, temperature, tau):
        """ transient response, see equ. 1 in
            M. Mengel et al., PNAS (2016) """

        if self.temp_anomaly_year is not None:
            # do not let temperature drive ice loss before first year of SL
            # observation
            temperature -= temperature[self.temp_anomaly_year]
            temperature[temperature.time < self.temp_anomaly_year] = 0.

        temperature = temperature.values
        discharge = np.zeros_like(temperature, dtype="float")
        discharge_rate = 0

        for t in np.arange(1, len(temperature), 1):

            discharge_rate = (self.equilibrium_sl(
                temperature[t]) - discharge[t - 1]) / tau
            discharge_rate = np.maximum(discharge_rate, 0.0)
            discharge[t] = discharge[t - 1] + self.dtime * discharge_rate

        return discharge
