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
# import sealevel as sl
# reload(sl)

project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
inputdatadir = os.path.join(project_dir, "data")

dtime = 1
area_ocean = 3.61132e14  # % m^2
area_antarctica = 0.14e14  # m^2


################################################

## equilibrium sea level functions for mountain glaciers

def gic_equi_func(equi_temp, a, b):
    """ assume exponential functional form. This is justified, because
    1) the function saturates for high temperatures.
    2) the function passes zero for zero temperature, as wanted for
       anthropogenic glacier contribution.
    """

    return a * (1 - np.exp(-1 * b * equi_temp))


def func_creator(a, b):

    def func(equi_temp):
        """ same as gic_equi_func with explicit values for a and b as set by curve_fit. """
        return a * (1 - np.exp(-1 * b * equi_temp))

    return func


gic_equi_coeffs = np.loadtxt(os.path.join(inputdatadir,"glacier_equi",
                             "glacier_equi_coefficients.csv"))
gic_equi_functions = []
for coeffs in gic_equi_coeffs:
    # print coeffs
    gic_equi_functions.append(func_creator(*coeffs))


class contribution(object):


    """ General class for sea level contributors. """

    def __init__(self, parameters, temp_anomaly_year=None):

        self.dtime = dtime
        self.alpha = parameters[0]
        self.tau = parameters[1]
        self.temp_anomaly_year = temp_anomaly_year

    def calc_contribution(self, delta_gmt, proj_period, tau=None):
        """ transient response, see equ. 1 in
            M. Mengel et al., PNAS (2016) """

        # for projection, tau needs not to be explicitely provided
        # as argument. in calibration it has to.
        if tau is None: tau = self.tau

        if self.temp_anomaly_year is not None:
            # do not let temperature drive ice loss before first year of SL
            # observation
            delta_gmt -= delta_gmt[self.temp_anomaly_year]
            delta_gmt[delta_gmt.time < self.temp_anomaly_year] = 0.

        delta_gmt = delta_gmt[proj_period].values
        sl_contrib = np.zeros_like(delta_gmt, dtype="float")
        sl_rate = 0.0

        for t in np.arange(1, len(delta_gmt), 1):

            sl_rate = (self.equilibrium_sl(
                delta_gmt[t - 1]) - sl_contrib[t - 1]) / tau
            # sl_rate = np.maximum(sl_rate,0.)
            sl_contrib[t] = sl_contrib[t - 1] + self.dtime * sl_rate

        return sl_contrib



class thermal_expansion(contribution):

    """ Thermal expansion equilibrium sea level response and transient response.
        See M. Mengel et al., PNAS (2016), Materials and Methods and equ. 2. """

    def equilibrium_sl(self, delta_temp):
        """ equilibrium response
            alpha is the equilibrium sl sensitivity in m/K """

        return self.alpha * delta_temp


class glaciers_and_icecaps(contribution):

    """ Glaciers and ice caps equilibrium sea level response and transient response.
        See M. Mengel et al., PNAS (2016), Materials and Methods and supplementary
        material. """

    def __init__(self, parameters, temp_anomaly_year=None):

        self.dtime = dtime
        self.modelno = int(parameters[0])
        self.tau = parameters[1]
        self.temp_anomaly_year = temp_anomaly_year

    def equilibrium_sl(self, delta_temp):
        """ glacier equilbrium response via custom functions, see
            M. Mengel et al., PNAS (2016), equation in supplementary material. """

        return gic_equi_functions[self.modelno](delta_temp)


class surfacemassbalance_gis(contribution):

    """ Greenland ice sheet surface mass balaance
        equilibrium sea level response and transient response.
        See M. Mengel et al., PNAS (2016), Materials and Methods
    """

    def __init__(self, parameters, temp_anomaly_year=None):

        self.dtime = dtime
        self.smb_coeff = parameters[0]
        self.tau = parameters[1]
        self.temp_anomaly_year = temp_anomaly_year

    def equilibrium_sl(self, delta_temp):
        """ equilibrium response as in equ.3.
            alpha is the equilibrium sl sensitivity in m/K^2 """

        return self.smb_coeff * np.sign(delta_temp) * delta_temp**2


class surfacemassbalance_ais(contribution):

    """ Antarctic ice sheet surface mass balance. This is not calibrated and
        we do not apply the pursuit curve method. We use scaling through the
        Clausius Clapeyron temperature dependence of water carrying capacity of air."""

    def __init__(self, parameters, temp_anomaly_year=None):

        self.dtime = dtime
        self.smb_coeff = parameters[0]
        # 25 per cent will be lost due to increased discharge.
        self.winkelmann_factor = 0.75

    def calc_contribution(self, delta_gmt, proj_period):
        """ Only scaling, no pursuit curve approach. See equ.5 in
            Materials and Methods."""

        delta_gmt = delta_gmt[proj_period].values
        sl_contrib = np.zeros_like(delta_gmt, dtype="float")
        sl_rate = 0.0

        for t in np.arange(1, len(delta_gmt), 1):

            sl_rate = self.winkelmann_factor * \
                self.smb_coeff * delta_gmt[t - 1]
            sl_contrib[t] = sl_contrib[t - 1] + self.dtime * sl_rate

        # to m slr
        return -sl_contrib * area_antarctica / area_ocean


class solid_ice_discharge_gis(contribution):

    """ Solid ice discharge from the Greenland ice sheet. There is no equilibrium
        estimate available for solid ice discharge. We therefore use a response
        function approach, see equ. 4 in Materials and Methods. """

    def __init__(self, parameters, temp_anomaly_year=None):

        self.alpha = parameters[0]
        self.prefactor = parameters[1]
        self.dtime = dtime
        self.temp_anomaly_year = temp_anomaly_year

    # @profile
    def calc_contribution(self, temperature, proj_period, prefactor=None):
        """ No pursuit curve, but response function approach. """

        # for projection, tau needs not to be explicitely provided
        # as argument. in calibration it has to.
        if prefactor is None: prefactor = self.prefactor

        oceantemp = temperature

        if self.temp_anomaly_year is not None:
            # do not let temperature drive ice loss before first year of SL
            # observation
            oceantemp -= oceantemp[self.temp_anomaly_year]
            oceantemp[oceantemp.time < self.temp_anomaly_year] = 0.

        otemp = oceantemp[proj_period].values
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


class solid_ice_discharge_ais(contribution):

    """ Antarctic ice sheet solid ice discharge,
        equilibrium sea level response and transient response.
        See M. Mengel et al., PNAS (2016), Materials and Methods
    """


    def equilibrium_sl(self, delta_temp):
        """ see equ. 6 in Materials and Methods. """
        return self.alpha * delta_temp



class antarctica_dp16(contribution):

    """ An emulator of the Antarctic sea level response as in
        Deconto & Pollard, Nature (2016). This emulator is the same
        as in Nauels et al. 2017, ERL, submitted.
    """

    def __init__(self, parameters, temp_anomaly_year=None):

        self.parameters = parameters
        self.temp_anomaly_year = temp_anomaly_year

        # maximum ice volume that can be lost
        # see github.com/matthiasmengel/fast_ant_sid for why
        # we use this hardcoded value.
        max_volume_to_lose = 17560. # in m
        self.initial_icesheet_vol = max_volume_to_lose


    def square(arg):

        """ quadratic temperature sensitvity """

        return np.sign(arg)*np.square(arg)


    def calc_solid_ice_discharge(self, forcing_temperature,
        initial_icesheet_vol, temp_sensitivity=square):

        """ Solid ice discharge as used in Nauels et al. 2017, ERL, submitted.
            Assuming yearly time step, or dtime=1.
        """

        sid_sens, fast_rate, temp0, temp_thresh = self.parameters

        def slow_discharge(volume, temperature, temp0, sid_sens):
            # negative rates mean ice volume loss
            return sid_sens*volume*temp_sensitivity(temperature-temp0)

        ## time spans forcing period
        time = np.arange(0,len(forcing_temperature),1)

        icesheet_vol = np.zeros_like(forcing_temperature)
        icesheet_vol[0] = initial_icesheet_vol
        slr_from_sid = np.zeros_like(forcing_temperature)

        # negative rates lead to sea level rise
        fast_sid = fast_rate*np.array(forcing_temperature > temp_thresh,
                                          dtype = np.float)

        for t in time[0:-1]:
            sid_slow = slow_discharge(icesheet_vol[t], forcing_temperature[t],
                                       temp0, sid_sens)

            # yearly discharge rate cannot be larger than remaining volume
            discharge = np.minimum(icesheet_vol[t], sid_slow+fast_sid[t])
            # positive sid rates mean ice volume loss
            icesheet_vol[t+1] = icesheet_vol[t] - discharge
            slr_from_sid[t+1] = initial_icesheet_vol - icesheet_vol[t+1]

        return slr_from_sid

    def calc_contribution(self, temperature, proj_period):

        # make the temperature relative to year 1850, as we calibrated
        # the fast_ant_sid contribution to gmt relative to 1850.
        temperature -= temperature[1850]
        # from dimarray to numpy array
        temperature = temperature[proj_period].values

        # return in meter
        return self.calc_solid_ice_discharge(temperature,
                   self.initial_icesheet_vol)/1.e3



contributor_functions = {"thermexp":thermal_expansion, "gic":glaciers_and_icecaps,
                         "gis_smb":surfacemassbalance_gis, "gis_sid":solid_ice_discharge_gis,
                         "ant_smb":surfacemassbalance_ais, "ant_sid":solid_ice_discharge_ais,
                         "ant_dp16":antarctica_dp16}