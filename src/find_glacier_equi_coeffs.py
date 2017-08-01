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


""" find coefficients for analytic glacier equilibrium functions from
    Radic & Hock (2010) and Marzeion et al, Cryosphere 2012 data.
    See 'Estimates of equilibrium contribution for glaciers and ice caps'
    in the supplementary material in Mengel et al., PNAS 2016.
    See also supplementary Fig. 4.
    This script just needs to be run if glacier equilibrium estimates change.
    Coefficients are stored in data/glacier_equi/glacier_equi_coefficients.csv
"""


import os
import scipy.io
import numpy as np
from scipy import optimize
import get_calibration_data as gcd
reload(gcd)
import get_gmt_data as ggd
reload(ggd)

basedir = os.path.dirname(
    os.path.dirname(os.path.realpath(__file__)))
inputdatadir = os.path.join(basedir, "data/input/")

######## Glacier equilibrium estimates ########

# Marzeion et al. (2012) equilibrium estimates
gic_equ_marzeion12 = scipy.io.loadmat(
    inputdatadir + "glaciers_equi/marzeion12_new.mat")
equi_exp_labels = [exp[0][0] for exp in gic_equ_marzeion12["experiments"]]
# this is relative to the 1961-1990 mean temperature
gic_temp_marzeion12 = np.array(
    np.squeeze(
        gic_equ_marzeion12["t_equi"]),
    dtype="float")
# make it relative to preindustrial
gic_temp_marzeion12 += ggd.preind_to_1961_1990
# remove the natural contribution from Marzeion equilibrium data.
# natural contribution is a temperature-independent value that can be approximated as
# the magnitude of natural contribution between 1851 and 2010, see green line in Fig. 1d, in
# http://www.sciencemag.org/content/345/6199/919.abstract
# we assume glaciers to be in equilibrium before 1851, as glacier volumes did not change much,
# see the leclerq 11 timeseries.
nat_model_mean = gcd.marzeion_gic_nat_up.mean(axis="model") / 1000.  # in m
natural_gic_offset = nat_model_mean[2012] - nat_model_mean[1851]
gic_equi_marzeion12 = gic_equ_marzeion12[
    "sle_equi"] / 1000. - natural_gic_offset

## Radic & Hock
gic_equi_radic = np.genfromtxt(
    inputdatadir +
    "glaciers_equi/Cmip3_sl_equ_gic.txt")
# see the Cmip3_sl_equ_gic_info.txt, this is (presumably) relative to
# preindustrial
gic_temp_radic = np.arange(0., 7., 1.)


def gic_equi_func(equi_temp, a, b):
    """ assume exponential functional form. This is justified, because
    1) the function saturates for high temperatures.
    2) the function passes zero for zero temperature, as wanted for
       anthropogenic glacier contribution.
    """

    return a * (1 - np.exp(-1 * b * equi_temp))


def get_equi_coefficients(temp, equi_estimate):
    """ fit exponential function and retrieve coefficients. """

    temp_nonan = temp[~np.isnan(equi_estimate)]
    equi_estimate_nonan = equi_estimate[~np.isnan(equi_estimate)]
    popt, pcov = optimize.curve_fit(gic_equi_func, temp_nonan,
                                    equi_estimate_nonan, p0=[1., 1.])
    return popt


gic_equi_coeffs = []
for gic_equi in gic_equi_marzeion12:
    gic_equi_coeffs.append(
        get_equi_coefficients(
            gic_temp_marzeion12,
            gic_equi))
for gic_equi in gic_equi_radic:
    gic_equi_coeffs.append(get_equi_coefficients(gic_temp_radic, gic_equi))

gic_equi_coeffs = np.array(gic_equi_coeffs)
np.savetxt(os.path.join(basedir, "data/glacier_equi/glacier_equi_coefficients.csv"),
           gic_equi_coeffs, fmt='%.8f', header="a        b")
