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
import dimarray as da
import scipy
import scipy.interpolate
import scipy.io
import collections
from scipy import optimize

project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
inputdatadir = os.path.join(project_dir, "data/input/")

contrib_ids = [
    "thermexp",
    "gic_anth",
    "gis_smb",
    "gis_sid",
    "ant_smb",
    "ant_sid"]

######## Church and White 2011 sea level data ########

# http://onlinelibrary.wiley.com/doi/10.1029/2011GL048794/abstract
sealevelbudget = np.loadtxt(
    inputdatadir +
    "Observations/SeaLevelBudgetChurch.txt")  # ,delimiter=" ")
budgetlbl = [
    "satellite",
    "tide gauge",
    "dummy",
    "thermexp",
    "gic",
    "gis",
    "ant",
    "",
    "terrestrial_storage"]
years_church = np.array(sealevelbudget[:, 0] - 0.5, dtype="int")

church_observed = {}
for j, lbl in enumerate(budgetlbl):
    church_observed[lbl] = da.DimArray(sealevelbudget[:, j] * 1e-3,
                                       axes=years_church, dims="time")  # in m


######## Church et al. 2011; 1880 - today historical sea-level curve ########

# Survey in Geophysics 2011,
# link.springer.com/article/10.1007%2Fs10712-011-9119-1
church_past_data = np.loadtxt(inputdatadir +
                              "church_white11/CSIRO_Recons_gmsl_yr_2011.txt")
years_church_past = church_past_data[:, 0] - 0.5
church_past_sl = da.DimArray(church_past_data[:, 1] * 0.001,
                             axes=years_church_past, dims="time")
church_past_sl = church_past_sl - church_past_sl[1950:1980].mean()
church_past_sl_std = da.DimArray(church_past_data[:, 2] * 0.001,
                                 axes=years_church_past, dims="time")


######## Hay et al, Nature January 2015; doi:10.1038/nature14093 ########

haydata = np.loadtxt(inputdatadir + "hay2015_sealevel/fig02_blueline.csv")
fn = scipy.interpolate.interp1d(haydata[:, 0], haydata[:, 1] * 1e-3)  # in mm
hay15_past_sl = da.DimArray(
    fn(np.arange(1901, 2011, 1)), axes=np.arange(1901, 2011, 1), dims="time")

haydata = np.loadtxt(
    inputdatadir +
    "hay2015_sealevel/hay15_lower.csv",
    delimiter=",")
fn = scipy.interpolate.interp1d(haydata[:, 0], haydata[:, 1] * 1e-3)  # in mm
hay15_lower = da.DimArray(fn(np.arange(1901, 2010, 1)),
                          axes=np.arange(1901, 2010, 1), dims="time")

haydata = np.loadtxt(
    inputdatadir +
    "hay2015_sealevel/hay15_upper.csv",
    delimiter=",")
fn = scipy.interpolate.interp1d(haydata[:, 0], haydata[:, 1] * 1e-3)  # in mm
hay15_upper = da.DimArray(fn(np.arange(1901, 2010, 1)),
                          axes=np.arange(1901, 2010, 1), dims="time")


######## Thermal expansion observations ########

# Catia Domingues 0-700m
# http://www.nature.com/nature/journal/v453/n7198/abs/nature07080.html
# downloaded from
# http://www.cmar.csiro.au/sealevel/thermal_expansion_ocean_heat_timeseries.html
domingues08 = np.loadtxt(inputdatadir + "thermosteric_katjad/catia.txt")
thermo_obs_domingues08 = da.DimArray(
    domingues08[
        :, 1], axes=domingues08[
            :, 0] - 0.5, dims="time") / 1000.

# Ishii & Kimoto 0-700m
# http://link.springer.com/article/10.1007/s10872-009-0027-7
# downloaded from http://amaterasu.ees.hokudai.ac.jp/~ism/pub/ProjD/v6.14
ishii09 = np.loadtxt(
    inputdatadir +
    "thermosteric_ishii2009/global_mean_thermosteric_0-700m.txt",
    skiprows=2)
thermo_obs_ishii09 = da.DimArray(
    ishii09[
        :, 1], axes=ishii09[
            :, 0], dims="time") / 1000.

# Levitus 2012
# http://onlinelibrary.wiley.com/doi/10.1029/2012GL051106/abstract
# downloaded from
# http://www.nodc.noaa.gov/OC5/3M_HEAT_CONTENT/basin_tsl_data.html
levitus12_700m = np.loadtxt(
    inputdatadir +
    "thermosteric_levitus2012/pent_a-mm-w0-700m.dat",
    skiprows=1)
levitus12_2000m = np.loadtxt(
    inputdatadir +
    "thermosteric_levitus2012/pent_a-mm-w0-2000m.dat",
    skiprows=1)
thermo_obs_levit12_700 = da.DimArray(
    levitus12_700m[
        :, 1], axes=levitus12_700m[
            :, 0] - 0.5, dims="time") / 1000.
thermo_obs_levit12_2000 = da.DimArray(
    levitus12_2000m[
        :, 1], axes=levitus12_2000m[
            :, 0] - 0.5, dims="time") / 1000.

# Purkey & Johnson 2010
# http://journals.ametsoc.org/doi/abs/10.1175/2010JCLI3682.1
# estimates for deep ocean below 2000m
# overlong time for trend, the length will be shortened by upper ocean
# observation length.
purkey10_below2000m = da.DimArray(
    np.arange(64) *
    0.113e-3,
    axes=np.arange(
        1951,
        2015,
        1),
    dims="time")
purkey10_below3000m = da.DimArray(
    np.arange(64) *
    0.097e-3,
    axes=np.arange(
        1951,
        2015,
        1),
    dims="time")

# Llovel et al 2014
# www.nature.com/articles/nclimate2387
llovel14_below2000m = da.DimArray(
    np.arange(
        len(years_church)) * -0.13e-3,
    axes=years_church,
    dims="time")

# trend 700-3000m, see paper
# time axis is artifically expanded here. will be restricted by upper
# ocean observations.
church_observed["te_700_3000"] = da.DimArray(
    np.arange(64) *
    0.07e-3,
    axes=np.arange(
        1951,
        2015,
        1),
    dims="time")
# trend below 3000
church_observed["te_below_3000"] = da.DimArray(
    np.arange(64) * 0.09e-3, axes=np.arange(1951, 2015, 1), dims="time")

church_observed_700m = church_observed[
    "thermexp"] - church_observed["te_700_3000"] - church_observed["te_below_3000"]
church_observed_700m_2000m = church_observed[
    "te_700_3000"] + church_observed["te_below_3000"] - purkey10_below2000m


# thermal expansion from cmip5 models, taken from fig 13.4 ch13 ipcc ar5 wg1
ipcc_te_data = np.loadtxt(
    inputdatadir +
    "ipcc_ar5/thermal_contribution_ch13_fig13p4.csv")
f = scipy.interpolate.interp1d(
    ipcc_te_data[
        :, 0], ipcc_te_data[
            :, 1] * 1e-3)  # in mm
ipcc_thermal_expansion = da.DimArray(
    f(np.arange(1902, 2010, 1)), axes=np.arange(1902, 2010, 1), dims="time")
ipcc_thermal_expansion -= ipcc_thermal_expansion[1986:2005].mean()

# thermal expansion range 1901-1990, table 13.1 ch 13 IPCC
# in mm
ipcc_te_lower = da.DimArray(np.linspace(
    0., 0.06, 1991 - 1901) * (1991 - 1901), axes=np.arange(1901, 1991, 1), dims="time")
ipcc_te_med = da.DimArray(np.linspace(0., 0.37, 1991 - 1901)
                          * (1991 - 1901), axes=np.arange(1901, 1991, 1), dims="time")
ipcc_te_upper = da.DimArray(np.linspace(
    0., 0.67, 1991 - 1901) * (1991 - 1901), axes=np.arange(1901, 1991, 1), dims="time")


######## Glacier observations ########

# all data in m SLE (sea level rise is positive), starting with 0 at beginning of
# respective observation
gls = scipy.io.loadmat(
    inputdatadir +
    "glaciers_marzeion/glacier_mass_change_for_matthias.mat")

# Paul Leclercq
# http://link.springer.com/article/10.1007/s10712-011-9121-7
glacier_obs_leclercq11 = da.DimArray(-1. * gls['dV_new_paul_smoothed'].squeeze() / 1.e3,
                                     axes=gls['time_paul_new'].squeeze(), dims="time")
glacier_obs_leclercq11 -= glacier_obs_leclercq11.ix[0]

# Graham Cogley
# http://www.ingentaconnect.com/content/igsoc/agl/2009/00000050/00000050/art00014
# turn pentadal rates [m/y] into yearly sea level [m]
glacier_obs_cogley09 = da.DimArray(-1. * gls['global_excl_aa_dv'].repeat(5).cumsum() / 3.61e2 / 1.e3,
                                   axes=np.arange(1961, 2011), dims="time")
glacier_obs_cogley09 -= glacier_obs_cogley09.ix[0]

# Ben Marzeion
# http://www.the-cryosphere.net/6/1295/2012/tc-6-1295-2012.pdf
glacier_obs_marzeion12 = da.DimArray(-1. * gls['glacier_mass_change'].squeeze() / 1.e3,
                                     axes=np.array(gls['time'].squeeze(), dtype=np.int), dims="time")
glacier_obs_marzeion12 -= glacier_obs_marzeion12.ix[0]


######## Anthropogenic fraction of glacier contribution ########

# Marzeion et al 2014, antropogenic fraction of glaciers.
# http://www.sciencemag.org/content/345/6199/919.abstract
# updated data, in mm

natfile = inputdatadir + "glaciers_marzeion14/global_cumulative_mass_loss_nat_rgi_v4.txt"
marzeion_gic_nat_up = np.loadtxt(natfile, skiprows=2)
models = np.genfromtxt(natfile, skip_header=1, dtype=None)[0, 1:]
marzeion_gic_nat_up = da.DimArray(marzeion_gic_nat_up[:, 1:],
                                  axes=[marzeion_gic_nat_up[:, 0], models], dims=["time", "model"])

fullfile = inputdatadir + \
    "glaciers_marzeion14/global_cumulative_mass_loss_full_rgi_v4.txt"
marzeion_gic_full_up = np.loadtxt(fullfile, skiprows=2)
models = np.genfromtxt(fullfile, skip_header=1, dtype=None)[0, 1:]
marzeion_gic_full_up = da.DimArray(marzeion_gic_full_up[:, 1:],
                                   axes=[marzeion_gic_full_up[:, 0], models], dims=["time", "model"])

natfile = inputdatadir + "glaciers_marzeion14/global_mean_specific_mb_nat_rgi_v4.txt"
marzeion_gic_nat_mb = np.loadtxt(natfile, skiprows=2)
models = np.genfromtxt(natfile, skip_header=1, dtype=None)[0, 1:]
marzeion_gic_nat_mb = da.DimArray(marzeion_gic_nat_mb[:, 1:],
                                  axes=[marzeion_gic_nat_mb[:, 0], models], dims=["time", "model"])

fullfile = inputdatadir + "glaciers_marzeion14/global_mean_specific_mb_full_rgi_v4.txt"
marzeion_gic_full_mb = np.loadtxt(fullfile, skiprows=2)
models = np.genfromtxt(fullfile, skip_header=1, dtype=None)[0, 1:]
marzeion_gic_full_mb = da.DimArray(marzeion_gic_full_mb[:, 1:],
                                   axes=[marzeion_gic_full_mb[:, 0], models], dims=["time", "model"])

anthropogenic_frac = (marzeion_gic_full_mb.mean(axis=1)
                      - marzeion_gic_nat_mb.mean(axis=1)) / marzeion_gic_full_mb.mean(axis=1)
anthropogenic_frac_perc = ((marzeion_gic_full_mb.T.to_pandas().quantile([0.05, 0.5, 0.95])
                            - marzeion_gic_nat_mb.T.to_pandas().quantile([0.05, 0.5, 0.95]))
                           / marzeion_gic_full_mb.T.to_pandas().quantile([0.05, 0.5, 0.95]))
anthropogenic_frac_perc = da.from_pandas(
    anthropogenic_frac_perc, dims=[
        "percentile", "time"])
# TODO: reference motivate linear increase of anth fraction
#anthropogenic_frac_lin = da.DimArray(np.linspace(0,1.,2056-1870),axes=np.arange(1870,2056,1),dims="time")

# linear fit to fraction since 1900
tlin = np.arange(1870, 2056, 1)
anthropogenic_frac_lin_1900 = da.DimArray(
    6.78606933e-03 *
    tlin -
    1.28443788e+01,
    axes=tlin,
    dims="time")


######## Greenland surface mass balance ########

def convert_raw_mb_to_sl(rawdata):

    fn = scipy.interpolate.interp1d(rawdata[:, 0], rawdata[:, 1])
    # in Gt/y
    mb = da.DimArray(fn(np.arange(1871, 2009, 1)),
                     axes=np.arange(1871, 2009, 1), dims="time")
    # assume "equilibrium" in years 1961-1990 or 1971-1988, see discussion in Hanna et al 2011
    # and in Rignot et al. 2008
    # doi:10.1029/2011JD016387
    # assume "equilibrium" before 1900, see the zero before 1900 in Fig. 7,
    # Box and Colgan 2013
    mb -= mb[1870:1900].mean()
    # convert to m SLE
    return mb.cumsum() / 360. / 1.e3

# Box and Colgan, Journal of Climate 2013, reconstruction of total
# Greenland mass balance
boxdata = np.loadtxt(
    inputdatadir +
    "greenland_box/box_colgan2013_fig6.csv",
    delimiter=",")
boxdata[:, 1] *= -1
# box_gis_tmb = convert_raw_mb_to_sl(boxdata)

# Box SMB; from Box 2013, Journal of Climate 2013
boxdata = np.loadtxt(
    inputdatadir +
    "greenland_box/box13_fig8_box.csv",
    delimiter=",")
fn = scipy.interpolate.interp1d(
    boxdata[
        :, 0], boxdata[
            :, 1] / 360. / 1.e3)  # in m/y SLE
box_smb_rate = da.DimArray(
    fn(np.arange(1871, 2009, 1)), axes=np.arange(1871, 2009, 1), dims="time")
box_smb_offset = box_smb_rate[
    1960:1990].mean() - box_smb_rate[1870:1900].mean()
# Hanna SMB; from Box 2013, Journal of Climate 2013
# exclude this from analysis, the largely increasing ice sheet mass before 1960 makes
# it very difficult to link the the GMT-driven parametrization.
# hannadata = np.loadtxt(inputdatadir+"greenland_box/box13_fig8_hanna.csv",delimiter=",")

# hanna_smb = hannadata
# add marine ice loss to raw smb data to yield total mass balance tmb (not used)
# hanna_tmb[:,1] = hanna_tmb[:,1] + marine_ice_loss(hanna_tmb[:,1])

box_gis_smb = convert_raw_mb_to_sl(boxdata)
# hanna_gis_smb = convert_raw_mb_to_sl(hannadata)

# assume 50% of church mass loss was due to smb,
# recent estimates, though for shorter time periods: andersen_senseng15: 61% smb
# csatho_schenk14 over period 2003-2009 Fig.4 and Table S5:
# Total  SMB: 763 +-53 SID: 698 +-94, so 52%
church_gis_smb = church_observed['gis'] * 0.5
church_gis_sid = church_observed['gis'] * 0.5

# Sasgen et al. 2012, http://dx.doi.org/10.1007/s10712-013-9261-z
rawdata = np.loadtxt(
    inputdatadir +
    "greenland_sasgen12/sasgen_vandenbroeke12_smb.csv",
    delimiter=",")
fn = scipy.interpolate.interp1d(
    rawdata[:, 0], -1. * rawdata[:, 1] / 360. / 1.e3)  # in m SLE
sasgen_smb = da.DimArray(fn(np.arange(1959, 2009, 1)),
                         axes=np.arange(1959, 2009, 1), dims="time")

# van Angelen et al. 2014, http://dx.doi.org/10.1016/j.epsl.2012.03.033
rawdata = np.loadtxt(inputdatadir + "greenland_vanangelen14/vanangelen_vandenbroeke14_fig11.csv",
                     delimiter=",")
fn = scipy.interpolate.interp1d(
    rawdata[:, 0], -1. * rawdata[:, 1] / 360.)  # in m SLE
angelen_smb_since_1992 = da.DimArray(
    fn(np.arange(1990, 2013, 1)), axes=np.arange(1990, 2013, 1), dims="time")
rawdata = np.loadtxt(inputdatadir + "greenland_vanangelen14/vanangelen_vandenbroeke14_fig10.csv",
                     delimiter=",")
fn = scipy.interpolate.interp1d(
    rawdata[:, 0], -1. * rawdata[:, 1] / 360. / 1000.)  # in m SLE
angelen_smbrate = da.DimArray(fn(np.arange(1960, 2013, 1)),
                              axes=np.arange(1960, 2013, 1), dims="time")
angelen_smbrate = angelen_smbrate - \
    angelen_smbrate[1960:1990].mean() + box_smb_offset
angelen14 = angelen_smbrate.cumsum()


######## Greenland solid ice discharge ########

def marine_ice_loss(smb):
    """ marine ice loss parametrization after Box and Colgan 2013
        input in Gt/y, output in Gt/y """
    return 1.088 * smb - 1068.

# Sasgen et al. 2012, http://dx.doi.org/10.1016/j.epsl.2012.03.033
rawdata = np.loadtxt(
    inputdatadir +
    "greenland_sasgen12/sasgen_vandenbroeke12_sid.csv",
    delimiter=",")
fn = scipy.interpolate.interp1d(
    rawdata[:, 0], -1. * rawdata[:, 1] / 360. / 1.e3)  # in m SLE
sasgen_sid = da.DimArray(fn(np.arange(1959, 2009, 1)),
                         axes=np.arange(1959, 2009, 1), dims="time")

# andresen straneo nature geoscience 2012, 1.5 degC (not used currently)
gis_ocean_warming_obs_1992_2008 = da.DimArray(
    np.linspace(
        0,
        1.5,
        2009 - 1992),
    axes=np.arange(
        1992,
        2009,
        1),
    dims="time")

# constructed greenland warming, 3 degree in 2 decades
greenland_warming_1992_2010 = da.DimArray(np.linspace(0, 3., 2011 - 1990),
                                          axes=np.arange(1990, 2011, 1), dims="time")


######## Antarctica solid ice discharge ########

# jacobs jenkins 2011
pine_island_ocean_warming_obs_1992_2008 = da.DimArray(np.linspace(0, 0.3, 2009 - 1970),
                                                      axes=np.arange(1970, 2009, 1), dims="time")

# schmidtko et al science 2014, trend for all regions.
trend = np.array([8.98, 22.58, -1.88, -5.06]).mean() * 1e-3  # degC/yr
# trend = np.array([8.98,22.58]).mean()*1e-3 # degC/yr
schmidtko_warming_1970_2010 = da.DimArray(np.linspace(0, trend, 2014 - 1970) * (2013 - 1970),
                                          axes=np.arange(1970, 2014, 1), dims="time")

# shepherd irvins 2012
gt_to_m = 1. / 360. * 1e-3
shepherd12_ais_lower = da.DimArray(np.linspace(0, 71. - 53., 2012 - 1992) * (2012 - 1992),
                                   axes=np.arange(1992, 2012, 1), dims="time") * gt_to_m
shepherd12_ais_upper = da.DimArray(np.linspace(0, 71. + 53., 2012 - 1992) * (2012 - 1992),
                                   axes=np.arange(1992, 2012, 1), dims="time") * gt_to_m

# mouginot rignot 14 ice flow, these are only the Amundsen glaciers.
mouginot_ase_data = np.loadtxt(
    inputdatadir +
    "mouginot_rignot14/ase_iceflow.txt")
f = scipy.interpolate.interp1d(
    mouginot_ase_data[
        :, 0], mouginot_ase_data[
            :, 1])
# integrate to yield sl. nearly in balance in 1974 (Rignot 2008)
contrib_to_slr = f(np.arange(1974, 2014))
contrib_to_slr -= contrib_to_slr[0]
mouginot_ase_sl = da.DimArray(
    np.cumsum(contrib_to_slr), axes=np.arange(
        1974, 2014), dims="time") * gt_to_m

# Harig and Simons 2015, total antarctic mass change in Gt/y
harig_data = np.loadtxt(
    inputdatadir +
    "harig_simons15/harig_simons15_fig2e_sm.csv",
    delimiter=",")
f = scipy.interpolate.interp1d(harig_data[:, 0], harig_data[:, 1])
harig_simons15 = da.DimArray(f(np.arange(2003, 2014)), axes=np.arange(
    2003, 2014), dims="time") * -1 * gt_to_m
# start at 0 for convencience.
harig_simons15 -= harig_simons15.ix[0]
