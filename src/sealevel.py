import os
import numpy as np
import pandas as pd

project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
inputdatadir = os.path.join(project_dir,"data/")

######## paramerters that need to be known for calibrations and projections ########

# add temperature offset to be used with the box & colgan 2013 data to fit past observations
# similar offset is used by e.g. Rahmstorf 2007, Science.
gis_colgan_temperature_offset = 0.5


######## equilibrium sea level functions ########

def gic_equi_func(equi_temp,a,b):

    """ assume exponential functional form. This is justified, because
    1) the function saturates for high temperatures.
    2) the function passes zero for zero temperature, as wanted for
       anthropogenic glacier contribution.
    """

    return a*(1-np.exp(-1*b*equi_temp))


def func_creator(a,b):

    def func(equi_temp):
        """ same as gic_equi_func with explicit values for a and b as set by curve_fit. """
        return a*(1-np.exp(-1*b*equi_temp))

    return func


gic_equi_coeffs = pd.read_csv(inputdatadir+"glacier_equi/glacier_equi_coefficients.csv")
gic_equi_functions = []
for coeffs in gic_equi_coeffs.values:
    gic_equi_functions.append(func_creator(*coeffs))


######## project function for Monte-Carlo sampling


def project(gmt, proj_period, calibdata, sample_number):

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
    sample_number : number for seed to be created to make sampling reproducible

    Returns
    -------
    contrib : timeseries of sea level contribution with length proj_period
    """

    np.random.seed(sample_number)

    try:
        gmt_ensemble_size = gmt.shape[1]
        gmt_choice = np.random.randint(gmt_ensemble_size)
        #print gmt_choice
        driving_temperature = gmt[proj_period,gmt_choice]
    except IndexError:
        ## this is the case if single gmt is supplied
        driving_temperature = gmt[proj_period]

    # print contrib_name, temp_anomaly_year
    # use one of the observational dataset
    obs_choice = np.random.choice(calibdata["params"].keys())
    params = calibdata["params"][obs_choice]
    temp_anomaly_year = params.temp_anomaly_year

    if obs_choice == "box_colgan13":
        driving_temperature += gis_colgan_temperature_offset

    # choose a random parameter set
    tuple_choice = np.random.randint(len(params.commitment_parameter))
    independent_param = params.commitment_parameter[tuple_choice]
    dependent_param   = params.calibrated_parameter[tuple_choice]
    # print "tuple",independent_param,dependent_param

    contributor = params.sl_contributor(independent_param, temp_anomaly_year)
    contrib = contributor.calc_contribution(driving_temperature,dependent_param)
    # print contrib
    return contrib
