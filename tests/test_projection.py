import os
import numpy as np
import pandas as pd
import settings
import sealevel.projection
import sealevel.get_ipcc_data
import sealevel.contributor_functions as cf


temp_anomaly_years = pd.read_csv(
    os.path.join(settings.calibfolder, "temp_anomaly_years.csv"), index_col=[0, 1]
)
temp_anomaly_years = temp_anomaly_years.where(pd.notnull(temp_anomaly_years), None)

# test only one scenario for now.
scen = "rcp85"
gmt = sealevel.get_ipcc_data.tas_data[scen]

# for testing let us use a small number
nrealizations = 10
realizations = np.arange(nrealizations)
test_path = os.path.dirname(__file__)
# test data without anomalization, so starting with slr=0 in first entry.
proj_testdata_n10 = pd.read_csv(
    os.path.join(test_path, "data/projection/" + scen + "n10.csv"), index_col=0
)


def get_percentiles(data):

    """ assumes year 2100 is last element in data. Return in mm. """
    return np.percentile(data[-1, :], [5, 50, 95]) * 1.0e3


def test_thermexp_projection():

    contrib_name = "thermexp"
    calibdata = pd.read_csv(
        os.path.join(settings.calibfolder, contrib_name + ".csv"), index_col=[0]
    )

    temp_anomaly_year = temp_anomaly_years.loc[contrib_name]
    sl_contributor = cf.contributor_functions[contrib_name]

    proj = np.zeros([len(settings.proj_period), nrealizations])

    for n in realizations:
        slr, gmt_n, obs_choice, params = sealevel.projection.project(
            gmt,
            settings.proj_period,
            calibdata,
            temp_anomaly_year,
            sl_contributor,
            n,
            contrib_name,
        )
        proj[:, n] = slr

    np.testing.assert_allclose(
        get_percentiles(proj), proj_testdata_n10.loc[contrib_name]
    )
