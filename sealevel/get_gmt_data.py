import os
import numpy as np
import dimarray as da

project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
inputdatadir = os.path.join(project_dir, "data/input/")


######## HADCRUT4 global mean temperature observations ########
try:
    hadcrut = np.loadtxt(inputdatadir + "hadcrut4/HadCRUT.4.6.0.0.annual_ns_avg.txt")
except IOError:
    raise IOError(
        "HADcrut v4.0 temperature data missing, "
        + "please run src/download_input_data.py first"
    )


hadcrut_temp = da.DimArray(hadcrut[:, 1], axes=hadcrut[:, 0], dims="time")
# define "preindustrial" as 1850:1860
preind_to_1951_1980 = hadcrut_temp[1951:1980].mean() - hadcrut_temp[1850:1860].mean()
preind_to_1961_1990 = hadcrut_temp[1961:1990].mean() - hadcrut_temp[1850:1860].mean()
hadcrut_temp -= hadcrut_temp[1850:1860].mean()


######## GISS global mean temperature observations ########
try:
    giss_landocean = np.loadtxt(
        inputdatadir + "gisstemp/giss_landocean_2013.txt", usecols=(0, 13)
    )
except IOError:
    raise IOError(
        "Giss global mean temperature data missing, "
        + "please run src/download_input_data.py first."
    )

giss_temp = giss_landocean[:, 1] * 0.01  # column 14 is the annual global value
giss_years = giss_landocean[:, 0]
giss_temp = da.DimArray(giss_temp, axes=giss_years, dims="time")  # [1900:]
## make giss global mean temperature relative to preindustrial, which is here defined
## as the 1850:1860 hadcrut mean. The global mean temperature anomaly defined this way
## is used throughout the code for calibration and projections.
giss_temp = giss_temp - giss_temp[1951:1980].mean()
giss_temp = giss_temp + preind_to_1951_1980
