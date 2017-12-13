""" Use this file to download publicly available data for calibration. """

import subprocess

datadir = "../data/input/"


def download(info, source, targetdir, target=None):

    print "##", info, "##"

    targetfolder = datadir + targetdir
    cmd = "mkdir -p " + targetfolder
    # print cmd
    subprocess.check_call(cmd, shell=True)

    if target is None:
        subprocess.check_call(
            "wget -c " +
            source +
            " -P " +
            targetfolder +
            "/",
            shell=True)
    else:
        subprocess.check_call(
            "wget -c " +
            source +
            " -O " +
            targetfolder +
            "/" +
            target,
            shell=True)


info = "Hadcrut4 global mean temperature timeseries"
source = "https://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/time_series/HadCRUT.4.6.0.0.annual_ns_avg.txt"
download(info, source, "hadcrut4")

info = "GISSTEMP global mean temperature series"
source = "http://data.giss.nasa.gov/gistemp/tabledata_v3/GLB.Ts+dSST.txt"
download(info, source, "gisstemp")

# remove all empty lines, all header lines starting with "Year" and all lines after 2013
# though years 2014 and 2015 are available now, we remove them to ensure consistency with
# data as used in the analysis of M. Mengel et al., PNAS (2016)
sedcmd = "sed '/^Year/ d' GLB.Ts+dSST.txt | sed '/^$/d' | sed '1,5d' | sed '135,$d'"
subprocess.check_call(
    "cd " +
    datadir +
    "gisstemp && " +
    sedcmd +
    " > giss_landocean_2013.txt",
    shell=True)


info = "IPCC AR5 global mean temperature data from CMIP5 models"
source = "http://www.ipcc.ch/report/ar5/wg1/docs/ar5_wg1_annexI_all.zip"
download(info, source, "ipcc_ar5")
cmd = "cd " + datadir + \
    "ipcc_ar5 && unzip ar5_wg1_annexI_all.zip **WGIAR5_FD_AnnexI_series_tas_modelmean_*_world_annual.txt"
subprocess.check_call(cmd, shell=True)
cmd = "cd " + datadir + \
    "ipcc_ar5 && mv -v WGIAR5_AnnexI_timeseries/tas/*/WGIAR5_FD_AnnexI_series_tas_modelmean_*_world_annual.txt ./"
subprocess.check_call(cmd, shell=True)
