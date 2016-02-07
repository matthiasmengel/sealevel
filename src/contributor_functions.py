import os
import numpy as np
import dimarray as da
import sealevel as sl; reload(sl)

project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
inputdatadir = os.path.join(project_dir,"inputdata/")

dtime = 1
area_ocean      = 3.61132e14 # % m^2
area_antarctica = 0.14e14# m^2

# Greenland subsurface ocean warming from yin et al 2011, for different basins
owarming_yin_2100 = np.array([1.61,1.16,1.25])
gmt_warming = 2.72
gis_gmt_subocean_ratio  = owarming_yin_2100.mean()/gmt_warming
# gis_gmt_subocean_ratio  = 0.2
## Antarctica
ais_gmt_subocean_ratio  = 0.2 # roughly from levermann_winkelmann14

gis_polar_amp_factor = 1.5


################################################

class thermal_expansion(object):


  def __init__(self,alpha,temp_anomaly_year=None):

    """ temp anomaly is unused here """
    self.dtime = dtime
    self.alpha = alpha

  def te_equilibrium_sl(self,delta_temp):

    """ alpha is the equilibrium sl sensitivity in m/K """

    return self.alpha * delta_temp


  def calc_contribution(self,delta_gmt,tau):

    # oceantemp = ocean_temperature.copy()
    # oceantemp[oceantemp.time<1992] = 0.
    delta_gmt  = delta_gmt.values
    sl_contrib = np.zeros_like(delta_gmt,dtype="float")
    sl_rate    = 0.0

    for t in np.arange(1,len(delta_gmt),1):

      sl_rate = (self.te_equilibrium_sl(delta_gmt[t-1])-sl_contrib[t-1])/tau
      # sl_rate = np.maximum(sl_rate,0.)
      sl_contrib[t] = sl_contrib[t-1]+self.dtime*sl_rate

    return sl_contrib


################################################

class glaciers_and_icecaps(object):


  def __init__(self,modelno,temp_anomaly_year=None):

    """ temp anomaly is unused here """

    self.dtime    = dtime
    # self.wl       = warming_levels
    # self.eq_sl_wl = equilibrium_sl_for_warming_levels
    self.modelno  = modelno
    self.temp_anomaly_year = temp_anomaly_year


  def gic_equilibrium_sl(self,delta_temp):

    return  sl.gic_equi_functions[self.modelno](delta_temp)


  def calc_contribution(self,delta_gmt,tau):


    # delta_gmt  -= delta_gmt[1880:1890].mean()
    # print delta_gmt.values

    if self.temp_anomaly_year != None:
      ## do not count temperature as a driver before first year of SL observation
      delta_gmt -= delta_gmt[self.temp_anomaly_year]
      delta_gmt[delta_gmt.time<self.temp_anomaly_year] = 0.

    delta_gmt  = delta_gmt.values
    sl_contrib = np.zeros_like(delta_gmt,dtype="float")
    sl_rate    = 0.0

    for t in np.arange(1,len(delta_gmt),1):
      sl_rate = (self.gic_equilibrium_sl(delta_gmt[t]) - sl_contrib[t-1])/tau
      # sl_rate = np.maximum(sl_rate,0.)
      sl_contrib[t] = sl_contrib[t-1]+self.dtime*sl_rate

    return sl_contrib


################################################

class surfacemassbalance_gis(object):

  def __init__(self,smb_coeff,temp_anomaly_year=None):

    self.dtime     = dtime
    self.smb_coeff = smb_coeff
    self.temp_anomaly_year = temp_anomaly_year
    # print self.temp_anomaly_year
    # self.calibrate = calibrate


  def equilibrium_sl(self,delta_temp):

    """ from fig.2 levermann commitment paper """

    return self.smb_coeff*np.sign(delta_temp)*delta_temp**2

  def calc_contribution(self,delta_gmt,tau):

    if self.temp_anomaly_year != None:
      ## do not count temperature as a driver before first year of SL observation
      delta_gmt -= delta_gmt[self.temp_anomaly_year]
      delta_gmt[delta_gmt.time<self.temp_anomaly_year] = 0.

    delta_gmt  = delta_gmt.values
    # print delta_gmt
    sl_contrib = np.zeros_like(delta_gmt,dtype="float")
    sl_rate    = 0.0

    for t in np.arange(1,len(delta_gmt),1):

      sl_rate       = (self.equilibrium_sl(delta_gmt[t])-sl_contrib[t-1])/tau
      sl_contrib[t] = sl_contrib[t-1]+self.dtime*sl_rate

    return sl_contrib


################################################

class surfacemassbalance_ais(object):

  """ no dogs curve style, just scaling slr with clausius clapeyron like factor. """

  def __init__(self,smb_coeff,temp_anomaly_year=None):

    self.dtime     = dtime
    self.smb_coeff = smb_coeff
    # 25 per cent will be lost due to increased discharge.
    self.winkelmann_factor = 0.75

  def calc_contribution(self,delta_gmt, dummy):

    """ dummy without influence """

    delta_gmt  = delta_gmt.values
    sl_contrib = np.zeros_like(delta_gmt,dtype="float")
    sl_rate    = 0.0

    for t in np.arange(1,len(delta_gmt),1):

      sl_rate       = self.winkelmann_factor * self.smb_coeff*delta_gmt[t-1]
      sl_contrib[t] = sl_contrib[t-1]+self.dtime*sl_rate
      # print sl_rate,delta_gmt[t-1]

    # to m slr
    return -sl_contrib*area_antarctica/area_ocean


################################################

class solid_ice_discharge_gis(object):

  def __init__(self,alpha,temp_anomaly_year=None):

    self.alpha = alpha
    self.dtime = dtime
    # tau in years; does only play a role in saturation, so not until 2100
    self.tau   = 40.
    # self.calibrate = calibrate
    self.temp_anomaly_year = temp_anomaly_year

  def calc_contribution(self,temperature,prefactor):

    # if self.calibrate:
    oceantemp = temperature
    # else:
    #   delta_gmt = temperature
    #   oceantemp = gis_gmt_subocean_ratio*delta_gmt
    # print self.temp_anomaly_year

    if self.temp_anomaly_year != None:
      ## do not count temperature as a driver before first year of SL observation
      oceantemp -= oceantemp[self.temp_anomaly_year]
      oceantemp[oceantemp.time<self.temp_anomaly_year] = 0.
      # print oceantemp.values

    # oceantemp  = gis_gmt_subocean_ratio*delta_gmt

    otemp = oceantemp.values
    # print otemp
    discharge = np.zeros_like(otemp,dtype="float")

    for t in np.arange(1,len(otemp),1):
      discharge_rate = 0
      # integrate over t to yield slr, see winkelmann response function paper p. 2581
      tp = np.arange(t)
      discharge_rate = prefactor*np.trapz(otemp[tp]*((t-tp))**self.alpha,dx=self.dtime)
      # print t, (t-tp)/self.tau, discharge_rate
      # discharge_rate = np.minimum(discharge_rate,(0.13 - discharge[t-1])/self.tau)
      discharge_rate = np.maximum(discharge_rate,0.0)
      discharge[t]   = discharge[t-1]+self.dtime*discharge_rate

    return discharge


################################################

class solid_ice_discharge_ais(object):

  def __init__(self,alpha,temp_anomaly_year=None):

    self.alpha     = alpha
    self.dtime     = dtime
    self.temp_anomaly_year = temp_anomaly_year

  def equilibrium_sl(self,delta_temp):

    """ from fig.2 levermann commitment paper """

    return self.alpha*delta_temp

  def calc_contribution(self,temperature,tau):

    if self.temp_anomaly_year != None:
      ## do not count temperature as a driver before first year of SL observation
      temperature -= temperature[self.temp_anomaly_year]
      temperature[temperature.time<self.temp_anomaly_year] = 0.


    temperature = temperature.values
    discharge = np.zeros_like(temperature,dtype="float")
    discharge_rate = 0

    for t in np.arange(1,len(temperature),1):

      discharge_rate = (self.equilibrium_sl(temperature[t]) - discharge[t-1])/tau
      discharge_rate = np.maximum(discharge_rate,0.0)
      discharge[t]   = discharge[t-1]+self.dtime*discharge_rate

    return discharge


################################################

class total_contribution(object):

  """ not used. """

  def __init__(self,te_alpha,gic_modelno,gis_smbcoeff,):

    self.dtime = dtime
    self.alpha = alpha

