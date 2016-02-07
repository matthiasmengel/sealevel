import numpy as np
from scipy.optimize import curve_fit,leastsq
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
            self.calibrated_parameter = np.zeros_like(self.commitment_parameter)
            pass
        else:
            # restrict calibration period to available observations.
            if observation_period == None:
                print "observation period none"
                cal_period_0 = np.max([self.gmt.time[0],sl_observation.time[0]])
                cal_period_1 = np.min([self.gmt.time[-1],sl_observation.time[-1]])

            else:
                # also restrict to provided observation period
                cal_period_0 = np.max([observation_period[0],self.gmt.time[0],
                    sl_observation.time[0]])
                cal_period_1 = np.min([observation_period[-1],self.gmt.time[-1],
                    sl_observation.time[-1]])

            self.cal_period = np.arange(cal_period_0,cal_period_1+1,1)



        # if temp_anomaly_year != None:
        #     ## do not apply temperature as a driver before first year of SL observation
        #     self.gmt -= self.gmt[temp_anomaly_year]
        #     self.gmt[self.gmt.time<temp_anomaly_year] = 0.

        # print self.cal_period
        # sl_observation -= sl_observation[observation_period[0]]
        # self.sl_observation = sl_observation

        # print self.sl_observation[self.obs_period]

    def calc_residuals(self,param,sl_function):

        """ Only residuals in calibration period play a role for fit.
            Residuals start always with zero in the first year of cal_period. """

        contrib   = da.DimArray(sl_function(self.gmt, param), axes=self.gmt.time,dims="time")
        residuals = da.DimArray(np.zeros_like(contrib.values), axes=self.gmt.time,dims="time")

        residuals[self.cal_period]  = self.sl_observation[self.cal_period] - contrib[self.cal_period]
        residuals[self.cal_period] -= residuals[self.cal_period[0]]
        # alternative, but untested, residuals relative to their mean, may reduce
        # dependency on the very first value.
        # residuals[self.cal_period] -= residuals[self.cal_period].mean()

        return residuals.values


    # def find_sl_gmt_overlap(self):

    #     sl_obs_on_gmt_axis = self.gmt.copy()
    #     sl_obs_on_gmt_axis[:] = 0.
    #     sl_obs_on_gmt_axis[self.obs_period] = self.sl_obervation[self.obs_period]
    #     return sl_obs_on_gmt_axis

    def calibrate(self):

        # sl_obs_on_gmt_axis = self.find_sl_gmt_overlap()

        tau0 = 100.
        optimal_tau = np.zeros_like(self.commitment_parameter,dtype="float")
        for i,cparam in enumerate(self.commitment_parameter):

            cl = self.sl_contributor(cparam,self.temp_anomaly_year)
            sl_function = cl.calc_contribution
            optimal_tau[i], pcov= leastsq(self.calc_residuals, tau0, args=(sl_function))

            print "cparam=",cparam," tau=", optimal_tau[i]

        self.calibrated_parameter = optimal_tau
        return [self.commitment_parameter,optimal_tau]



    # def calc_contribution(self,gmt_anomaly):

    #     """ calculates sea level response based on fitted parameters and
    #         external temperature input. calibrate() has to be run before.
    #     """

    #     for j,cparam in enumerate(self.commitment_parameter):
    #         tau = self.calibrated_parameter[j]

    #         cl = self.sl_contributor(cparam,self.temp_anomaly_year)
    #         sl_calculated = da.DimArray(cl.calc_contribution(gmt_anomaly,tau),
    #                             axes=gmt_anom.time,dims="time")

    #     return