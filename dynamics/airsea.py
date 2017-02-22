# -*- coding: utf-8 -*-
"""
Atlantis dynamics: air-sea interaction
======================================

Atlantis is a Python library for atmospheric, oceanographic and
hydrographic data analysis and visualization.

This is a set of function and classes to help with air-sea interaction.

TODO
----


Disclaimer
----------
This software may be used, copied, or redistributed as long as it is not
sold and this copyright notice is reproduced on each copy made. This
routine is provided as is without any express or implied warranties
whatsoever.

Authors
-------
Sebastian Krieger (sebastian.krieger@usp.br)

Revision history
----------------
2 (2016-05-30 23:17 -0300)
1 (2013-06-05 14:24 -0300)

"""
from __future__ import division

__version__ = '$Revision: 2 $'
# $Source$

from numpy import (arctan, convolve, ones, exp, log, pi, nan, nanmin, nanmean,
    nanmax, hanning)

from atlantis import data


class BulkFluxes():
    """
    Air-sea bulk turbulent fluxes computation.

    Definitions
    -----------
    Tau : Momentum, wind stress, in N \cdot m^{-2}, or
          kg \cdot m \cdot s^{-2} \cdot m^{-2}
    Q_H : Sensible heat, in W \cdot m^{-2}, or
          J \cdot s^{-1} \cdot m^{-2}
    Q_E : Latent heat, in W \cdot m^{-2}, or
          J \cdot s^{-1} \cdot m^{-2}
    E   : Evaporation, in mg \cdot s^{-1} \cdot m^{-2}


    Given
    -----
    U_O : Observed vector wind at height z_U
    T_O : Observed temperature at height z_T
    q_O : Observed specific humidity at height z_q
    SLP : Sea level pressure
    SST : Sea surface temperature
    SSU : Surface vector current


    TODO
    ----
    * Calculate precipitation fluxes;
    * Calculate long wave radiation fluxes;


    References
    ----------
    Large, W. G. and S. Yeager (2004). Diurnal to decadal global forcing
    for ocean and sea-ice models: The data sets and flux climatologies.
    Technical note NCAR/TN-460+STR, NCAR.

    Large, W. G. (2006). Surface fluxes for practioners of global ocean
    data assimilattion. In E. Chassignet and J. Verron (Eds.), Ocean
    weather and forecasting, pp. 229–270. Heidelberg: Springer.

    """
    ###########################################################################
    # CONSTANTS AND PARAMETERS
    ###########################################################################
    GAMMA = 0.01   # adiabaic lapse rate for dry air
    R_gas = 287.04 # dry air gas constant [J/kg/K]
    KAPPA = 0.4    # von Karman constant
    LAMBDA = 2.5e6 # latent heat of vaporization [J/kg]
    g = 9.816      # gravitational acceleration [m/s**2]
    cp = 1003.5    # specific heat capacity of dry air [J/kg/K]


    ###########################################################################
    # FUNCTIONS
    ###########################################################################
    def __init__(self, dat=None, **kwargs):
        # Resets fields dictionary.
        self.fields = dict()
        # If data is given, updates flux container
        if dat is not None:
            self.update(dat, **kwargs)
        #
        return None


    def __getitem__(self, key):
        return self.fields[key]


    def __setitem__(self, key, value):
        # Assigns value to class' list of items (fields)
        self.fields[key] = value


    def __delitem__(self, key):
        del self.fields[key]


    def __iter__(self):
        return iter(self.fields)


    def update(self, dat, append=False, smooth=False, report=True, **kwargs):
        """
        Updates bulk fluxes.

        Parameters
        ----------
        dat : odv.Sequence()
            ODV data object.
        append : boolean, optional
            Sets whether calculated parameters are appended to data
            sequence. Default is `False`.
        smooth : boolean, optional
            If 'True' (default), smoothes data using a boxcar moving
            average.
        report : boolean, optional

        Returns
        -------
        Nothing

        """
        # Checks data fields and uses the following wind components in order
        # of preference:
        #     . ['alongchannel_wind', 'crosschannel_wind'];
        #     . ['eastward_wind', 'northward_wind']
        fields = dat.fields.keys()
        if ('alongchannel_wind' in fields) & ('crosschannel_wind' in fields):
            _wind = ('alongchannel_wind', 'crosschannel_wind')
        elif ('eastward_wind' in fields) & ('northward_wind' in fields):
            wind = ('eastward_wind', 'northward_wind')
        else:
            raise ValueError('Noooo!!!!!')
        # Assings data to some variables.
        U_O = dat[_wind[0]].data + 1j * dat[_wind[1]].data
        T_O = dat['air_temperature'].data
        SLP = dat['air_pressure_at_sea_level'].data
        SST = dat['sea_water_temperature'].data
        SSU = 0
        try:
            q_O = dat['specific_humidity'].data
        except:
            q_O = self.relative2specific_humidity(
                dat['relative_humidity'].data, SST, SLP)
        # Makes some units conversion
        if dat['air_pressure_at_sea_level'].units != 'hPa':
            raise ValuError('Invalid unit for pressure.')
        #
        units = dict(
            U_O=dat[_wind[0]].units,
            T_O=dat['air_temperature'].units,
            q_O='kg kg-1',
            SLP='hPa',
            SST=dat['sea_water_temperature'].units,
            SSU=dat[_wind[0]].units
        )
        # Smoothes data using moving average.
        if smooth:
            window = ones(29.) #hanning(29.)
            window /= window.sum()
            U_O = convolve(U_O.real, window, mode='same') +  \
                1j * convolve(U_O.imag, window, mode='same')
            T_O = convolve(T_O, window, mode='same')
            q_O = convolve(q_O, window, mode='same')
            SLP = convolve(SLP, window, mode='same')
            SST = convolve(SST, window, mode='same')
            #SSU = convolve(SSU, window, mode='same')
            print 'Smooooooth'
        # Calculates bulk_fluxes
        Tau, QH, QE, dU10, dtheta10, dq10, L, zeta = \
            self.bulk_fluxes(U_O, T_O, q_O, SLP, SST, SSU, units=units,
            result='all', **kwargs)
        # Apppends data to class items.
        self.__setitem__('Tau_x', data.get_standard_variable(
            'surface_downward_northward_stress', data=Tau.imag,
            string_format='{:.4f}'))
        self.__setitem__('Tau_y', data.get_standard_variable(
            'surface_downward_eastward_stress', data=Tau.real,
            string_format='{:.4f}'))
        self.__setitem__('QH',  data.get_standard_variable(
            'surface_downward_sensible_heat_flux', data=QH))
        self.__setitem__('QE', data.get_standard_variable(
            'surface_downward_latent_heat_flux', data=QE))
        #
        CD = self.C_D(dU10)
        self.__setitem__('CD', data.get_standard_variable(
            'surface_drag_coefficient_for_momentum_in_air', data=CD))
        CH = self.C_H(dU10, zeta=zeta)
        self.__setitem__('CH', data.get_standard_variable(
            'surface_drag_coefficient_for_heat_in_air', data=CH))
        CE = self.C_E(dU10)
        self.__setitem__('CE', data.Variable(
            standard_name='surface_drag_coefficient_for_evaporation_in_air',
            data=CE, units='1'))
        #
        if append:
            dat['surface_downward_northward_stress'] = self.__getitem__('Tau_x')
            dat['surface_downward_eastward_stress'] = self.__getitem__('Tau_y')
            dat['surface_downward_sensible_heat_flux'] = self.__getitem__('QH')
            dat['surface_downward_latent_heat_flux'] = self.__getitem__('QE')
            #
            dat['surface_drag_coefficient_for_momentum_in_air'] = \
                self.__getitem__('CD')
            dat['surface_drag_coefficient_for_heat_in_air'] = \
                self.__getitem__('CH')
            dat['surface_drag_coefficient_for_evaporation_in_air'] =\
                self.__getitem__('CE')
        #
        if report:
            print r'\Delta \theta', nanmin(dtheta10), nanmean(dtheta10), nanmax(dtheta10)
        #
        return -1


    def relative2specific_humidity(self, RH, Ta, Pa=101800):
        """Converts relative humidity to specific humidity."""
        return RH / 100 * self.q_sat(Ta, Pa)


    def saturated_water_pressure(self, Ta, Pa=1018.):
        """
        Computes saturated water vapor pressure in air.

        Parameters
        ----------
        Ta : array like
            Air temperature [degC]
        Pa : array like, optional
            Air pressure [hPa]

        Returns
        -------
        ew : array like
            Saturated water pressure [hPa].

        """
        # This uses the computation in the COARE code.
        return (6.1121 * exp(17.502 * Ta / (Ta + 240.97)) *
            (1.0007 + 3.46e-6 * Pa))


    def q_sat(self, Ta, Pa=1018.):
        """
        Computes specific humidity at saturation.

        Parameters
        ----------
        Ta : array like
            Air temperature [degC]
        Pa : array like, optional
            Pressure [Pa].

        Returns
        -------
        q_sat : array like
            Saturation specific humidity [kg kg-1]

        """
        # Calculates saturation water pressure
        ew = self.saturated_water_pressure(Ta, Pa)
        # Returns specific humidity. Remember pressure is in [Pa] and that
        # specific humidity is in [kg kg-1].
        return 0.62197 * (ew / (Pa - 0.378 * ew));


    def pot_temperature(self, T, z):
        """
        Potential temperature.

        Parameters
        ----------
        T : array like
            Temperature.
        z : array like
            Height of temperature.

        Returns
        -------
        theta : array like
            Potential temperature, in degrees Celsius

        """
        return T + self.GAMMA * z


    def virtual_pot_temperature(self, theta, q):
        """
        Virtual potential air temperature.

        Parameters
        ----------
        theta : array like
            Potential temperature.
        q : array like
            Specific humidity.

        Returns
        -------
        theta_V : array like
            Virtual air temperature

        """
        return theta * (1 + 0.608 * q)


    def density_air(self, SLP, theta_V, SLP_unit='hPa', T_unit='degC') :
        """
        Density of air.

        Parameters
        ----------
        SLP : array like
            Sea level pressure.
        theta_V : array like
            Virtual temperature.
        SLP_unit : string, optional
            Sets wether sea level pressure is given in bar ('bar'),
            hectopascal ('hPa', default), millibar ('mbar'),
            kilopascal ('kPa') or standard atmosphere ('atm').
        T_unit : string, optional
            Sets wether sea surface temperature is given in degrees
            Celsius ('degC', default) or Kelvin ('K').

        Returns
        -------
        rho_air : array like :
            Density of air.

        """
        if SLP_unit == 'bar':
            # Converts SLP from bar to Pa
            SLP = SLP * 1e5
        elif SLP_unit in ['hPa', 'mbar']:
            # Converts SLP from hPa to Pa
            SLP = SLP * 1e2
        elif SLP_unit == 'kPa':
            # Converts SLP from kPa to Pa
            SLP = SLP * 1e3
        elif SLP_unit == 'atm':
            # Converts SLP from atm to Pa
            SLP = SLP * 1.01325e5
        elif SLP_unit == 'Pa':
            SLP = SLP * 1.
        else:
            raise Warning, 'Pressure unit %s not implemented yet.' % (SLP_unit)

        if T_unit == 'degC':
            # Converts temperature from degrees Celsius to Kelvin.
            theta_V = theta_V + 273.15

        return SLP / (self.R_gas * theta_V)


    def humidity_sat(self, SST, rho, q1=0.98, q2=640380, q3=-5107.4,
        SST_unit='degC'):
        """
        Parameterized saturated humidity.

        Parameters
        ----------
        SST : array like
            Sea surface temparature.
        rho : array like
            Air density, in kg*m**(-3).
        q1, q2 : float, optional
            Specicic coefficients for sea-water.
        SST_unit : string, optional
            Sets wether sea surface temperature is given in degrees
            Celsius ('degC', default) or Kelvin ('K').

        Returns
        -------
        q_sat : array like
            Saturated humidity over seawater, in kg/kg

        """
        if SST_unit == 'degC':
            # Converts SST from degrees Celsius to Kelvin.
            SST = SST + 273.15

        return q1 * q2 / rho  * exp(q3 / SST)


    def u_star(self, dU, CD):
        """
        Calculates friction velocity u*.

        Paramaters
        ----------
        dU : array like
            Difference between wind velocity and sea surface velocity,
            in m/s.
        CD : array like
            Drag coefficient.

        Returns
        -------
        u_star : array like
            Friction velocity.

        """
        return (CD ** 0.5) * dU


    def t_star(self, dtheta, CD, CH):
        """
        Calculates turbulent fluctuations of potential temperature.

        Parameters
        ----------
        dtheta : array like
            Difference between potential temperature and sea surface
            temperature, in degrees Celsius.
        CD : array like
            Drag coefficient.
        CH : array like
            Stanton number.

        Returns
        -------
        theta_star : array like
            Turbulent fluctuations of potential temperature.

        """
        return CH / (CD ** 0.5) * dtheta


    def q_star(self, dq, CD, CE):
        """
        Calculates turbulent fluctuations of specific humidity.

        Parameters
        ----------
        dq : array like
            Difference between specific humidity and saturated humidity
            at sea surface, in mg/kg.
        CD : array like
            Drag coefficient.
        CE : array like
            Dalton number.

        Returns
        -------
        q_star : array like
            Turbulent fluctuations of specific humidity.

        """
        return CE / (CD ** 0.5) * dq


    def stability_atmosphere(self, z, L):
        """
        Atmospheric stability.

        Parameters
        ----------
        z : array like
            Height, in m.
        L : array like
            Monin-Obhukhov length

        Returns
        -------
        zeta : array like
            Atmospheric stability.

        """
        return z / L


    def Monin_Obukhov_Length(self, u_star, theta_star, q_star, theta_V, q):
        """
        Calculates the Monin-Obukhov length.

        The length is used to describe the effects of buoyancy on turbulent
        flows, particularly in the lower tenth of the atmospheric boundary
        layer.

        Parameters
        ----------
        u_star : array like
            Frinction velocity u*, in m/s.
        theta_star : array like
            Scaling temperature, in degrees Celsius.
        q_star : array like
            Scaling specific humidity, in mg/kg.
        theta_V : array like
            Virtual potential temperature, in degrees Celsius.
        q : array like
            Specific humidity, in mg/kg.

        Returns
        -------
        L : array like
            Monin-Obukhov length

        References
        ----------
        Large, W. G. (2006). Surface fluxes for practioners of global
        ocean data assimilattion. In E. Chassignet and J. Verron (Eds.),
        Ocean weather and forecasting, pp. 229–270. Heidelberg:
        Springer.

        http://en.wikipedia.org/wiki/Monin%E2%80%93Obukhov_length

        """
        # Converts virtual potential temperature from degres Celsius to Kelvin
        theta_V = theta_V + 273.15
        B0 = self.g * (theta_star / theta_V + q_star / (q + 0.608**(-1)))
        return u_star**2 / (self.KAPPA * B0)


    def Psi(self, zeta, result='both'):
        """
        Empirical stability functions.

        The returned values of Psi_M and Psi_S are used to bring
        observational measurements of wind speed, potential temperature and
        humidity from non-netural profiles to neutral profiles. They are the
        integrals of the dimensionless flux profiles of momentum,
        \Psi_M(\zeta), and of the scalars heat and moisture, \Psi_S(\zeta).

        Parameters
        ----------
        zeta : array like
            Atmospheric stability.
        result : string
            Determines if either the 'momentum', 'scalar' or 'both' (default)
            dimensionless flux profiles are returned.

        Returns
        -------
        Psi_M, Psi_S : array like
            According to the 'result' parameters, function returns
            either \Psi_M(\zeta), and/or \Psi_S(\zeta).

        Reference
        ---------
        Paulson, C. A., 1970: The Mathematical Representation of Wind
        Speed and Temperature Profiles in the Unstable Atmospheric
        Surface Layer. J. Appl. Meteor., 9, 857–861.

        """
        mask_stable = (zeta >= 0)
        mask_unstable = (zeta < 0)
        # Initializes variables
        Psi_M = zeta * 0
        Psi_S = zeta * 0
        # Calculates Psi_M and Psi_S for the stable case (zeta > 0)
        Psi_M += -5 * zeta * mask_stable
        Psi_S += -5 * zeta * mask_stable
        # Calculates Psi_M and Psi_S for the unstable case (zeta < 0). It is
        # important to note that X = (1 - 16 * zeta) ** 0.25, but since zeta < 0,
        # we use absolute values and then mask them out to ensure proper
        # calculations.
        X = (1 + 16 * abs(zeta)) ** 0.25
        Psi_M += 2 * log((1 + X**2) / 2) * mask_unstable
        Psi_S += (log((1 + X**2) / 2) + log((1 + X) / 2) - 2 * arctan(X) +
            0.5 * pi) * mask_unstable
        if result == 'both':
            return Psi_M, Psi_S
        elif result == 'momentum':
            return Psi_M
        elif result == 'scalar':
            return Psi_S
        else:
            raise ValueError('Invalid result type `{}`.'.format(result))


    def Q_E(self, E):
        """
        Converts evaporation to latent heat of evaporation

        Parameters
        ----------
        E : array like
            Evaporation, in mg \cdot s^{-1} \cdot m^{-2}.

        Returns
        -------
        Q_E : array like
            Latent heat of evaporation, in W \cdot m^{-2}.

        """
        return self.LAMBDA * E


    def law_of_the_wall(self, u_star, z, z0, Psi=0):
        """
        Law of the wall for wind, temperature or humidity profiles.

            DU(z) = (u_star / KAPPA) * ln(z / z0 - Psi)

        Parameters
        ----------
        u_star : float
            Depending on the application, this parameter may be the
            friction velocity (in m/s), scaling temperature (in K or
            degC), scaling specific humidity (in mg/kg).
        z : array like
            Height above sea level (in m).
        z0 : float
            Roughness length (in m).
        Psi : float
            Empirical correction for stability. Neutral profiles have
            Psi=0.

        Returns
        -------
        DU : array like
            Either wind velocity, temperature or specific humidity at
            height z.

        """
        return u_star / self.KAPPA * (log(z / z0) - Psi)


    def C_D(self, dU):
        """
        Calculates the drag coefficient C_D using multiple regression
        parameters.

        Parameters
        ----------
        dU : array like
            Difference between wind and sea surface velocities,
            in m/s.

        Returns
        -------
        CD : array like
            Drag coefficient.

        """
        # CD = u_star**2 / dU ** 2
        a1 = 0.00270
        a2 = 0.000142
        a3 = 0.0000764
        return a1 / dU + a2 + a3 * dU


    def C_E(self, dU):
        """
        Calculates the Dalton number.

        Parameters
        ----------
        dU : array like
            Difference between wind and sea surface velocities, in m/s.

        Returns
        -------
        CE : array like
            Dalton number.

        """
        return 0.0346 * self.C_D(dU) ** 0.5


    def C_H(self, dU, zeta=0):
        """
        Calculates the Stanton number.

        Parameters:
        dU : array like
            Difference between wind and sea surface velocities, in m/s.
        zeta : array like
            Atmospheric stability.

        Returns
        -------
        CH : array like
            Stanton number

        """
        return (0.0180 * self.C_D(dU) ** 0.5 * (zeta >= 0) +
            0.0327 * self.C_D(dU) ** 0.5 * (zeta < 0))


    def bulk_fluxes(self, U_O, T_O, q_O, SLP, SST, SSU, z_U=10, z_T=10,
        z_q=10, N=5, result='fluxes', units=dict()):
        """
        Computes bulk air-sea fluxes.

        Parameters
        ----------
        U_O : array like
            Observed vector wind at height z_U, in m/s. Accepts
            velocity given in complex notation.
        T_O : array like
            Observed temperature at height z_T, in degrees Celsius
        q_O : array like
            Observed specific humidity at height z_q, in kg/kg.
        SLP : array like
            Sea level pressure, in hectopascal (hPa).
        SST : array like
            Sea surface temperature, in degrees Celsius.
        SSU : array like
            Surface vector current, in m/s. Accepts velocity given in
            complex notation.
        z_U, z_T, z_q : float, array like, optional
            Height of the observed vector wind, temperature and
            specific humidity, in m. Default value is 10 m.
        N : integer, optional
            Number of iterations to calculate bulk parameters.
        result : string, optional
            Determines if either only 'fluxes' (default) or if 'all'
            bulk calculations are returned.
        units : dictionary, optional
            Sets the units for the variables.

        Returns
        -------
        Tau, Q_H, Q_E : array like
            Momentum, sensible and latent heat fluxes.
        dU10, dtheta10, dq10, L, zeta: array like, optional :
            If `result` is set to 'all', additionally to the fluxes, ...

        """
        # STEP ZERO: Checks for proper units.
        default_units = dict(U_O='m s-1', T_O='degC', q_O='kg kg-1', SLP='hPa',
            SST='degC', SSU='m s-1', z_U='m', z_T='m', z_q='m')
        units = units or default_units
        #
        if units['U_O'] != 'm s-1' or units['SSU'] != 'm s-1':
            raise Warning, 'Wind speed or surface current speed unit invalid.'
        if  units['T_O'] == 'K':
            T_O = T_O - 273.15
        if units['SST'] == 'K':
            SST = SST - 273.15
        if units['q_O'] == 'g kg-1':
            q_O = q_O / 1000
        if units['SLP'] == 'bar':
            SLP = SLP * 1e3             # From bar to hPa
        elif units['SLP'] == 'Pa':
            SLP = SLP *  1e-2           # From Pa to hPa
        elif units['SLP'] == 'kPa':
            SLP = SLP *  10             # From kPa to hPa
        elif units['SLP'] == 'atm':
            SLP = SLP *  1.01325e3      # From atm to hPa
        units = default_units

        # FIRST STEP: Assume theta(z_u) = theta(z_theta) and q(z_u) = q(z_q),
        # compute potential temperature and virtual potential temperature, sea
        # surface humidity, difference between wind speed and surface current
        # speed, difference between potential temperature and sea surface
        # temperature, observed specific humidity and sea surface humidity.
        theta_O = self.pot_temperature(T_O, z_T)
        theta_O_v = self.virtual_pot_temperature(theta_O, q_O)
        SSq = self.humidity_sat(SST, self.density_air(SLP, theta_O_v))
        DU = U_O - SSU
        dU = abs(DU)
        dtheta = theta_O - SST
        dq = q_O - SSq

        # SECOND STEP: Assume neutral stability (\zeta = 0) and observations at
        # z=10m. Calculate the transfer coefficients and initial turbulent scales.
        zeta = 0
        dU10, dtheta10, dq10 = dU, dtheta, dq
        for i in range(N):
            # Some checks first:
            try:
                dU10[dU10 < 1] = 1.
            except:
                pass

            CD = self.C_D(dU10)
            CH = self.C_H(dU10, zeta=zeta)
            CE = self.C_E(dU10)

            US = self.u_star(dU10, CD)
            thetaS = self.t_star(dtheta10, CD, CH)
            qS = self.q_star(dq10, CD, CE)

            q = SSq + dq10
            theta_v = SST + dtheta10 + 0.608 * q
            L = self.Monin_Obukhov_Length(US, thetaS, qS, theta_v, q)
            zeta = self.stability_atmosphere(10., L)

            dU10 = dU - self.law_of_the_wall(US, z_U, 10.,
                Psi=self.Psi(zeta*z_U/10., result='momentum'))
            dtheta10 = dtheta - self.law_of_the_wall(thetaS, z_T, 10.,
                Psi=self.Psi(zeta*z_T/10., result='scalar'))
            dq10 = dq - self.law_of_the_wall(qS, z_q, 10.,
                Psi=self.Psi(zeta*z_q/10., result='momentum'))

            #print US, thetaS, qS, dU10, dtheta10, dq10, L, zeta

        # THIRD STEP: Compute the bulk turbulent fluxes
        rho = self.density_air(SLP, theta_v)
        tau = rho * US ** 2
        Tau = tau * DU / dU
        QH = rho * self.cp * US * thetaS
        E = rho * US * qS
        E[E > 0] = nan
        QE = self.Q_E(E)
        if result == 'fluxes':
            return Tau, QH, QE
        elif result == 'all':
            return Tau, QH, QE, dU10, dtheta10, dq10, L, zeta
        else:
            raise ValueError('Invalid result type `{}`.'.format(result))
