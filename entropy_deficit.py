#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Lingwei Meng (lingweim@princeton.edu)
import numpy as np, xarray as xr
from meteocalc import mixing_ratio_from_q,\
                      mixing_ratio_from_e,\
                      saturation_vapor_pressure,\
                      vapor_pressure_by_mixing_ratio


# physical parameters
c_p = 1005.7    # specific heat at constant pressure for dry air (J kg**-1 K**-1)
Rd = 287.05     # gas constant for dry air (J kg**-1 K**-1)
Rv = 461.51     # gas constant for water vapor (J kg**-1 K**-1)
Lv = 2.555e6    # latent heat of vaporization (J kg**-1)
epsilon = Rd/Rv # ~ 0.622


def pseudoadiabatic_entropy(T, p, RH=None, q=None):
    '''
    Calculate the pseudoadiabatic moist entropy from Bryan (2008).
    The formulation is: s = c_p*log(T) - Rd*log(p-e) + Lv*r_v/T - Rv*r_v*log(RH)

    Inputs
    ----------
    T : air temperature (in Kelvin)
    p : air pressure (in Pa)
    RH: relative humidity (unitless)
    q : specific humidity (in kg kg**-1)

    Note: At least one of the two variables of RH and q must be provided.

    Required calculations
    -------
    e  : water vapor pressure (in Pa)
    r_v: mixing ratio (in kg kg**-1)

    Returns
    -------
    s : pseudoadiabatic moist entropy (in J kg**-1 K**-1)

    '''
    if RH is None:
        # q is given by the inputs
        assert q is not None, 'At least one of the two variables must be provided: relative humidity or specific humidity.'
        r_v = mixing_ratio_from_q(q)
        e = vapor_pressure_from_mixing_ratio(p, r_v)
        RH = e / saturation_vapor_pressure(T)
    else:
        # RH is given by the inputs
        e = saturation_vapor_pressure(T) * RH
        r_v = mixing_ratio_from_e(p, e)

    s = c_p*np.log(T) - Rd*np.log(p-e) + Lv*r_v/T - Rv*r_v*np.log(RH)
    s.attrs['long_name'] = 'Pseudoadiabatic moist entropy from Bryan (2008).'
    s.attrs['units'] = 'J/K/kg'

    return s


def entropy_deficit(slp, sst, Tb, RHb, Tm, RHm, p_b = 950e2, p_m = 600e2):
    '''
    Calculate the entropy deficit chi_m from Tang and Emanuel (2012)

    Parameters
    ----------
    slp : sea level pressure (in Pa)
    sst : sea surface temperature (in Kelvin)
    p_b : the pressure of the boundary layer (in Pa).
          the default value is 95000 Pa.
    Tb  : the temperature of the boundary layer (in Kelvin)
    RHb : the relative humidity of the boundary layer (unitless)
    p_m : the pressure of the midlevel (in Pa).
          the default value is 60000 Pa (as in Tang and Emanuel (2012)).
    Tm  : the temperature of the midlevel (in Kelvin)
    RHm : the relative humidity of the midlevel (unitless)

    Returns
    -------
    chi_m : the entropy deficit (unitless)

    '''
    # Saturation moist entropy at mid-level (600hPa) in the inner core of the TC
    s_m_sat   = pseudoadiabatic_entropy(T=Tm, p=p_m, RH=1)
    # Moist entropy at mid-level (600hPa) in the environment
    s_m       = pseudoadiabatic_entropy(T=Tm, p=p_m, RH=RHm)
    # Saturation moist entropy at sea surface temperature
    s_sst_sat = pseudoadiabatic_entropy(T=sst,p=slp, RH=1)
    # Moist entropy for boundary layer (950 hPa)
    s_b       = pseudoadiabatic_entropy(T=Tb, p=p_b, RH=RHb)

    # Entropy deficit
    chi_m = ((s_m_sat - s_m) / (s_sst_sat - s_b)).pipe(lambda x: x.where(x>0)).pipe(lambda x: x.where(x<2))
    chi_m.attrs['long_name'] = 'entropy deficit from Tang and Emanuel (2012)'

    return chi_m
