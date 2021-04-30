#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Lingwei Meng (lingweim@princeton.edu)
import numpy as np

# physical parameters
c_p = 1005.7    # specific heat at constant pressure for dry air (J * kg^-1* k^-1)
Rd = 287.05     # gas constant for dry air (J * kg^-1 * k^-1)
Rv = 461.51     # gas constant for water vapor (J * kg^-1 * k^-1)
Lv = 2.555e6    # latent heat of vaporization (J * kg^-1)
epsilon = Rd/Rv # ~ 0.622


def mixing_ratio_from_q(q):
    ''' 
    Calculate maxing ratio (in kg * kg^-1) 
    given specific humidity (q; in kg * kg^-1). 
    
    '''
    return q / (1-q)

def mixing_ratio_from_e(p, e):
    '''
    Calculate mixing ratio
    given air pressure (p; in Pa) and water vapor pressure (e; in Pa).
    
    '''
    return epsilon * e / (p-e)


def saturation_vapor_pressure(T):
    '''
    Calculate the saturation vapor pressure (in Pa) 
    given temperature (T; in Kelvin) 
    from the Clausius–Clapeyron relation (the August–Roche–Magnus formula).
    
    '''
    return 610.94 * np.exp(17.625 * (T-273.15) / (T-273.15 + 243.04))


def vapor_pressure_by_mixing_ratio(p, r_v):
    '''
    Calculate water vapor pressure (in Pa) 
    given air pressure (in Pa) and mixing ratio (in kg kg**-1).
    
    '''
    return p * r_v / (epsilon+r_v)

