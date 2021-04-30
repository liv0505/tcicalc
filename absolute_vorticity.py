#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Lingwei Meng (lingweim@princeton.edu)
import xarray as xr
import numpy as np 

def absolute_vorticity(vor, lat):
    '''
    calculate absolute vorticity (in s**-1)
    given the relative vorticity (vor; in s**-1) 
    and the corresponding latitudes (degrees_north)
    
    '''
    eta = vor + 2 * (2*np.pi/(24*3600)) * np.sin(lat*np.pi/180)
    eta.attrs['long_name'] = 'Absolute vorticity at 850hPa.'
    eta.attrs['units'] = 's**-1'

    return eta