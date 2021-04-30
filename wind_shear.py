#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Lingwei Meng (lingweim@princeton.edu)
import numpy as np
import xarray as xr

def wind_shear(u850, v850, u200, v200):
    '''
    Calculate 200-850 hPa Vertical Wind Shear (vws; in m s**-1).

    '''
    vws = ((u200 - u850)**2 + (v200 - v850)**2) ** 0.5

    vws.attrs['long_name'] = 'Vertical wind shear between 200hPa and 850hPa'
    vws.attrs['units'] = 'm/s'

    return vws

