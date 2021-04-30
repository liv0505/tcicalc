#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Lingwei Meng (lingweim@princeton.edu)
from tcpi import pcmin, pcmin3_kflag, cape
import numpy as np, xarray as xr
from meteocalc import mixing_ratio_from_q


def potential_intensity(sst, slp, plevel, T, q, reverse_p=True):
    '''
    Calculate the potential intensity by calling pcmin.f90 adopted from Prof. Emanuel.
    See Emanuel's original Fortran versions: ftp://texmex.mit.edu/pub/emanuel/TCMAX/

    Parameters
    ----------
    sst    : sea surface temperature (in Kelvin);
    slp    : sea level pressure (in Pa);
    plevel : pressure levels ;
    T      : air temperature (in Kelvin);
    q      : specific humidity (in kg kg**-1);


    Returns
    -------
    vmax : the maximum wind speed achievable in tropical cyclones (in m s**-1)
    pmin : the minimum central pressure achievable in tropical cyclones (in hPa)
    iflag: a value of 1 means OK;
           a value of 0 indicates no convergence (hypercane);
           a value of 2 means that the CAPE routine failed.
    '''
    r_v = mixing_ratio_from_q(q)

    if reverse_p:
        # In pcmin.f90, the arrays MUST be arranged with increasing index
        # corresponding to decreasing pressure.
        plevel = plevel.isel(level=slice(-1,None,-1))
        T = T.isel(level=slice(-1, None, -1))
        r_v = r_v.isel(level=slice(-1, None, -1))

    pmin, vmax, iflag = xr.apply_ufunc(pcmin3_kflag,
        sst, slp, plevel, T, r_v,
        input_core_dims=[sst.dims, slp.dims, plevel.dims, T.dims, r_v.dims],
        output_core_dims=[sst.dims, sst.dims, []],
        vectorize=True)

    pmin.attrs['long_name'] = 'mininum central pressure'
    pmin.attrs['units'] = 'hPa'
    vmax.attrs['long_name'] = 'maximum surface wind speed'
    vmax.attrs['units'] = 'm/s'
    iflag = iflag.astype('int32')
    iflag.attrs['long_name'] = '1: OK; 0: no convergence; 2: CAPE routine failed.'

    # wrap the pmin, vmax and iflag into an xarray.Dataset
    PI = xr.Dataset(dict(pmin=pmin, vmax=vmax, iflag=iflag))

    return PI
