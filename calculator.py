#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Lingwei Meng (lingweim@princeton.edu)
import xarray as xr, numpy as np
import os.path
from entropy_deficit import entropy_deficit
from potential_intensity import potential_intensity
from wind_shear import wind_shear
from absolute_vorticity import absolute_vorticity



def calc_PI(year, datadir=None, outdir=None):
    # ------------------------------------------------------------------
    # the potential intensity PI
    # ------------------------------------------------------------------

    datadir = './' if datadir is None else datadir
    outdir = './' if outdir is None else outdir

    # [load sst]
    ifile = os.path.join(datadir, f'era5_sst_monthly_{year:04d}.nc')
    sst0 = xr.open_dataarray(ifile)

    # [load sea level pressure]
    ifile = os.path.join(datadir, f'era5_slp_monthly_{year:04d}.nc')   
    slp0 = xr.open_dataarray(ifile)

    # [load or make the landmask]
    ofile = os.path.join(datadir, f'era5_oceanarea.nc')
    if os.path.exists(ofile):
        ocean_area = xr.open_dataarray(ofile)
    else:
        ocean_area = sst0.isel(time=0).drop('time').pipe(lambda x: x*0==0)
        ocean_area.to_dataset(name='ocean_area').to_netcdf(ofile)
        print('The landmask for ERA5 dataset is saved to:', ofile)

    # [load temperature]
    ifile = os.path.join(datadir, f'era5_plevels_monthly_{year:04d}.nc')
    T0 = xr.open_dataset(ifile).t

    # [load specfic humidity]
    q0 = xr.open_dataset(ifile).q

    print('*** Calculating the potential intensity ***')
    PI_list = []
    ntime = len(sst0.time)
    for it in range(ntime):
        print(f'Time:{it+1:02d}/{ntime}')
        iPI = potential_intensity(
            sst = sst0.isel(time = it),
            slp = slp0.isel(time = it).where(ocean_area),
            plevel = T0.level,
            T = T0.isel(time = it).where(ocean_area),
            q = q0.isel(time = it).where(ocean_area),
            reverse_p = True
            )
        PI_list.append(iPI)
    PI = xr.concat(PI_list, dim='time')

    # save the results to the output directory.
    ofile = os.path.join(outdir, f'PI.monthly.{year}.nc' )
    encoding = {dname:{'dtype': 'float32', 'zlib': True, 'complevel': 1}
            for dname in ('pmin', 'vmax')}
    encoding['iflag'] = {'dtype': 'int32'}
    PI.to_netcdf(ofile, encoding=encoding, unlimited_dims='time')

    print('The output for PI has been saved to:', ofile)





def calc_chim(year, p_b=950e2, p_m=600e2, datadir=None, outdir=None):
    # ------------------------------------------------------------------
    # the entropy deficit (chim; unitless) (Tang and Emanuel 2012):
    # (s_m_star - s_m) / (s_sst_star - s_b)
    # ------------------------------------------------------------------

    datadir = './' if datadir is None else datadir
    outdir = './' if outdir is None else outdir

    # [load sst]
    ifile = os.path.join(datadir, f'era5_sst_monthly_{year:04d}.nc')
    sst0 = xr.open_dataarray(ifile)

    # [load sea level pressure]
    ifile = os.path.join(datadir, f'era5_slp_monthly_{year:04d}.nc')   
    slp0 = xr.open_dataarray(ifile)

    # [load or make the landmask]
    ofile = os.path.join(datadir, f'era5_oceanarea.nc')
    if os.path.exists(ofile):
        ocean_area = xr.open_dataarray(ofile)
    else:
        ocean_area = sst0.isel(time=0).drop('time').pipe(lambda x: x*0==0)
        ocean_area.to_dataset(name='ocean_area').to_netcdf(ofile)
        print('The landmask for ERA5 dataset is saved to:', ofile)
        
    # [load temperature]
    ifile = os.path.join(datadir, f'era5_plevels_monthly_{year:04d}.nc')
    T0 = xr.open_dataset(ifile).t

    # [load relative humidity]
    RH0 = xr.open_dataset(ifile).r / 100
    print('*** Calculating the entropy deficit ***')
    chim_list = []
    ntime = len(sst0.time)
    for it in range(ntime):
        print(f'Time:{it+1:02d}/{ntime}')
        ichim = entropy_deficit(
            slp = slp0.isel(time = it).where(ocean_area),
            sst = sst0.isel(time = it),
            Tb = T0.isel(time = it).interp(level=p_b/100).where(ocean_area),
            RHb = RH0.isel(time = it).interp(level=p_b/100).where(ocean_area),
            Tm = T0.isel(time = it).interp(level=p_m/100).where(ocean_area),
            RHm = RH0.isel(time = it).interp(level=p_m/100).where(ocean_area),
            p_b = p_b, p_m = p_m
            )
        chim_list.append(ichim)
    chim = xr.concat(chim_list, dim='time')

    # save the results to the output directory.
    ofile = os.path.join(outdir, f'chim.monthly.{year}.nc' )
    encoding = {'chim':{'dtype': 'float32', 'zlib': True, 'complevel': 1}}
    chim.to_dataset(name='chim') \
        .to_netcdf(ofile, encoding=encoding, unlimited_dims='time')
    print('The output for entropy deficit has been saved to:', ofile)




def calc_vws(year, datadir=None, outdir=None):
    # ------------------------------------------------------------------
    # 850 and 200hPa vertical wind shear (vws; in m/s)
    # ------------------------------------------------------------------

    datadir = './' if datadir is None else datadir
    outdir = './' if outdir is None else outdir

    # [load or make the landmask]
    ofile = os.path.join(datadir, f'era5_oceanarea.nc')
    if os.path.exists(ofile):
        ocean_area = xr.open_dataarray(ofile)
    else:
        ifile = os.path.join(datadir, f'era5_sst_monthly_{year:04d}.nc')
        sst0 = xr.open_dataarray(ifile)
        ocean_area = sst0.isel(time=0).drop('time').pipe(lambda x: x*0==0)
        ocean_area.to_dataset(name='ocean_area').to_netcdf(ofile)
        print('The landmask for ERA5 dataset is saved to:', ofile)

    # [load u and v]
    ifile = os.path.join(datadir, f'era5_winds_monthly_{year:04d}.nc')
    u850 = xr.open_dataset(ifile).u.sel(level=850).where(ocean_area)
    u200 = xr.open_dataset(ifile).u.sel(level=200).where(ocean_area)
    v850 = xr.open_dataset(ifile).v.sel(level=850).where(ocean_area)
    v200 = xr.open_dataset(ifile).v.sel(level=200).where(ocean_area)

    print('*** Calculating the vertical wind shear ***')
    vws = wind_shear(u850, v850, u200, v200)

    # save the results to the output directory.
    ofile = os.path.join(outdir, f'vws.monthly.{year}.nc' )
    encoding = {'vws': {'dtype': 'float32', 'zlib': True, 'complevel': 1}}
    vws.to_dataset(name = 'vws') \
       .to_netcdf(ofile, encoding=encoding, unlimited_dims='time')
    print('The output for vertical wind shear has been saved to:', ofile)




def calc_VI(year, datadir=None, outdir=None):
    # ------------------------------------------------------------------
    # Ventilation index (VI; unitless) (Tang and Emanuel 2012):
    # vws * chim / vmax
    # ------------------------------------------------------------------

    datadir = './' if datadir is None else datadir
    outdir = './' if outdir is None else outdir

    # [load or calculate the vertical wind shear]
    ifile = os.path.join(outdir, f'vws.monthly.{year:04d}.nc')
    if not os.path.exists(ifile):
        calc_vws(year=year, datadir=datadir, outdir=outdir)
    vws = xr.open_dataarray(ifile)

    # [load or calculate the entropy deficit]
    ifile = os.path.join(outdir, f'chim.monthly.{year:04d}.nc')
    if not os.path.exists(ifile):
        calc_chim(year=year, p_b=950e2, p_m=600e2, datadir=datadir, outdir=outdir)
    chim = xr.open_dataarray(ifile)

    # [load or calculate the potential intensity]
    ifile = os.path.join(outdir, f'PI.monthly.{year:04d}.nc')
    if not os.path.exists(ifile):
        calc_PI(year=year, datadir=datadir, outdir=outdir)
    vmax = xr.open_dataset(ifile).vmax

    print('*** Calculating the ventilation index ***')
    VI = (vws * chim / vmax).pipe(lambda x: x.where(x>0)) # only keep the positive values

    # save the results to the output directory.
    ofile = os.path.join(outdir, f'VI.monthly.{year}.nc')
    encoding = {'VI': {'dtype': 'float32', 'zlib': True, 'complevel': 1}}
    VI.to_dataset(name='VI').to_netcdf(ofile, encoding=encoding,unlimited_dims='time')
    print('The output for ventilation index has been saved to:', ofile)




def calc_eta(year, datadir=None, outdir=None):
    # ------------------------------------------------------------------
    # absolute vorticity (eta; in s**-1)
    # ------------------------------------------------------------------

    datadir = './' if datadir is None else datadir
    outdir = './' if outdir is None else outdir

    # [load or make the landmask]
    ofile = os.path.join(datadir, f'era5_oceanarea.nc')
    if os.path.exists(ofile):
        ocean_area = xr.open_dataarray(ofile)
    else:
        ifile = os.path.join(datadir, f'era5_sst_monthly_{year:04d}.nc')
        sst0 = xr.open_dataarray(ifile)
        ocean_area = sst0.isel(time=0).drop('time').pipe(lambda x: x*0==0)
        ocean_area.to_dataset(name='ocean_area').to_netcdf(ofile)
        print('The landmask for ERA5 dataset is saved to:', ofile)

    # [load the relative vorticity]
    ifile = os.path.join(datadir, f'era5_vor_monthly_{year:04d}.nc')
    vor = xr.open_dataarray(ifile).where(ocean_area)

    print('*** Calculating the absolute vorticity ***')
    eta = absolute_vorticity(vor=vor, lat=vor.latitude)

    # save the results to the output directory.
    ofile = os.path.join(outdir, f'eta850.monthly.{year}.nc')
    encoding = {'eta': {'dtype': 'float32', 'zlib': True, 'complevel': 1}}
    eta.to_dataset(name='eta').to_netcdf(ofile, encoding=encoding,unlimited_dims='time')
    print('The output for absolute vorticity at 850hPa has been saved to:', ofile)



def calc_GPI04(year, datadir=None, outdir=None):
    # ------------------------------------------------------------------
    # GPI04 (Emanuel and Nolan 2004):
    # |10**5*eta|**(3/2) * (RH/50)**3 * (vmax/70)**3 * (1+0.1*vws)**(-2)
    # ------------------------------------------------------------------

    datadir = './' if datadir is None else datadir
    outdir = './' if outdir is None else outdir

    # [load the relative humidity at 700hPa]
    ifile = os.path.join(datadir, f'era5_plevels_monthly_{year:04d}.nc')
    RH = xr.open_dataset(ifile).r.sel(level = 700).drop('level')

    # [load or calculate the absolute vorticity at 850hPa]
    ifile = os.path.join(outdir, f'eta850.monthly.{year:04d}.nc')
    if not os.path.exists(ifile):
        calc_eta(year=year, datadir=datadir, outdir=outdir)
    eta = xr.open_dataarray(ifile)

    # [load or calculate the potential intensity]
    ifile = os.path.join(outdir, f'PI.monthly.{year:04d}.nc')
    if not os.path.exists(ifile):
        calc_PI(year=year, datadir=datadir, outdir=outdir)
    vmax = xr.open_dataset(ifile).vmax

    # [load or calculate the vertical wind shear]
    ifile = os.path.join(outdir, f'vws.monthly.{year:04d}.nc')
    if not os.path.exists(ifile):
        calc_vws(year=year, datadir=datadir, outdir=outdir)
    vws = xr.open_dataarray(ifile)

    print('*** Calculating the GPI (Emanuel and Nolan 2004) ***')
    GPI04 = (1e5 * np.absolute(eta) )**(3/2) \
        * (RH/50)**3 \
        * (vmax/70)**3 \
        * (1+0.1*vws)**(-2)
    GPI04.attrs['long_name'] = 'Genesis Potential Index (Emanuel and Nolan 2004).'
    GPI04.attrs['history'] = '|10**5*eta|**(3/2) * (RH/50)**3 * (vmax/70)**3 * (1+0.1*vws)**(-2)'

    # save the results to the output directory.
    ofile = os.path.join(outdir, f'GPI04.monthly.{year}.nc')
    encoding = {'GPI04': {'dtype': 'float32', 'zlib': True, 'complevel': 1}}
    GPI04.to_dataset(name='GPI04') \
         .to_netcdf(ofile, encoding=encoding, unlimited_dims='time')
    print('The output for GPI (Emanuel and Nolan 2004) has been saved to:', ofile)



def calc_GPI10(year, datadir=None, outdir=None):
    # ------------------------------------------------------------------
    # GPI10 (m**-2 s**-1)(Emanuel 2010):
    # |eta|**3 * chi**(-4/3) * max((vmax-35m/s),0)**2 * (25m/s+vws)**(-4)
    # ------------------------------------------------------------------

    datadir = './' if datadir is None else datadir
    outdir = './' if outdir is None else outdir

    # [load or calculate the absolute vorticity at 850hPa]
    ifile = os.path.join(outdir, f'eta850.monthly.{year:04d}.nc')
    if not os.path.exists(ifile):
        calc_eta(year=year, datadir=datadir, outdir=outdir)
    eta = xr.open_dataarray(ifile)

    # [load or calculate the entropy deficit]
    ifile = os.path.join(outdir, f'chim.monthly.{year:04d}.nc')
    if not os.path.exists(ifile):
        calc_chim(year=year, p_b=950e2, p_m=600e2, datadir=datadir, outdir=outdir)
    chim = xr.open_dataarray(ifile)

    # [load or calculate the potential intensity]
    ifile = os.path.join(outdir, f'PI.monthly.{year:04d}.nc')
    if not os.path.exists(ifile):
        calc_PI(year=year, datadir=datadir, outdir=outdir)
    vmax = xr.open_dataset(ifile).vmax

    # [load or calculate the vertical wind shear]
    ifile = os.path.join(outdir, f'vws.monthly.{year:04d}.nc')
    if not os.path.exists(ifile):
        calc_vws(year=year, datadir=datadir, outdir=outdir)
    vws = xr.open_dataarray(ifile)

    print('*** Calculating the GPI (Emanuel 2010) ***')
    GPI10 = np.absolute(eta)**3 \
        * chim**(-4/3) \
        * (vmax - 35).clip(min=0)**2 \
        * (25 + vws)**(-4)
    GPI10.attrs['long_name'] = 'Genesis Potential Index (Emanuel 2010).'
    GPI10.attrs['history'] = '|eta|**3 * chi**(-4/3) * max((vmax-35),0)**2 * (25+vws)**(-4)'

    # save the results to the output directory.
    ofile = os.path.join(outdir, f'GPI10.monthly.{year}.nc')
    encoding = {'GPI10': {'dtype': 'float32', 'zlib': True, 'complevel': 1}}
    GPI10.to_dataset(name='GPI10') \
         .to_netcdf(ofile, encoding=encoding, unlimited_dims='time')
    print('The output for GPI (Emanuel 2010) has been saved to:', ofile)


def main():

    year = 2020

    # The directory where you store your ERA5 reanalysis data;
    # the default value (datadir=None) is the current folder
    datadir = '/Users/serene.meng/TCindex/tcicalc/ERA5_Data/'

    # The directory to store the output TC indexes;
    # the default value (outdir=None) is the current folder
    outdir = '/Users/serene.meng/TCindex/tcicalc/TC_example/'


    # calculate and save the TC indexes.
    calc_PI(year=year, datadir=datadir, outdir=outdir)
    calc_chim(year=year, p_b=950e2, p_m=600e2, datadir=datadir, outdir=outdir)
    calc_vws(year=year, datadir=datadir, outdir=outdir)
    calc_VI(year=year, datadir=datadir, outdir=outdir)
    calc_eta(year=year, datadir=datadir, outdir=outdir)
    calc_GPI04(year=year, datadir=datadir, outdir=outdir)
    calc_GPI10(year=year, datadir=datadir, outdir=outdir)
    print('*** The calculations have been completed successfully. ***')

if __name__ == '__main__':
    main()
