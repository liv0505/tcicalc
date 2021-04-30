#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Lingwei Meng (lingweim@princeton.edu)
import cartopy.crs as ccrs
import xarray as xr
import os.path
from cartopy.mpl.geoaxes import GeoAxes
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
import matplotlib.colors as colors



def TCI_plot(da, vmin, vmax, norm, cmap, title, ofig=None):
    projection = ccrs.PlateCarree(central_longitude = 180)
    axes_class = (GeoAxes,
                  dict(map_projection=projection))
    fig = plt.figure(figsize = (15, 12),dpi = 300)
    axgr = AxesGrid(fig, 111, axes_class=axes_class,
                    nrows_ncols=(4, 3),
                    axes_pad=0.4,
                    cbar_mode='None',
                    label_mode='')
    for i, ax in enumerate(axgr):
        dap = np.maximum(da.isel(time = i), vmin)
        dap = np.minimum(da.isel(time = i), vmax)
        scale='110m'
        #ax.stock_img()
        ax.coastlines(scale,linewidth = 0.3)
        ax.set_xticks(np.linspace(-180, 180, 5), crs=projection)
        ax.set_yticks(np.linspace(-90, 90, 5), crs=projection)
        ax.tick_params(axis = 'x', labelsize = 7, colors = 'gray')
        ax.tick_params(axis = 'y', labelsize = 7, colors = 'gray')
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)
        p = dap.plot(ax = ax, transform = ccrs.PlateCarree(),
                  norm=norm,
                  cmap=cmap,
                  add_colorbar = False)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.title.set_fontsize(9)
    cb_ax = fig.add_axes([0.92, 0.2, 0.02, 0.6])
    cbar = fig.colorbar(p, cax=cb_ax)
    cbar.set_label(title, fontsize=10)
    if ofig is not None:
        plt.savefig(ofig, dpi=300)
        print('The figure for ' + title + ' has been saved to:', ofig)


def main():
    outdir = '/Users/serene.meng/TCindex/tcicalc/TC_example/'
    year = 2020

    ifile = os.path.join(outdir, f'VI.monthly.{year:04d}.nc')
    ofile = os.path.join(outdir, f'ventilation_index_{year:04d}.png')
    vi = xr.open_dataarray(ifile)
    TCI_plot(vi, vmin=1e-3, vmax=5,
                 norm=colors.LogNorm(vmin=1e-3, vmax=5),
                 cmap='PuBu_r',
                 title='Ventilation Index',
                 ofig=ofile)

    ifile = os.path.join(outdir, f'PI.monthly.{year:04d}.nc')
    ofile = os.path.join(outdir, f'potential_intensity_{year:04d}.png')
    pi = xr.open_dataset(ifile).vmax
    TCI_plot(pi, vmin=0, vmax=120,
                 norm=colors.Normalize(vmin=0, vmax=120),
                 cmap='RdYlBu_r',
                 title='Potential Index (m/s)',
                 ofig=ofile)

    ifile = os.path.join(outdir, f'GPI04.monthly.{year:04d}.nc')
    ofile = os.path.join(outdir, f'GPI04_{year:04d}.png')
    gpi04 = xr.open_dataarray(ifile)
    TCI_plot(gpi04, vmin=1e-5, vmax=300,
                 norm=colors.LogNorm(vmin=1e-5, vmax=300),
                 cmap= 'RdYlBu_r',
                 title='Genesis Potential Index (Emanuel and Nolan 2004)',
                 ofig=ofile)

    ifile = os.path.join(outdir, f'GPI10.monthly.{year:04d}.nc')
    ofile = os.path.join(outdir, f'GPI10_{year:04d}.png')
    gpi10 = xr.open_dataarray(ifile)
    TCI_plot(gpi10, vmin=1e-17, vmax=1e-14,
                 norm=colors.LogNorm(vmin=1e-17, vmax=1e-14),
                 cmap='Reds',
                 title=r'Genesis Potential Index ($m^{-2}s^{-1}$)(Emanuel 2010)',
                 ofig=ofile)

if __name__ == '__main__':
    main()
