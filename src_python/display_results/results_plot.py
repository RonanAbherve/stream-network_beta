# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 09:16:58 2022

@author: ronan
"""

#%% Librairies

import fiona
import os
import pandas as pd
import numpy as np
from glob import glob
import geopandas as gpd
import imageio
from osgeo import gdal
import rasterio
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from matplotlib.font_manager import FontProperties
import shapely
shapely.speedups.disable()
from os.path import dirname, abspath
import sys
from decimal import Decimal
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from matplotlib_scalebar.scalebar import ScaleBar
import os
from glob import glob
import geopandas as gpd
from osgeo import gdal
import rasterio
from mpl_toolkits.axes_grid1 import make_axes_locatable

root_dir = dirname(abspath(__file__))
sys.path.append(root_dir)
DIR = dirname(dirname(dirname(abspath(__file__))))
sys.path.append(DIR)

import whitebox
wbt = whitebox.WhiteboxTools()
wbt.set_verbose_mode(False)
from src_python.display_results import toolbox

#%% Map

def display_results_map(dem_path, data_path, out_path, site, colorbar):
    
    fontprop = toolbox.plot_params(8,15,18,20) # small, medium, interm, large

    geol_s = gpd.read_file(data_path+'GEO1M.shp')
    obs = 'streams'
            
    fig, axs = plt.subplots(1, 2, figsize=(5,5), dpi=300)
    axs = axs.ravel()
    
    watershed_name = site
    obs_folder = out_path+'/'+watershed_name+'/obs/'
    gis_folder = out_path+'/'+watershed_name+'/gis/'
    simulations_folder = out_path+'/'+watershed_name+'/'
    
    df = pd.read_csv(simulations_folder+watershed_name+'_dichotomy.csv', sep='\t', header=0)
    kroptim = df.iloc[-1]['Kr'].round(3)
    koptim = df.iloc[-1]['K'].round(3)
    doptim =  int(((df.iloc[-1]['Oflow'] + df.iloc[-1]['Sflow'].round(3))/2).round(0))
    scan = sorted(glob(simulations_folder+'/'+'dic*'), key=os.path.getmtime)
    for ids, j in enumerate(scan):
        split = j.split('\\')[-1].split('_')[5]
        if split==str(kroptim):
            optimcase = j
            
    streams = gpd.read_file(obs_folder+obs+'.shp')
    polyg = gpd.read_file(gis_folder+'watershed.shp')
    contour = gpd.read_file(gis_folder+'watershed_contour.shp')
    bounds = contour.geometry.total_bounds
    xlim = ([bounds[0], bounds[2]])
    ylim = ([bounds[1], bounds[3]])   
    dem = gis_folder+'watershed_extent.tif'
    gdal.Translate(dem, gdal.Open(dem_path), projWin=[xlim[0],ylim[1],xlim[1],ylim[0]], noData=-99999)
    hill = gis_folder+'watershed_extent_hill.tif'
    wbt.hillshade(dem, hill, azimuth=315.0, altitude=45.0, zfactor=2)    
    dem = rasterio.open(gis_folder+'watershed_extent.tif')
    hill = rasterio.open(gis_folder+'watershed_extent_hill.tif')
    img = imageio.imread(gis_folder+'watershed_extent.tif')
    
    simflow = gpd.read_file(optimcase+'/simflow.shp')
    raster = rasterio.open(optimcase+'/simflow.tif')
    
    ax=axs[0] 
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.set_title(watershed_name.upper(), fontproperties=fontprop)
    ax.set(aspect='equal')
    scalebar = AnchoredSizeBar(ax.transData, 2000, '2 km', 'lower left', 
                               pad=0.2, color='k', frameon=False, size_vertical=1)
    ax.add_artist(scalebar)
    
    image_hidden = ax.imshow(np.ma.masked_where(dem.read(1) < 0, dem.read(1)), cmap='terrain')
    mnt = rasterio.plot.show(np.ma.masked_where(dem.read(1) < 0, dem.read(1)), ax=ax, transform=dem.transform, cmap='terrain', alpha=1, zorder=2, aspect="auto")
    hil = rasterio.plot.show(np.ma.masked_where(hill.read(1) < 0, hill.read(1)), ax=ax, transform=dem.transform, cmap='Greys_r', alpha=0.5, zorder=2, aspect="auto")
    streams.plot(ax=ax, lw=1, color='navy', zorder=3, edgecolor='none')
    contour.plot(ax=ax, lw=1.5, color='k', zorder=6)
    
    if colorbar == True:
        divider = make_axes_locatable(ax)
        cax = divider.new_vertical(size="2.5%", pad=0.05, pack_start=True)
        fig.add_axes(cax)
        cbar = fig.colorbar(image_hidden, cax=cax, orientation="horizontal")
        ticklabels = cbar.ax.get_ymajorticklabels()
        ticks = list(cbar.get_ticks())
        val = np.ma.masked_where(dem.read(1) < 0, dem.read(1))
        minVal =  int(round(np.min(val[np.nonzero(val)],0)))
        maxVal =  int(round(np.max(val[np.nonzero(val)],0)))
        meanVal = int(round(minVal+((maxVal-minVal)/2),0))
        cbar.set_ticks([minVal, meanVal, maxVal])
        cbar.set_ticklabels([minVal, meanVal, maxVal])
        cbar.mappable.set_clim(minVal, maxVal)
        cbar.ax.tick_params(labelsize=10)    
        cbar.ax.yaxis.set_ticks_position('left')
        cbar.ax.tick_params(size=0)
    
    ax=axs[1]
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.set_title('K = '+('%.1E'%Decimal(koptim/24/3600))+' m/s'+'  -  '+'D = '+str(doptim)+' m',
                 fontproperties=fontprop)
    xlims = ax.get_xlim()[1] - ax.get_xlim()[0]
    ylims = ax.get_ylim()[1] - ax.get_ylim()[0]
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    width *= fig.dpi
    height *= fig.dpi
    
    hil = rasterio.plot.show(np.ma.masked_where(hill.read(1) < 0, hill.read(1)), ax=ax, transform=dem.transform, cmap='Greys_r', alpha=0.5, zorder=2)
    geol_s.plot(ax=ax, color=list(geol_s['hex']),alpha=0.3, edgecolor='dimgrey', zorder=0) 
    streams.plot(ax=ax, lw=1, color='navy', zorder=3, edgecolor='none')
    contour.plot(ax=ax, lw=1.5, color='k', zorder=5)    
    simflow.plot(ax=ax, alpha=1, column='VALUE1', cmap="RdYlGn_r", 
                  marker='s', markersize=5, lw=0.1, edgecolor='none', scheme="User_Defined", 
                  classification_kwds=dict(bins=[150, 450, 750]), zorder=4)
    
    dem = None
    
    plt.tight_layout()
    fig.tight_layout()
    fig.savefig(out_path+watershed_name+'/Mapping_results'+'.png', dpi=300, bbox_inches='tight', transparent=False)


#%% Graph

def display_results_graph(data_path, out_path, site):
    
    fontprop = toolbox.plot_params(8,15,18,20) # small, medium, interm, large
    
    watershed_name = site
    obs_folder = out_path+'/'+watershed_name+'/obs/'
    gis_folder = out_path+'/'+watershed_name+'/gis/'
    simulations_folder = out_path+'/'+watershed_name+'/'
    
    df = pd.read_csv(simulations_folder+watershed_name+'_dichotomy.csv', sep='\t', header=0)
    
    fig, axs = plt.subplots(1,1, figsize=(4.5,4), sharex=True, sharey=True)
    (ax1) = axs
    toplot = df.sort_values('Kr')
    toplot=toplot.mask(toplot==0).fillna(0.1)
    toplot['Ratio'] = toplot.Sflow / toplot.Oflow
    
    n = len(df)
    col1 = pl.cm.RdYlGn(np.linspace(0,1,n))
    col2 = pl.cm.RdYlBu_r(np.linspace(0,1,n))
        
    def color(s, classif):
        lower = classif[0]
        middle = classif[1]
        upper = classif[2]
        s25 = np.ma.masked_where(s > lower, s)
        s50 = np.ma.masked_where((s < lower) | (s > middle), s)
        s75 = np.ma.masked_where((s < middle) | (s > upper), s)
        s100 = np.ma.masked_where(s < upper, s)
        return s25, s50, s75, s100
    test = color(toplot.Oflow, [75, 500, 1000])
    
    xcross = df.iloc[-1]['Kr']
    ycross = (df.iloc[-1]['Sflow'] + df.iloc[-1]['Oflow']) / 2
    
    for t in range(len(df)):
    
        ax = ax1
        # ax.plot(toplot.Kr.iloc[t], toplot.Sflow.iloc[t], 'o', ms=3, lw=0.1, zorder=3, c=col1[t], mfc='k', 
        #         mec='k', mew=1.5)
        # ax.plot(toplot.Kr[t], toplot.Oflow[t], 'o', ms=3, lw=0.1, zorder=3,  c=col2[t], mfc='gray', 
        #         mec='gray', mew=1.5)
        ax.plot(toplot.Kr, toplot.Oflow, c='gray', ls='-', lw=2, zorder=1, label='$D_{os}$')
        ax.plot(toplot.Kr, toplot.Sflow, c='k', ls='-', lw=2, zorder=1, label='$D_{so}$')
        
        # ax.set_xlim(1000,10000)
        ax.set_xscale('log')
        ax.set_xlabel('K/R [-]')
        ax.set_ylabel('Average distance [m]')
        if t == 0:
            ax.legend(loc='lower center', frameon=True)
        
        y1 = np.array(df.Sflow)
        y2 = np.array(df.Oflow)
        ax.axhline(y=ycross, ls='--', lw=2, color='dodgerblue', zorder=0, label='Distance Ɛ')
        ax.axvline(x=xcross, ls='--', lw=2, color='darkmagenta', zorder=0, label='Optimal K/R')
        # ax.scatter(xcross, ycross, s=100, marker='+', color='red', lw=2, zorder=5)
        if t == 0:
            ax.legend(loc='lower center', frameon=True)
    
    plt.tight_layout()
    fig.savefig(out_path+watershed_name+'/Graphical_results'+'.png', dpi=300, bbox_inches='tight', transparent=False)