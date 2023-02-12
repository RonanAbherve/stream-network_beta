# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 08:05:41 2021

@author: Ronan Abhervé
"""

#%% Library necessary to import in this script

import os
import pandas as pd
from os.path import dirname, abspath
import sys

#%% Librairies necessary to install for others scripts


import flopy

import glob
import numpy as np
from osgeo import gdal, osr
import shutil

import geopandas as gpd
import imageio
import rasterio
import rasterio as rio
import shapely

import rasterio.plot
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
from matplotlib.font_manager import FontProperties

from decimal import Decimal

#%% Modules necessary to run the script

root_dir = dirname(abspath(__file__))
sys.path.append(root_dir)
DIR = dirname(dirname(dirname(abspath(__file__))))
sys.path.append(DIR)

cwd = os.getcwd()
import whitebox
wbt = whitebox.WhiteboxTools()
wbt.set_verbose_mode(False)

from src_python.watershed_extraction import watershed_process, watershed_topography
from src_python.groundwaterflow_model import modflow_model, post_processing
from src_python.calibration_method import launch_dichotomy, objective_function
from src_python.display_results import toolbox, results_plot

os.chdir(cwd)

#%% Changing the path of the parent folder downloaded

work_folder = "D:/.../StreamNetwork/"
work_folder = "D:/Users/abherve/GITHUB/StreamNetwork/"

#%% Paths fixed from the folder downloaded on the GitHub page

perso_path = work_folder + "stream-network_beta/"

data_path = perso_path + "/example_data/"
bin_path = perso_path + "/external_bin/"
out_path = work_folder + "/outputs_results/"

#%% Example of input from the provided data

outlet_caract = pd.read_csv(perso_path + "/test_outlet.txt", sep='\t', header=0, engine='python')

dem_path = data_path + 'BDALTI_Brittany_75m.tif'
obs = 'streams'

fid = outlet_caract.FID.values[0] 
site = outlet_caract.Site.values[0]
snap = outlet_caract.Snap.values[0]
buffer = outlet_caract.Buffer.values[0]

thick = 30
lay = 1
first_kr = 1
last_kr = 10000
gap = 1
exe = 'mfnwt.exe'
rech = 200 / 12 / 1000 # mm/year to m/months

#%% Catchment extraction from the outlet coordinates

print('### WATERSHED '+site.upper()+' ###')
    
launch_dichotomy.delimit_size(dem_path=dem_path,
                              outlet=outlet_caract,
                              watershed=site,
                              snap_dist=snap,
                              buff_percent=buffer,
                              save_gis=True,
                              type_obs=obs,
                              data_path=data_path,
                              tmp_path=out_path+'_tmp/',
                              out_path=out_path)

#%% Calibration from stream network based on dichotomy approach

print('### DICHOTOMY '+site.upper()+' ###')

launch_dichotomy.dichotomy_loop(first=first_kr,
                                last=last_kr,
                                gap=gap,
                                df=pd.DataFrame(),
                                watershed=site,
                                climatic=[rech],
                                lay_number=lay,
                                thick=thick,
                                type_obs=obs,
                                data_path=data_path,
                                out_path=out_path,
                                exe=bin_path+exe)

#%% Plots of major results

# Map of the best simulation with the K optimal estimated
results_plot.display_results_map(dem_path, data_path, out_path, site, False)

# Graph of Dos and Dso definig the best K/R value
results_plot.display_results_graph(data_path, out_path, site)
