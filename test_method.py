# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 08:05:41 2021

@author: Ronan Abhervé
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

#%% Modules

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

#%% Changing path

work_folder = "D:/Users/abherve/GITHUB/StreamNetwork/"

#%% Foxed paths

perso_path = work_folder + "stream-network_beta/"

data_path = perso_path + "/example_data/"
bin_path = perso_path + "/external_bin/"
out_path = work_folder + "/outputs_results/"

#%% Inputs

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

#%% Watershed

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

#%% Dichotomy

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

#%% Plots

results_plot.display_results_map(dem_path, data_path, out_path, site, False)
results_plot.display_results_graph(data_path, out_path, site)

#%% Notes
