# -*- coding: utf-8 -*-

import os
import sys
import geopandas as gpd
import numpy as np
import imageio
'''os.path.dirname(os.getcwd())'''
sys.path.append(os.getcwd())

### Method 1
import whitebox
wbt = whitebox.WhiteboxTools()
wbt.set_verbose_mode(False)
### Method 2
# from WBT.whitebox_tools import WhiteboxTools
# wbt = WhiteboxTools()

from src_python.watershed_extraction import watershed_process, watershed_topography

class extract_observed:
    def __init__(self, watershed='name', type_obs='streams', data_path=os.path.dirname(os.getcwd())+'\\data\\',
                 out_path=os.path.dirname(os.getcwd())+'\\output\\'):
        
        self.ws = os.getcwd()
        self.watershed = watershed
        self.type_obs = type_obs
        self.data_path = data_path
        self.out_path = out_path
        
        self.out_fold = self.out_path + self.watershed + '/'
        
        self.gis_path = self.out_fold + '/gis/'
        self.obs_path = self.out_fold + '/obs/'
        if not os.path.exists(self.obs_path):
            os.makedirs(self.obs_path)
        
        self.watershed_buff = self.gis_path + 'buff.shp'
        self.watershed_shp = self.gis_path + 'watershed.shp'
        self.watershed_fill = self.gis_path + 'watershed_fill.tif'
        self.watershed_buff_fill = self.gis_path + 'watershed_buff_fill.tif'
        
        self.streams = self.data_path + self.type_obs + '.shp'
        
        self.clip_observed()
        
    def clip_observed(self):
        
        ### WATSHD
        self.clip_streams = self.obs_path + self.type_obs + '.shp'
        wbt.clip(self.streams, self.watershed_shp, self.clip_streams)
        if self.watershed == 'Cheze':
            shp = gpd.read_file(self.clip_streams)
            shp = shp[shp.gid != 119811]
            shp.to_file(self.clip_streams)
        if self.watershed == 'Laizon':
            shp = gpd.read_file(self.clip_streams)
            shp = shp[shp.gid != 82698]
            shp.to_file(self.clip_streams)
        forma = gpd.read_file(self.clip_streams)
        forma = forma.geom_type[0] 
        if (forma != 'MultiPolygon') | (forma != 'Polygon'):
            self.tif_streams = self.obs_path + self.type_obs + '.tif'
            wbt.vector_lines_to_raster(self.clip_streams, self.tif_streams, field="FID", base=self.watershed_fill)
        else:
            self.tif_streams = self.obs_path + self.type_obs + '.tif'
            wbt.vector_polygons_to_raster(self.clip_streams, self.tif_streams, field="FID", base=self.watershed_fill)
        self.pt_streams = self.obs_path + self.type_obs + '_pt.shp'
        wbt.raster_to_vector_points(self.tif_streams, self.pt_streams)
        
        # BUFF
        self.clip_streams_buff = self.obs_path + self.type_obs + '_buff.shp'
        wbt.clip(self.streams, self.watershed_buff, self.clip_streams_buff)
        forma = gpd.read_file(self.clip_streams_buff)
        forma = forma.geom_type[0] 
        if (forma != 'MultiPolygon') | (forma != 'Polygon'):
            self.tif_streams_buff = self.obs_path + self.type_obs + '_buff.tif'
            wbt.vector_lines_to_raster(self.clip_streams_buff, self.tif_streams_buff, field="FID", base=self.watershed_buff_fill)
        else:
            self.tif_streams_buff = self.obs_path + self.type_obs + '_buff.tif'
            wbt.vector_polygons_to_raster(self.clip_streams_buff, self.tif_streams_buff, field="FID", base=self.watershed_buff_fill)
        self.pt_streams_buff = self.obs_path + self.type_obs + '_pt_buff.shp'
        wbt.raster_to_vector_points(self.tif_streams_buff, self.pt_streams_buff)
                
        return self, os.chdir(self.ws)

class generate_distances:
    def __init__(self, watershed='name', type_obs='streams', type_time='s', sim_id='identify', data_path = os.path.dirname(os.getcwd())+'\\data\\',
                 out_path=os.path.dirname(os.getcwd())+'\\output\\'):
        
        self.ws = os.getcwd()
        self.watershed = watershed
        self.type_obs = type_obs
        self.sim_id = sim_id
        self.data_path = data_path
        self.out_path = out_path
        
        self.out_fold = self.out_path + self.watershed + '/'
        
        self.gis_path = self.out_fold + '/gis/'
        self.obs_path = self.out_fold + '/obs/'
        if not os.path.exists(self.obs_path):
            os.makedirs(self.obs_path)
        
        self.sim_fold = self.out_fold + self.sim_id + '/'
        
        self.watershed_shp = self.gis_path + 'watershed.shp'
        self.watershed_fill = self.gis_path + 'watershed_fill.tif'
        self.watershed_buff_fill = self.gis_path + 'watershed_buff_fill.tif'
        self.watershed_direc = self.gis_path + 'watershed_direc.tif'
        
        self.tif_obs = self.obs_path + self.type_obs + '.tif'
        self.pt_obs = self.obs_path + self.type_obs + '_pt.shp'
        
        self.tif_obs_buff = self.obs_path + self.type_obs + '_buff.tif'
        self.pt_obs_buff = self.obs_path + self.type_obs + '_pt_buff.shp'

        self.seep_sim = self.sim_fold + 'seepage.tif'
        self.seep_sim_mask = self.sim_fold + 'mask_seepage.tif'
        self.drn_sim = self.sim_fold + 'outflow.tif'
        self.drn_sim_mask = self.sim_fold + 'mask_outflow.tif'
        self.wt_sim = self.sim_fold + 'watertable.tif'
        self.wt_sim_mask = self.sim_fold + 'mask_watertable.tif'
                
        self.clip_sim()
        self.sim_to_obs()
        self.obs_to_sim()

    def clip_sim(self):
        wbt.clip_raster_to_polygon(self.seep_sim, self.watershed_shp, self.seep_sim_mask)
        wbt.clip_raster_to_polygon(self.drn_sim, self.watershed_shp, self.drn_sim_mask)
        wbt.clip_raster_to_polygon(self.wt_sim, self.watershed_shp, self.wt_sim_mask)
        return self, os.chdir(self.ws)

    def sim_to_obs(self):
        self.dist_sim_obs = self.sim_fold + 'dist_sim_obs.tif'
        wbt.downslope_distance_to_stream(self.watershed_buff_fill, self.tif_obs_buff, self.dist_sim_obs)
        self.sim_shp = self.sim_fold + 'sim.shp'
        wbt.raster_to_vector_points(self.seep_sim_mask, self.sim_shp)
        self.sim_flow = self.sim_fold + 'simflow.tif'
        wbt.trace_downslope_flowpaths(self.sim_shp, self.watershed_direc, self.sim_flow)
        self.pt_sim_flow = self.sim_fold + 'simflow.shp'
        wbt.raster_to_vector_points(self.sim_flow, self.pt_sim_flow)
        wbt.add_point_coordinates_to_table(self.pt_sim_flow)
        wbt.extract_raster_values_at_points(self.dist_sim_obs, self.pt_sim_flow)
        wbt.vector_points_to_raster(self.pt_sim_flow, self.sim_fold + 'simflow_raster.tif',
                                    field="VALUE1", base=self.sim_fold + 'dist_sim_obs.tif')
        wbt.length_of_upstream_channels(self.watershed_direc, self.sim_fold + 'simflow.tif',
                                        self.sim_fold + 'simflow_length.tif')
        wbt.raster_streams_to_vector(self.sim_fold + 'simflow.tif',
                             self.watershed_direc, 
                             self.sim_fold + 'simflow_to_vec.shp')

        return self, os.chdir(self.ws)
                
    def obs_to_sim(self):
        self.dist_obs_sim = self.sim_fold + 'dist_obs_sim.tif'
        wbt.downslope_distance_to_stream(self.watershed_fill, self.sim_flow, self.dist_obs_sim)        
        self.obs_flow = self.sim_fold + 'obsflow.tif'
        wbt.trace_downslope_flowpaths(self.pt_obs, self.watershed_direc, self.obs_flow)
        self.pt_obs_flow = self.sim_fold + 'obsflow.shp'
        wbt.raster_to_vector_points(self.obs_flow, self.pt_obs_flow)
        wbt.add_point_coordinates_to_table(self.pt_obs_flow)
        wbt.extract_raster_values_at_points(self.dist_obs_sim, self.pt_obs_flow)
        wbt.vector_points_to_raster(self.pt_obs_flow, self.sim_fold + 'obsflow_raster.tif',
                            field="VALUE1", base=self.sim_fold + 'dist_obs_sim.tif')
        wbt.length_of_upstream_channels(self.watershed_direc, self.sim_fold + 'obsflow.tif',
                                        self.sim_fold + 'obsflow_length.tif')
        wbt.raster_streams_to_vector(self.sim_fold + 'obsflow.tif',
                             self.watershed_direc, 
                             self.sim_fold + 'obsflow_to_vec.shp')
        return self, os.chdir(self.ws)

class store_dataframe:
    def __init__(self, watershed='name', type_obs='streams', type_time='s', sim_id='identify',
                 out_path=os.path.dirname(os.getcwd())+'\\output\\'):
        self.ws = os.getcwd()
        self.watershed = watershed
        self.type_time = type_time
        self.sim_id = sim_id
        self.out_path = out_path
                
        self.out_fold = self.out_path + self.watershed + '\\'
        
        self.gis_path = self.out_fold + '/gis/'
        self.obs_path = self.out_fold + '/obs/'

        self.sim_fold = self.out_fold + self.sim_id + '\\'

        self.watershed_fill = self.gis_path + 'watershed_fill.tif'
        self.dem = watershed_topography.dem(self.watershed_fill)
        
        self.pt_obs_flow = self.sim_fold + 'obsflow.shp'
        self.pt_sim_flow = self.sim_fold + 'simflow.shp'
        
        self.drn_sim_mask = self.sim_fold + 'mask_outflow.tif'
        
        self.mean_distances()
        self.total_length()
        self.mean_outflow()
        self.fuzzy()
        
    def mean_distances(self):
        self.obs_to_sim = gpd.read_file(self.pt_obs_flow)
        self.obs_to_sim = self.obs_to_sim.rename(columns={'VALUE':'count', 'VALUE1':'distance'})
        self.obs_to_sim = self.obs_to_sim[self.obs_to_sim['distance'] >= 0]
        self.obs_to_sim_mean = np.nanmean(self.obs_to_sim['distance'])
        self.sim_to_obs = gpd.read_file(self.pt_sim_flow)
        self.sim_to_obs = self.sim_to_obs.rename(columns={'VALUE':'count', 'VALUE1':'distance'})
        self.sim_to_obs = self.sim_to_obs[self.sim_to_obs['distance'] >= 0]
        self.sim_to_obs_mean = np.nanmean(self.sim_to_obs['distance'])
        return self, os.chdir(self.ws)
    
    def total_length(self):
        # self.obs_length = imageio.imread(self.sim_fold + 'obsflow_length.tif')
        # self.obs_length = np.nanmax(self.obs_length)
        # self.sim_length = imageio.imread(self.sim_fold + 'simflow_length.tif')
        # self.sim_length = np.nanmax(self.sim_length)
        
        self.obs_length = gpd.read_file(self.sim_fold + 'obsflow_to_vec.shp')
        self.obs_length['length'] = self.obs_length.geometry.length
        self.obs_length = self.obs_length[(self.obs_length.STRM_VAL>2) | (self.obs_length.length>110)]
        self.obs_length.to_file(self.sim_fold + 'obsflow_to_vec.shp')
        self.obs_length = self.obs_length.length.sum()
        
        self.sim_length = gpd.read_file(self.sim_fold + 'simflow_to_vec.shp')
        self.sim_length['length'] = self.sim_length.geometry.length
        self.sim_length = self.sim_length[(self.sim_length.STRM_VAL>2) | (self.sim_length.length>110)]
        self.sim_length.to_file(self.sim_fold + 'simflow_to_vec.shp')
        self.sim_length = self.sim_length.length.sum()
        
    def mean_outflow(self):
        self.flux = imageio.imread(self.drn_sim_mask) # L/T
        self.flux = np.ma.masked_array(self.flux, mask=(self.dem.data==-99999))
        self.cell = self.flux.count()
        self.outflow = (np.nansum(self.flux) / (self.cell * self.dem.pixel**2)) # M/T
        return self, os.chdir(self.ws)
    
    def fuzzy(self):
        maskdata = imageio.imread(self.watershed_fill)
        obs = imageio.imread(self.sim_fold + 'obsflow.tif')
        obs[obs>=0] = 1
        obs[obs<=0] = 0
        obs[maskdata==-99999] = np.nan
        
        raw = gpd.read_file(self.sim_fold + 'sim.shp')
        raw = len(raw)
        
        sim = imageio.imread(self.sim_fold + 'simflow.tif')
        sim[sim>=0] = 2
        sim[sim<=0] = 0
        sim[maskdata==-99999] = np.nan
        
        self.Sraw = raw
        
        self.So = np.count_nonzero(obs==1)
        self.No = np.count_nonzero(obs==0)
        self.Sm = np.count_nonzero(sim==2)
        self.Nm = np.count_nonzero(sim==0)
        
        self.diff = obs + sim
        self.Si = np.count_nonzero(self.diff==2) # incorrect pred (saturated but not)
        self.Ni = np.count_nonzero(self.diff==1) # incorrect pred (unsaturated but not)
        self.Sc = np.count_nonzero(self.diff==3) # correct pred (saturated)
        self.Nc = np.count_nonzero(self.diff==0) # correct pred (unsaturated)
        
        self.Ea = (1 - ((np.abs(self.Sm - self.So))/self.So))
        self.Sa = 1 - (self.Si / self.So)
        self.Na = 1 - (self.Ni  / self.No)
        self.E = self.Ea * (self.Sa) * (self.Na)
        
        return self, os.chdir(self.ws)
