# coding:utf-8

import os

import geopandas as gpd
from osgeo import gdal, osr
import numpy as np

### Method 1
import whitebox
wbt = whitebox.WhiteboxTools()
wbt.set_verbose_mode(False)

class extract_watershed:
	def __init__(self, dem_path, outlet, snap_dist=150, buff_percent=10,
                save_gis=True,
                tmp_path=os.path.dirname(os.getcwd())+'\\tmp\\', 
                out_path=os.path.dirname(os.getcwd())+'\\output\\'):
    
		self.ws = os.getcwd()
		self.dem_path = dem_path 
		self.out_path = out_path
		self.tmp_path = tmp_path
		self.save_gis = save_gis
		
		self.outlet = outlet
		self.snap_dist = snap_dist
		self.buff_percent = buff_percent
		
		self.fill = self.tmp_path + 'fill.tif'
		self.direc = self.tmp_path + 'direc.tif'
		self.acc = self.tmp_path + 'acc.tif'
     
		if self.save_gis == True :
			self.gis_path = self.out_path + self.outlet.Site.values[0] + '/gis/'
		else:
			self.gis_path = self.tmp_path

		self.outlet_shp = self.gis_path + 'outlet.shp'
		self.outlet_snap_shp = self.gis_path + 'outlet_snap.shp'
		self.watershed = self.gis_path + 'watershed.tif'
		self.watershed_shp = self.gis_path + 'watershed.shp'
		self.watershed_contour_shp = self.gis_path + 'watershed_contour.shp'		
		self.watershed_dem = self.gis_path + 'watershed_dem.tif'
		self.watershed_fill = self.gis_path + 'watershed_fill.tif'
		self.watershed_direc = self.gis_path + 'watershed_direc.tif'
		self.buffer = self.gis_path + 'buff.shp'
		self.watershed_buff_dem = self.gis_path + 'watershed_buff_dem.tif'
		self.watershed_buff_fill = self.gis_path + 'watershed_buff_fill.tif'
		self.watershed_buff_direc = self.gis_path + 'watershed_buff_direc.tif'
        
		self.generate_watershed_dem()

	def generate_watershed_dem(self):
        
		if not os.path.exists(self.gis_path):
				os.makedirs(self.gis_path)
                
		if not os.path.exists(self.tmp_path):
  				os.makedirs(self.tmp_path)
        
		self.dem = gdal.Open(self.dem_path)
		geodata = self.dem.GetGeoTransform()
		data = self.dem.GetGeoTransform()
		proj = osr.SpatialReference(wkt=self.dem.GetProjection())
		self.crs = 'EPSG:'+str(proj.GetAttrValue('AUTHORITY',1))
		wbt.fill_depressions(self.dem_path, self.fill)
		wbt.d8_pointer(self.fill, self.direc, esri_pntr=False)
		wbt.d8_flow_accumulation(self.fill, self.acc, log=True)
		df = self.outlet
		df.columns = ['FID','Site','X','Y','Snap','Buffer']
		gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['X'], df['Y']), crs=self.crs)
		gdf.to_file(self.outlet_shp)
		wbt.snap_pour_points(self.outlet_shp, self.acc, self.outlet_snap_shp, self.snap_dist)
		wbt.watershed(self.direc, self.outlet_snap_shp, self.watershed, esri_pntr=False)
		wbt.raster_to_vector_polygons(self.watershed, self.watershed_shp)
		wbt.polygons_to_lines(self.watershed_shp, self.watershed_contour_shp)
        
		wbt.polygon_area(self.watershed_shp)
		area = gpd.read_file(self.watershed_shp).AREA[0]/1000000
		area = np.abs(area)
		buff_raw = (area) * (self.buff_percent/100) * 1000
		buff_raw = int(round(buff_raw))
		dist = np.linspace(0,buff_raw,buff_raw+1)*np.abs(geodata[1])
		self.buff_dist = dist[np.abs(dist-buff_raw).argmin()]
        
		site_polyg = gpd.read_file(self.watershed_shp)
		site_polyg.to_file(self.watershed_shp)
		site_polyg['geometry'] = site_polyg.geometry.buffer(self.buff_dist)
		site_polyg.to_file(self.buffer)
        
		wbt.clip_raster_to_polygon(self.dem_path,self.buffer,self.watershed_buff_dem)
		wbt.clip_raster_to_polygon(self.fill,self.buffer,self.watershed_buff_fill)
		wbt.clip_raster_to_polygon(self.direc,self.buffer,self.watershed_buff_direc)
        
		wbt.clip_raster_to_polygon(self.dem_path,self.watershed_shp,self.watershed_dem)
		wbt.clip_raster_to_polygon(self.fill,self.watershed_shp,self.watershed_fill)
		wbt.clip_raster_to_polygon(self.direc,self.watershed_shp,self.watershed_direc)
        
		return self, os.chdir(self.ws)
	