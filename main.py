import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature 
from netCDF4 import Dataset
import numpy as np
import metpy.plots as mpplots 
from shapely.geometry import Polygon
from libs.segmentation import * 
from geopy.distance import geodesic 

# Convert pixel distance to actual distance
# Look up geopy calculate distance

# Create a base map to display Watershed and 
def new_map(fig, lon, lat):
  # Create projection centered on the radar. Allows us to use x and y relative to the radar
  proj = ccrs.LambertConformal(central_longitude = lon, central_latitude = lat)

  # New axes with the specified projection
  ax = fig.add_axes([0.02, 0.02, 0.96, 0.96], projection = proj)

  # Add coastlines and states
  ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth = 2)
  ax.add_feature(cfeature.STATES.with_scale('50m'))

  return ax

## ** Main Function where everything happens **
def main():
  # Load watershed
  watershed = gpd.read_file("santa_ana_r_a.geojson")

  # Load NEXRAD data from netcdf4 file
  ncfile = '/media/jntp/D2BC15A1BC1580E1/NCFRs/20170217_18.nc'
  nexdata = Dataset(ncfile, mode = 'r')
  print(nexdata)

  # Get data from netcdf file
  lons = nexdata['Longitude'][:][:]
  lats = nexdata['Latitude'][:][:]
  refs = nexdata['Reflectivity'][20] 
  ref_rows, ref_cols = refs.shape  

  # Find pixel dimensions
  pixel_length, pixel_width = get_pixel_dimensions(np.max(lats), np.max(lons), np.min(lats), \
          np.min(lons), ref_rows, ref_cols) 

  # Specify a central longitude and latitude (i.e. reference point)
  central_lon = -117.636
  central_lat = 33.818

  # Create a new figure and map 
  fig = plt.figure(figsize = (10, 10))
  ax = new_map(fig, central_lon, central_lat) # -117.636, 33.818 

  # Set limits in lat/lon space
  ax.set_extent([-121, -114, 32, 36]) # SoCal

  # Get color table and value mapping info for the NWS Reflectivity data
  ref_norm, ref_cmap = mpplots.ctables.registry.get_with_steps('NWSReflectivity', 5, 5) 

  # Transform to this projection
  use_proj = ccrs.LambertConformal(central_longitude = central_lon, central_latitude = central_lat)

  # Transfer lats, lons matrices from geodetic lat/lon to LambertConformal
  out_xyz = use_proj.transform_points(ccrs.Geodetic(), lons, lats)
  
  # Separate x, y from out_xyz
  x = out_xyz[:, :, 0] 
  y = out_xyz[:, :, 1]

  new_refs = np.empty((1336, 1506)) 
  new_refs[:] = np.nan

  # Find a numpy function that gives you the VALUES based off of specified indices
  indices = np.where(refs >= 5)
  i = indices[0]
  j = indices[1]
  results = refs[refs >= 5]

  # Create new matrix of indices where reflectivity >= 5 that is retrofitted for numpy.put
  full_indices = i * 1506 + j 

  # Take the reflectivity values that meet the standard and put in new_refs matrix
  np.put(new_refs, full_indices, results)

  # Add watershed geometry 
  ax.add_geometries(watershed.geometry, crs = ccrs.PlateCarree(), zorder = 1, facecolor = 'red', edgecolor = 'red')

  # Identify NCFR; use reflectivity >= 45 dbz as the ncfr threshold
  ncfr = np.empty((1336, 1506))
  ncfr[:] = np.nan 
  ncfr_indices = np.where(new_refs >= 45)
  m = ncfr_indices[0] # size 83887
  n = ncfr_indices[1] # size 83887
  ncfr_refs = new_refs[new_refs >= 45]
  ncfr_full_indices = m * 1506 + n
  np.put(ncfr, ncfr_full_indices, ncfr_refs)

  # Obtain the x and y coordinates from the NCFR indices
  x_ncfr = x[m, n] 
  y_ncfr = y[m, n] 

  # Segmentation 
  # conv_cells = find_convective_cells(refs) 
  # narrow_conv_cells = remove_wide_cells(refs, conv_cells)
  # connect_cores(refs, 15) # wrong

  # Add colormesh (radar reflectivity) 
  ax.pcolormesh(x, y, new_refs, cmap = ref_cmap, norm = ref_norm, zorder = 2) 

  # plt.show()

if __name__ == '__main__':
  main() 

# Next step identify NCFR and the intersection
# Clean code later
