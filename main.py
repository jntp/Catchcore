import geopandas
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature 
from netCDF4 import Dataset
import numpy as np
import metpy.plots as mpplots 

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
  watershed = geopandas.read_file("santa_ana_r_a.geojson")

  # Load NEXRAD data from netcdf4 file
  ncfile = '/media/jntp/D2BC15A1BC1580E1/NCFRs/20170217_18.nc'
  nexdata = Dataset(ncfile, mode = 'r')
  print(nexdata)

  # Get data from netcdf file
  lons = nexdata['Longitude'][:][:]
  lats = nexdata['Latitude'][:][:]
  refs = nexdata['Reflectivity'][0]

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

  # To find coordinate, take i * 1506 + j 

  # How to make an array of the results via incremental i and j (no for loops please...)?

  testmatrix.put([testindex1, testindex2], [24, -56])

  # Find a more efficient way for this? this is too slow***
  # for i, column in enumerate(refs):
    # for j, value in enumerate(column):
      # if refs[i][j] >= 5:
        # test[i][j] = refs[i][j]

  # Add watershed geometry 
  ax.add_geometries(watershed.geometry, crs = ccrs.PlateCarree(), zorder = 1, facecolor = 'red', edgecolor = 'red')

  # Add colormesh (radar reflectivity) 
  ax.pcolormesh(x, y, test, cmap = ref_cmap, norm = ref_norm, zorder = 2) 
 
  # plt.show()

if __name__ == '__main__':
  main() 

# Address Github Deprecation Notice!!!
# First clean, optimize, and annotate code
# Next step identify NCFR and the intersection
