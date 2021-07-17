import geopandas
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime, timedelta
from siphon.radarserver import RadarServer
import metpy.plots as mpplots
import numpy as np

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

# Create an instance of RadarServer and query data from the NEXRAD server based on given time and GPS coordinates
def nexrad_query(lon, lat, year, month, day, hour):
  rs = RadarServer('http://tds-nexrad.scigw.unidata.ucar.edu/thredds/radarServer/nexrad/level2/S3/')
  query = rs.query()
  dt = datetime(year, month, day, month) # time in UTC
  query.lonlat_point(lon, lat).time_range(dt, dt + timedelta(hours = 1))
  catalog = rs.get_catalog(query)
  ds = catalog.datasets[0]
  data = ds.remote_access()

  return catalog, data

# Convert data to float values
def raw_to_masked_float(var, data):
  # Convert unsigned values to range [-127, 128] to [0, 255]
  if var._Unsigned:
    data = data & 255

    # Mask missing points
    data = np.ma.array(data, mask = data == 0)

    # Convert to float using the scale and offset 
    return data * var.scale_factor + var.add_offset 

# Convert polar coordinates to Cartesian (x, y)
def polar_to_cartesian(az, rng):
  # x = r*cos*theta, y = r*sin*theta where theta is in radians
  az_rad = np.deg2rad(az)[:, None]
  x = rng * np.sin(az_rad) 
  y = rng * np.cos(az_rad)

  return x, y

# Create a list that contains a color mesh for each time stamp
def create_mesh(catalog, ref_cmap, ref_norm):
  # Add for loop later
  # mesh = []
  # Pull dataset object out of each list of item
  data = catalog.datasets[0].remote_access() # change later 

  # Pull data as well as variables for azimuth and range
  sweep = 0
  ref_var = data.variables['Reflectivity_HI']
  rng = data.variables['distanceR_HI'][:]
  az = data.variables['azimuthR_HI'][sweep]

  # Convert raw data to floating point values and polar coordinates to Cartesian
  ref = raw_to_masked_float(ref_var, ref_var[sweep])
  x, y = polar_to_cartesian(az, rng)

  # Plot the data
  # mesh = ax.pcolormesh(x, y, ref, cmap = ref_cmap, norm = ref_norm, zorder = 0)
  # Add text later

  return ref, x, y # later should be return mesh  

# Still working on this lol
def identifyNCFR(ref):
  # Create array of zeros (boolean) which will be used to locate where the NCFR is on x, y, ref matrices
  locator = np.zeros((720, 1832), dtype = bool)  

  for i, column in enumerate(ref):
    for j, value in enumerate(column):
      if value > 45:
        locator[i][j] = True

  return locator 

## ** Main Function where everything happens **
def main():
  # Load watershed
  watershed = geopandas.read_file("santa_ana_r_a.geojson")

  # Query data from NEXRAD server
  catalog, data = nexrad_query(-117.636, 33.818, 2017, 2, 18, 1)

  # Create a new figure and map
  fig = plt.figure(figsize = (10, 10))
  ax = new_map(fig, -117.636, 33.818)

  # Set limits in lat/lon space
  ax.set_extent([-121, -114, 32, 36]) # SoCal 

  # Get color table and value mapping info for the NWS Reflectivity data 
  ref_norm, ref_cmap = mpplots.ctables.registry.get_with_steps('NWSReflectivity', 5, 5)

  # Add watershed geometry
  ax.add_geometries(watershed.geometry, crs = ccrs.PlateCarree(), zorder = 1, facecolor = 'red', edgecolor = 'red')

  # Add colormesh (radar reflectivity)
  ref, x, y = create_mesh(catalog, ref_cmap, ref_norm) 
  ax.pcolormesh(x, y, ref, cmap = ref_cmap, norm = ref_norm, zorder = 2) # call the mesh function later

  # Test
  # locator = identifyNCFR(ref)
  # gdf['boundary'] = watershed.boundary
  test = watershed.centroid[0]
  print(test)

  plt.show()

if __name__ == '__main__':
  main()

# Hard part.... how to acquire data from geodataframe? Then how should we identify the "intersection?"
# 1st way is to create polygon out of locator?
# 2nd way is to go the other direction... extrapolate "Cartesian" data from watershed polygon
# 3rd way is to experiment with using centroids, boundaries, etc. as possible metrics
# Will have to find a way to integrate reflectivity data of the different sitees
