import geopandas
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature 

# Create Map
def new_map(fig, lon, lat, gdf):
  # Create projection centered on the radar. Allows us to use x and y relative to the radar
  proj = ccrs.LambertConformal(central_longitude = lon, central_latitude = lat)
  
  # Convert into a 'proj4' string/dict compatible with GeoPandas
  crs_proj4 = proj.proj4_init
  df_ae = gdf.to_crs(crs_proj4) 

  # Add new geometries here
  # new_geometries = [proj.project_geometry(ii, src_crs = crs) for ii in df_ae['geometry'].values]

  # fig, ax = plt.subplots(subplot_kw = {'projection': proj}

  # New axes with the specified projection
  ax = fig.add_axes([0.02, 0.02, 0.96, 0.96], projection = proj)

  # Add coastlines and states
  ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth = 2)
  ax.add_feature(cfeature.STATES.with_scale('50m'))

  ax.add_geometries(df_ae['geometry'], crs = proj)

  return ax

def main():
  # Load watershed
  watershed = geopandas.read_file("santa_ana_r_a.geojson")

  # Create new map based off of CRS Projection
  proj = ccrs.LambertConformal(central_longitude = -117.636, central_latitude = 33.818)

  # Create a new figure
  fig = plt.figure(figsize = (10, 10))

  # New axes with the specified projection
  ax = fig.add_axes([0.02, 0.02, 0.96, 0.96], projection = proj)

  # Set limits in lat/lon space
  ax.set_extent([-121, -114, 32, 36])

  # Add coastlines andd states
  ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth = 2)
  ax.add_feature(cfeature.STATES.with_scale('50m'))

  # Add watershed geometry
  ax.add_geometries(watershed.geometry, crs = ccrs.PlateCarree(), zorder = 5, edgecolor = 'blue')

  plt.show()

if __name__ == '__main__':
  main()

# Next step... Clean up code
# Also next big step would be figuring out to incorporate NCFR core as "polygons"
