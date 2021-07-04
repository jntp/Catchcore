import geopandas
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature 

# Create Map
def new_map(fig, lon, lat):
  # Create projection centered on the radar. Allows us to use x and y relative to the radar
  proj = ccrs.LambertConformal(central_longitude = lon, central_latitude = lat)

  # New axes with the specified projection
  ax = fig.add_axes([0.02, 0.02, 0.96, 0.96], projection = proj)

  # Add coastlines and states
  ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth = 2)
  ax.add_feature(cfeature.STATES.with_scale('50m'))

  return ax

gdf = geopandas.read_file("santa_ana_r_a.geojson")
gdf.to_crs(epsg = 102009) # North America Lambert Conformal Conic

gdf["area"] = gdf.area

gdf.plot("area", legend = True)
plt.show() 

# Test phase
# You left off at meausring distance/making maps of geopandas introduction
# Now overlay watershed onto world map... may need to change projection? 
