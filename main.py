import math 
import geopandas as gpd
import matplotlib
import matplotlib.pyplot as plt  
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature 
from netCDF4 import Dataset
import numpy as np
import metpy.plots as mpplots
import pprint # test
from shapely.geometry import Point 
from libs.segmentation import *
from libs.intersection import *
from skimage.measure import find_contours # temp
from shapely.geometry import shape # temp
from shapely.geometry import Polygon 
from shapely.geometry.polygon import LinearRing
from shapely.ops import polygonize_full # test 

# Create a base map to display Watershed and radar imagery
def new_map(fig, lon, lat):
  # Create projection centered on the radar. Allows us to use x and y relative to the radar
  proj = ccrs.LambertConformal(central_longitude = lon, central_latitude = lat)

  # New axes with the specified projection
  ax = fig.add_axes([0.02, 0.02, 0.96, 0.96], projection = proj)

  # Add coastlines and states
  ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth = 2)
  ax.add_feature(cfeature.STATES.with_scale('50m'))

  return ax

# Creates a new reflectivity matrix that only displays dbZ > 5 and creates transparent background
def new_reflectivity(refs):
  # Create empty matrix, same size as ref
  new_refs = np.empty((1336, 1506)) 
  new_refs[:] = np.nan # Set all equal to nan initially, creates "transparent" background 

  # Find a numpy function that gives you the VALUES based off of specified indices
  indices = np.where(refs >= 5)
  i = indices[0]
  j = indices[1]
  results = refs[refs >= 5]

  # Create new matrix of indices where reflectivity >= 5 that is retrofitted for numpy.put
  full_indices = i * 1506 + j 

  # Take the reflectivity values that meet the standard and put in new_refs matrix
  np.put(new_refs, full_indices, results)

  return new_refs

# Runs the functions needed to delineate the boundaries of the NCFR and its cores
def segmentation(refs, core_buffer = 30, conv_buffer = 3):
  """
    core_buffer - size of search radius to connect cells (default: 30px) 
    conv_buffer - size of "holes" to fill in convective cells and clusters (default: 3px) 
  """
  # Segmentation Steps
  # Step 1 - Extract convective cells
  # Step 2 - Close holes in convective cells
  # Step 3 - Remove wide cells
  # Step 4 - Remove cells horizontally adjacent from suspected NCFR core 
  # Step 5 - Connect cells within a certain search radius
  # Step 6 - Check if the length of connected cells fail to meet NCFR criteria, remove
  # Step 7 - Extract NCFR cores from labeled NCFR given convective criteria
  conv_cells = find_convective_cells(refs) 
  closed_cells = close_holes(conv_cells, conv_buffer)
  narrow_conv_cells = remove_wide_cells(refs, closed_cells)
  narrow_conv_cells2 = remove_adjacent_cells(refs, narrow_conv_cells)
  merged_cells = connect_cells(narrow_conv_cells2, core_buffer) 
  labeled_ncfr = check_axis(refs, merged_cells)  
  labeled_cores = extract_cores(refs, labeled_ncfr, conv_buffer)

  return labeled_ncfr, labeled_cores 

# Plots a single image
def plot_single(ax, x, y, ref_cmap, ref_norm, new_refs, labeled_image):
  # Plot the NCFR "slices"
  ax.contour(x, y, 1 * (labeled_image > 0), colors = ['k',], linewidths = .5, linestyles = 'solid', \
      zorder = 5) 
 
  # Add colormesh (radar reflectivity) 
  ax.pcolormesh(x, y, new_refs, cmap = ref_cmap, norm = ref_norm, zorder = 2)

  plt.show() 


## ** Main Function where everything happens **
def main():
  # Load watershed
  watershed = gpd.read_file("santa_ana_r_a.geojson")
  print(watershed.crs)
  # Left off at looking at projected vs. geographic coordinate system
  # watershed = watershed.to_crs("+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m no_defs")
  # print(watershed.crs)

  # Load NEXRAD data from netcdf4 file
  ncfile = '/media/jntp/D2BC15A1BC1580E1/NCFRs/20170217_18.nc'
  nexdata = Dataset(ncfile, mode = 'r')
  print(nexdata)

  # Get data from netcdf file
  lons = nexdata['Longitude'][:][:]
  lats = nexdata['Latitude'][:][:]
 
  ref_refs = nexdata['Reflectivity'][:] # reference reflectivity
  years = nexdata['Year'][:]
  months = nexdata['Month'][:]
  days = nexdata['Day'][:]
  hours = nexdata['Hour'][:]
  minutes = nexdata['Minute'][:] 
 
  # Find pixel dimensions
  ref_rows, ref_cols = ref_refs[0].shape 
  pixel_length, pixel_width = get_pixel_dimensions(np.max(lats), np.max(lons), np.min(lats), \
          np.min(lons), ref_rows, ref_cols) 

  # Specify a central longitude and latitude (i.e. reference point)
  central_lon = -117.636
  central_lat = 33.818

  # Create a new figure and map 
  fig = plt.figure(figsize = (10, 10))
  ax = new_map(fig, central_lon, central_lat) # -117.636, 33.818 

  # Set limits in lat/lon 
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

  # Test - Extract points of watershed and convert to LambertConformal x, y coordinates
  test_lons = watershed.geometry
  watershed_lon, watershed_lat = watershed.geometry[1].exterior.coords.xy
 
  watershed_lon = np.array(watershed_lon) 
  watershed_lat = np.array(watershed_lat)

  water_xyz = use_proj.transform_points(ccrs.Geodetic(), watershed_lon, watershed_lat)

  x_water = water_xyz[:, 0]
  y_water = water_xyz[:, 1]

  water_points = []
  water_points_test = []
  
  for i, x_val in enumerate(x_water):
    water_points.append((x_water[i], y_water[i]))
    water_points_test.append((watershed_lon[i], watershed_lat[i])) 
 
  # Make linear ring of watershed; Test, reorganize code later 
  watershed_rings = []
  watershed_rings_test = []
  watershed_ring = LinearRing(water_points)
  watershed_ring_test = LinearRing(water_points_test) 
  
  watershed_rings.append(watershed_ring)
  watershed_rings.append(Point(0, 0))
  watershed_rings_test.append(Point(0, 0))

  # Add watershed geometry 
  ax.add_geometries(watershed_rings_test, crs = ccrs.PlateCarree(), zorder = 1, facecolor = 'red', edgecolor = 'red') 

  # Test... finding the linewidth of the segmented contours
  ref_ref = ref_refs[34] 
  labeled_ncfr, labeled_cores = segmentation(ref_ref)
  core_centroids = extract_core_centroids(labeled_cores, ref_ref)
  test_contours = find_contours(labeled_cores, 0)
   
  # It works!
  # Now convert x, y to lat, lon
  geom_proj = ccrs.PlateCarree() # previously in PlateCarree
  shapely_contours = []
  shapely_kontours = [] # test
  shapely_kantours = [] # test

  for contour in test_contours:
    i = np.array(contour[:, 0], dtype = int)
    j = np.array(contour[:, 1], dtype = int) 

    x_contour = x[i, j]
    y_contour = y[i, j]
    
    # Test
    lons_contour = lons[i, j]
    lats_contour = lats[i, j]  

    # Check if this conversion is even necessary; not necessary! Use Lambert Conformal for both systems 
    out_latlon = geom_proj.transform_points(use_proj, x_contour, y_contour) 
  
    contour_lon = out_latlon[:, 0]
    contour_lat = out_latlon[:, 1] 
 
    contour_points = []
    contour_xy = []
    kantour_points = [] # test

    for i, lon in enumerate(contour_lon):
      contour_points.append((contour_lon[i], contour_lat[i]))
      contour_xy.append((x_contour[i], y_contour[i]))
      kantour_points.append((lons_contour[i], lats_contour[i])) # test

    # print(kantour_points)
    # Make the shapely geometry of the labelled core
    shapely_contour = LinearRing(contour_points)
    shapely_contours.append(shapely_contour)
    
    # Test 
    shapely_kontour = LinearRing(contour_xy)
    shapely_kontours.append(shapely_kontour)
    shapely_kantour = Polygon(kantour_points)
    shapely_kantours.append(shapely_kantour) 

    # Plot contour_lon and contour_lat to check for accuracy
    # ax.plot(contour_lon, contour_lat, linewidth = 2) 

    # ax.plot(x_contour, y_contour, linewidth = 2)

  # Plot the shapely contours????
  ax.add_geometries(shapely_kantours, crs = ccrs.PlateCarree(), zorder = 2, facecolor = 'green', edgecolor = 'green')

  # Now do the intersection
  # Why isn't this working...
  # Coordinates of linearring don't seem to agree???

  # geom_contours = [shape(feat["geometry"]) for feat in shapely_contours]
  # print(geom_contours) 

  testing = watershed.iloc[1] 
  testing2 = testing['geometry']
  # testing22 = testing2.exterior # use this as shapely geometry for intersection
  # print(testing2.exterior)
  xpol, ypol = testing2.exterior.coords.xy

  # Test
  watershed_poly_coords = []
  for i, x_val in enumerate(xpol):
    watershed_poly_coords.append((xpol[i], ypol[i]))
  watershed_poly = Polygon(watershed_poly_coords) # Note this produces repeating coordinates
  # print(watershed_poly)

  # Create geoseries (test) 
  watershed_geoms = []
  for watershed_geom in watershed.geometry:
    watershed_geoms.append(watershed_geom) 
  watershed_gs = gpd.GeoSeries(watershed_geoms)

  cores_gs = gpd.GeoSeries(shapely_kantours)
  # cores_gs = cores_gs.set_crs(4326)
  # cores_gs = cores_gs.to_crs(4326) 

  # Test; use this for intersection!!!
  test_boundary = watershed.geometry[1].boundary
  test_series = gpd.GeoSeries(test_boundary)
  ax.add_geometries(test_boundary, crs = ccrs.PlateCarree(), zorder = 1, facecolor = 'red', edgecolor = 'red')
  polygons, dangles, cuts, invalids = polygonize_full(test_boundary) # Left off here... find area of polygons now

  test_intersections = cores_gs.intersection(polygons)
  print(test_intersections) # It returned something??? Plot this 

  ax.add_geometries(test_intersections, crs = ccrs.PlateCarree(), zorder = 3, facecolor = 'yellow', edgecolor = 'yellow')
  plt.show() 

  # new_refs = new_reflectivity(ref_ref)
 

  # plot_single(ax, x, y, ref_cmap, ref_norm, new_refs, labeled_cores)

  ## Animate the Plot
  def animate_contour(i):
    """
      Animates the segmented NCFR radar plot. Called using Matplotlib's FuncAnimation.
      Returns a contour plot with the mesh grid (radar reflectivity) and text attached as artists.
    """
    # Clear axes to prevent contours from clobbering up
    ax.clear()  

    ## Initiation
    # Create transparent background for reflectivity >= 5 dbZ
    new_refs = new_reflectivity(ref_refs[i])

    # Run segmentation algorithm for specific time frame 
    labeled_ncfr, labeled_cores = segmentation(ref_refs[i])

    # Extract "singular" time from time variables
    year = str(int(years[i]))
    month = str(int(months[i]))
    day = str(int(days[i]))
    hour = int(hours[i])
    minute = int(minutes[i])

    # Format the hour and minute so it will display properly
    if hour < 10:
      hour = "0" + str(hour) 
    else:
      hour = str(hour)

    if minute < 10:
      minute = "0" + str(minute)
    else:
      minute = str(minute)

    ## Bring back the ax features
    # Add coastlines and states
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth = 2)
    ax.add_feature(cfeature.STATES.with_scale('50m'))
    
    # Set limits in lat/lon space
    ax.set_extent([-121, -114, 32, 36]) # SoCal

    ## Add contours, mesh, and text
    # Plot the NCFR "slices"
    contour = ax.contour(x, y, 1 * (labeled_cores > 0), colors = ['k',], linewidths = .5, linestyles = 'solid', \
        zorder = 5) 

    # Add colormesh (radar reflectivity) 
    mesh = ax.pcolormesh(x, y, new_refs, cmap = ref_cmap, norm = ref_norm, zorder = 2) 

    # Add text
    date_string = year + "-" + month + "-" + day + " " + hour + ":" + minute + "Z"
    text = ax.text(0.7, 0.02, date_string, transform = ax.transAxes, fontdict = {'size': 16})

    # Add colormesh and text as artists
    ax.add_artist(mesh)
    ax.add_artist(text)

    return contour


  # Call animate function
  # ani = FuncAnimation(fig, animate_contour, interval = 100, frames = len(ref_refs))
  # ani.save("./plots/20170217_18.gif", writer = PillowWriter(fps = 1))  

if __name__ == '__main__':
  main() 

# Yes!!
# Now find area of intersecting polygons and reorganize/clean code
