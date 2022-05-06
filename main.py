import os
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
import pandas as pd 
import metpy.plots as mpplots 
from libs.segmentation import *
from libs.intersection import *
from libs.propagation_stats import *

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
def segmentation(refs, i = 0, core_buffer = 30, conv_buffer = 3):
  """
    i - timestep of the radar loop, only needed for the animation (default: 0)
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

  # Check if the algorithm returns nothing and that the timestep is later in the animation
  if not any(map(any, labeled_cores)) and i >= 40:
    # Run segmentation procedure again but stop at Step 4 and remove small cells
    conv_cells = find_convective_cells(refs)
    closed_cells = close_holes(conv_cells, conv_buffer)
    narrow_conv_cells = remove_wide_cells(refs, closed_cells, 120) # increase width of cells
    narrow_conv_cells2 = remove_adjacent_cells(refs, narrow_conv_cells) 
    labeled_cores = remove_small_cells(refs, narrow_conv_cells2)

  return labeled_ncfr, labeled_cores

def intersection(shapely_cores, watershed, threshold = 0.03, code = 0):
  if code == 0:
    watershed_boundary = watershed.geometry[1].boundary 
    watershed_polygon, watershed_intersections = find_intersection(shapely_cores, watershed_boundary)
    watershed_proportion, watershed_cross = check_area_intersections(watershed_intersections, watershed_polygon, threshold)

    return watershed_polygon, watershed_intersections, watershed_proportion, watershed_cross
  elif code == 1:
    watershed_boundary = watershed.geometry[1].boundary 
    water_polygon, watershed_intersections = find_intersection(shapely_cores, watershed_boundary, 1)
    watershed_polygon = gpd.GeoSeries(water_polygon) # turn water_polygon into geoseries since it is currently 1 polygon and not iterable
    watershed_proportion, watershed_cross = check_area_intersections(watershed_intersections, water_polygon, threshold)

    return watershed_polygon, watershed_intersections, watershed_proportion, watershed_cross

def plot_watershed_polygons(ax, *polygons):
  for polygon in polygons:
    ax.add_geometries(polygon, crs = ccrs.PlateCarree(), zorder = 1, facecolor = 'red', edgecolor = 'red') 

def plot_intersections(ax, *intersecting_polygons):
  for intersecting_polygon in intersecting_polygons:
     ax.add_geometries(intersecting_polygon, crs = ccrs.PlateCarree(), zorder = 3, facecolor = 'yellow', edgecolor = 'yellow') 

# Plots a single image
def plot_single(ax, x, y, ref_cmap, ref_norm, new_refs, labeled_image):
  # Plot the NCFR "slices"
  ax.contour(x, y, 1 * (labeled_image > 0), colors = ['k',], linewidths = .5, linestyles = 'solid', \
      zorder = 5) 
 
  # Add colormesh (radar reflectivity) 
  ax.pcolormesh(x, y, new_refs, cmap = ref_cmap, norm = ref_norm, zorder = 2)

  plt.show() 

def format_time(years_i, months_i, days_i, hours_i, minutes_i):
  # Extract "singular" time from time variables
  year = str(int(years_i))
  month = str(int(months_i))
  day = str(int(days_i))
  hour = int(hours_i)
  minute = int(minutes_i)

  # Format the hour and minute so it will display properly
  if hour < 10:
    hour = "0" + str(hour) 
  else:
    hour = str(hour)

  if minute < 10:
    minute = "0" + str(minute)
  else:
    minute = str(minute)

  return year, month, day, hour, minute

# Initialize the intersection and/or propagation statistics process
def initiation(refs, ts, lons, lats):
  ref_ref = refs[ts]
  labeled_ncfr, labeled_cores = segmentation(ref_ref, ts)
  shapely_contours = get_core_contours(labeled_cores, lons, lats)

  return ref_ref, labeled_ncfr, labeled_cores, shapely_contours

# Code 0 for default output, code 1 for NCFR propagation statistics
def out_to_csv(df, ref_year, ref_month, ref_day, timestep, code = 0):
  # Create a reference date for the file name 
  ref_date = str(int(ref_year)) + str(int(ref_month)) + str(int(ref_day))

  # Directory path
  out_dir = "./outputs/" + ref_date

  # Check if the directory exists; create new one if not
  isdir = os.path.isdir(out_dir)
  if not isdir:
    os.mkdir(out_dir)

  # Create output file path directory
  if code == 0:
    out_file = out_dir + "/" + ref_date + "_" + str(timestep) + ".csv"
  elif code == 1:
    out_file = out_dir + "/" + ref_date + "_" + "Propagation_Stats.csv"

  # Output dataframe to csv file 
  df.to_csv(out_file)


## ** Main Function where everything happens **
def main():
  # Load watersheds
  sepulveda = gpd.read_file("LA_R_A_Sepulveda_Dam.geojson")
  whittier = gpd.read_file("Rio_Hondo_AB_Whittier_Narrows_Dam.geojson") 
  santa_ana = gpd.read_file("santa_ana_r_a.geojson")
  san_diego = gpd.read_file("SD_R_A_Fashion_Valley.geojson") 

  # Load NEXRAD data from netcdf4 file
  ncfile = '/media/jntp/D2BC15A1BC1580E1/NCFRs/20170217_18.nc'
  nexdata = Dataset(ncfile, mode = 'r')
  print(nexdata)

  # Get data from netcdf file
  lons = nexdata['Longitude'][:][:]
  lats = nexdata['Latitude'][:][:] 
 
  ref_refs = nexdata['Reflectivity'][:] # reflectivities, used for animation
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

  "---------------Below is the code for running steps separately---------------"

  ## Segmentation - Delineate the NCFR cores; comment out if animating
  # Initiation
  # ts = 42
  # ref_ref, labeled_ncfr, labeled_cores, shapely_contours = initiation(ref_refs, ts, lons, lats) 

  ## Intersection - Find Intersections between cores and watershed; comment out if animating
  # Get polygon, intersections, proportion of intersection, and "cross" variable that indicates "true" intersection
  # sepulveda_polygon, sepulveda_intersections, sepulveda_proportion, sepulveda_cross = intersection(shapely_contours, sepulveda, 0.08)
  # whittier_polygon, whittier_intersections, whittier_proportion, whittier_cross = intersection(shapely_contours, whittier, 0.1)
  # santa_ana_polygon, santa_ana_intersections, santa_ana_proportion, santa_ana_cross = intersection(shapely_contours, santa_ana)
  # san_diego_polygon, san_diego_intersections, san_diego_proportion, san_diego_cross = intersection(shapely_contours, san_diego, 0.02, 1)

  # Plot the shapely geometries and intersections
  # ax.add_geometries(shapely_contours, crs = ccrs.PlateCarree(), zorder = 2, facecolor = 'green', edgecolor = 'green')
  # plot_watershed_polygons(ax, sepulveda_polygon, whittier_polygon, santa_ana_polygon, san_diego_polygon)
  # plot_intersections(ax, sepulveda_intersections, whittier_intersections, santa_ana_intersections, san_diego_intersections)

  ## Propagation Statistics
  # Initiation of another timeframe for comparison
  # ref_ref0, labeled_ncfr0, labeled_cores0, shapely_contours0 = initiation(ref_refs, 30, lons, lats) 

  # Manually track a core and obtain two points in different timeframes (spaced 12 timesteps or 1 hr apart)
  # Points should be an estimation of the tracked core's centroid
  # Then get the exact coordinates of the two centroids
  # centroid1 = get_closest_centroid(ref_ref0, labeled_cores0, (-117.698, 33.071), lats, lons)
  # centroid2 = get_closest_centroid(ref_ref, labeled_cores, (-117.233, 32.940), lats, lons)
  
  # Calculate the distance traveled, forward azimuth, and speed from the two coordinates
  # distance_km, fwd_azimuth, speed_m_s = calculate_stats(centroid1, centroid2)

  ## Maximum reflectivity
  # Get the maximum reflectivity
  # max_ref = get_max_ref(ref_ref, labeled_cores)

  "---------------Below is running everything altogether---------------"

  ## Save the Statistics
  # Create empty lists to store stat variables; don't comment this one out!
  date = []
  time = []
  sepulveda_prop = []
  sepulveda_crossing = []
  whittier_prop = []
  whittier_crossing = []
  santa_ana_prop = []
  santa_ana_crossing = []
  san_diego_prop = []
  san_diego_crossing = []
  max_reflectivity = []
  timestep = []
 
  ## Output the statistics to csv file
  def run_output_stats(i):
    ## Initiate and run segmentation algorithm for specific time frame 
    ref_ref, labeled_ncfr, labeled_cores, shapely_contours = initiation(ref_refs, i, lons, lats) 

    ## Obtain polygons and intersections
    # Get polygon, intersections, proportion of intersection, and "cross" variable that indicates "true" intersection
    sepulveda_polygon, sepulveda_intersections, sepulveda_proportion, sepulveda_cross = intersection(shapely_contours, sepulveda, 0.08)
    whittier_polygon, whittier_intersections, whittier_proportion, whittier_cross = intersection(shapely_contours, whittier, 0.1)
    santa_ana_polygon, santa_ana_intersections, santa_ana_proportion, santa_ana_cross = intersection(shapely_contours, santa_ana)
    san_diego_polygon, san_diego_intersections, san_diego_proportion, san_diego_cross = intersection(shapely_contours, san_diego, 0.02, 1)

    ## Maximum reflectivity
    # Get the maximum reflectivity
    max_ref = get_max_ref(ref_ref, labeled_cores)
 
    # Format time variables and convert to string
    year, month, day, hour, minute = format_time(years[i], months[i], days[i], hours[i], minutes[i]) 

    # Create/concatenate a date and time variable
    date_str = year + "-" + month + "-" + day
    time_str = hour + ":" + minute

    # Append stat variables to lists
    date.append(date_str)
    time.append(time_str) 
    sepulveda_prop.append(sepulveda_proportion)
    sepulveda_crossing.append(sepulveda_cross)
    whittier_prop.append(whittier_proportion) 
    whittier_crossing.append(whittier_cross)
    santa_ana_prop.append(santa_ana_proportion)
    santa_ana_crossing.append(santa_ana_cross)
    san_diego_prop.append(san_diego_proportion)
    san_diego_crossing.append(san_diego_cross)
    max_reflectivity.append(max_ref)
    timestep.append(i)
 
    # Create dataframe from list and then convert it to csv file
    df = pd.DataFrame({'Date': date, \
      'Time (Z)': time, \
      'Sepulveda Proportion': sepulveda_prop, \
      'Sepulveda Crossing': sepulveda_crossing, \
      'Whittier Proportion': whittier_prop, \
      'Whittier Crossing': whittier_crossing, \
      'Santa Ana Proportion': santa_ana_prop, \
      'Santa Ana Crossing': santa_ana_crossing, \
      'San Diego Proportion': san_diego_prop, \
      'San Diego Crossing': san_diego_crossing, \
      'Maximum Reflectivity (dbZ)': max_reflectivity, \
      'Time Step Index': timestep})

    # Output the statistics to csv file
    out_to_csv(df, years[0], months[0], days[0], i)

  ## Run the run_output_stats function
  # Define and starting and ending timestep
  start_ts = 0
  end_ts = 60

  # Check whether to run the function once or multiple times
  if start_ts == end_ts:
    run_output_stats(start_ts) # run the function once
  else:
    for i in range(start_ts, end_ts):
      run_output_stats(i)

  "--------------Below is the code for outputting NCFR Propagation Statistics---------------"

  ## Output NCFR Propagation Stats to csv file
  # Create another group of empty lists to store the NCFR Propagation Stats (storm speed, azimuth, direction)
  start_timestep = []
  start_centroid = []
  end_timestep = []
  end_centroid = []
  distances_km = []
  fwd_azimuths = []
  speeds_m_s = []

  # Specify a starting and ending timestep (normally should be spaced 12 timesteps apart or a 1 hour period)
  ts_beg = 30
  ts_end = 42 # this is the last timestep; no need to "add 1"

  # NCFR Propagation Statistics - run the steps for calculating the statistics
  def run_stats(coord0, coord1):
    """
      Parameters:
      coord0, coord1 - tuple in the form of (lon, lat), user estimated coordinates of a tracked core's centroid
    """
    # Initiate the segmetation process
    ref_ref0, labeled_ncfr0, labeled_cores0, shapely_contours0 = initiation(ref_refs, ts_beg, lons, lats)
    ref_ref1, labeled_ncfr1, labeled_cores1, shapely_contours1 = initiation(ref_refs, ts_end, lons, lats)

    # Manually track a core and obtain two points in different timeframes (spaced 12 timesteps or 1 hr apart)
    # Points should be an estimation of the tracked core's centroid
    # Then get the exact coordinates of the two centroids
    centroid0 = get_closest_centroid(ref_ref0, labeled_cores0, coord0, lats, lons)
    centroid1 = get_closest_centroid(ref_ref1, labeled_cores1, coord1, lats, lons)
  
    # Calculate the distance traveled, forward azimuth, and speed from the two coordinates
    distance_km, fwd_azimuth, speed_m_s = calculate_stats(centroid0, centroid1)

    # Append to lists
    start_timestep.append(ts_beg)
    start_centroid.append(centroid0)
    end_timestep.append(ts_end)
    end_centroid.append(centroid1) 
    distances_km.append(distance_km)
    fwd_azimuths.append(fwd_azimuth)
    speeds_m_s.append(speed_m_s)

    # Create dataframe from list and then convert it to csv file
    df = pd.DataFrame({'Starting Timestep': start_timestep, \
        'Centroid Start Position': start_centroid, \
        'Ending Timestep': end_timestep, \
        'Centroid End Position': end_centroid, \
        'Distance (km)': distances_km, \
        'Azimuth (deg)': fwd_azimuth, \
        'Speed (m/s)': speed_m_s})

    # Output Propagation Statistics to csv file
    out_to_csv(df, years[0], months[0], days[0], ts_beg, 1) # code = 1 for NCFR Propagation Statistics

  # Specify two coordinates (tuples) that our estimation of the tracked core's centroid over time
  coord0 = (-117.698, 33.071)
  coord1 = (-117.233, 32.940) 

  # Call the run_stats function using coord0 and coord1
  run_stats(coord0, coord1)

  "---------------Below is the Code for Plotting---------------"

  ## Plot a single image (comment out if animating)
  new_refs = new_reflectivity(ref_ref)
  plot_single(ax, x, y, ref_cmap, ref_norm, new_refs, labeled_cores)

  # plt.show()

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
    labeled_ncfr, labeled_cores = segmentation(ref_refs[i], i)

    # Format time variables and convert to string
    year, month, day, hour, minute = format_time(years[i], months[i], days[i], hours[i], minutes[i]) 

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
    date_string = year + "-" + month + "-" + day + " " + hour + ":" + minute + "Z" + " " + "ts=" + i
    text = ax.text(0.7, 0.02, date_string, transform = ax.transAxes, fontdict = {'size': 16})

    # Add colormesh and text as artists
    ax.add_artist(mesh)
    ax.add_artist(text)

    return contour

  ## Animate polygons of watersheds, intersections, and cores
  def animate_geometries(i):
    # Clear axes to prevent contours from clobbering up
    ax.clear()  

    ## Initiation
    # Create transparent background for reflectivity >= 5 dbZ
    new_refs = new_reflectivity(ref_refs[i])

    # Run segmentation algorithm for specific time frame 
    labeled_ncfr, labeled_cores = segmentation(ref_refs[i], i)

    # Obtain shapely geometric contours of the labeled cores
    shapely_contours = get_core_contours(labeled_cores, lons, lats)

    # Format time variables and convert to string
    year, month, day, hour, minute = format_time(years[i], months[i], days[i], hours[i], minutes[i]) 

    ## Bring back the ax features
    # Add coastlines and states
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth = 2)
    ax.add_feature(cfeature.STATES.with_scale('50m'))
    
    # Set limits in lat/lon space
    ax.set_extent([-121, -114, 32, 36]) # SoCal

    ## Obtain polygons and intersections
    # Get polygon, intersections, proportion of intersection, and "cross" variable that indicates "true" intersection
    sepulveda_polygon, sepulveda_intersections, sepulveda_proportion, sepulveda_cross = intersection(shapely_contours, sepulveda, 0.08)
    whittier_polygon, whittier_intersections, whittier_proportion, whittier_cross = intersection(shapely_contours, whittier, 0.1)
    santa_ana_polygon, santa_ana_intersections, santa_ana_proportion, santa_ana_cross = intersection(shapely_contours, santa_ana)
    san_diego_polygon, san_diego_intersections, san_diego_proportion, san_diego_cross = intersection(shapely_contours, san_diego, 0.02, 1)

    ## Plot cores, polygons, intersections, and text
    # Plot the shapely geometries and intersections
    ax.add_geometries(shapely_contours, crs = ccrs.PlateCarree(), zorder = 2, facecolor = 'green', edgecolor = 'green')
    plot_watershed_polygons(ax, sepulveda_polygon, whittier_polygon, santa_ana_polygon, san_diego_polygon)
    plot_intersections(ax, sepulveda_intersections, whittier_intersections, santa_ana_intersections, san_diego_intersections)

    # Add text
    date_string = year + "-" + month + "-" + day + " " + hour + ":" + minute + "Z" + " " + "ts=" + i
    text = ax.text(0.7, 0.02, date_string, transform = ax.transAxes, fontdict = {'size': 16})

    return text
    
  # Call animate function
  ani = FuncAnimation(fig, animate_geometries, interval = 100, frames = len(ref_refs))
  ani.save("./plots/20170217_18_polygon.gif", writer = PillowWriter(fps = 1))  

if __name__ == '__main__':
  main() 

