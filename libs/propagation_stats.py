import math
import pyproj
import numpy as np
from scipy.ndimage.measurements import label
from skimage.measure import regionprops

## "Auxiliary functions" used to shorten code
# Converts (row, col) point to (lon, lat) point
def convert_latlon(point, lats, lons):
  row, col = point
  i = math.floor(row)
  j = math.floor(col) 

  lon = lons[i, j]
  lat = lats[i, j]

  new_point = (lon, lat)

  return new_point

# Get the closest centroid based on a given point
def get_closest_centroid(refs, labeled_cores, point, lats, lons):
  centroids = []
  distances = []

  labeled_image, num_features = label(1 * (labeled_cores > 0), np.ones((3, 3)))
  regions = regionprops(labeled_image, intensity_image = refs)

  for region in regions:
    # Convert to lat lon coordinates
    centroid = convert_latlon(region.centroid, lats, lons)
    centroids.append(centroid)

    distance = math.dist(point, centroid)
    distances.append(distance)

  # Find the index of the minimum distance 
  k = np.argmin(distances) # index of the minimum distance
  
  # Find the centroid with the minimum distance
  closest_centroid = centroids[k]

  return closest_centroid

def get_current_stats(refs, labeled_cores, lats, lons):
  # Get stats of each region
  # center_pt_lon
  # center_pt_lat 
  # time_index
  labeled_image, num_features = label(1 * (labeled_cores > 0), np.ones((3, 3)))
  regions = regionprops(labeled_image, intensity_image = refs)
  for region in regions:
    print(region.label)
    print(region.centroid)
    centroid_lonlat = convert_latlon(region.centroid, lats, lons)
    print(centroid_lonlat)

  print(refs.shape)

def calculate_stats(centroid1, centroid2):
  # Obtain individual latitude and longitude values from the centroids
  lon1, lat1 = centroid1
  lon2, lat2 = centroid2

  # Define default ellipsoid calculations
  geo = pyproj.Geod(ellps = 'WGS84')

  # Get the forward and backward azimuths plus distance (in m) via an inverse transformation
  



# Find centroid with closest latitude and longitude coordinate 

# Next steps... use math.floor on row, col of each centroid
# Use np.where to obtain lon, lat values... then store in list
# Get max reflectivity of each region

# Make another function that calculates the average speed, max reflectivity, azimuth, etc.

