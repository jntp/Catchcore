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
  # Create lists to find minimum value later
  centroids = []
  distances = []

  # Get the labeled image and regions
  labeled_image, num_features = label(1 * (labeled_cores > 0), np.ones((3, 3)))
  regions = regionprops(labeled_image, intensity_image = refs)

  # Obtain the distances between the point and each of the core's centroid
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

## Propagation Statistics Steps

# time_minutes is how much time has elapsed between beginning and end of tracking period
def calculate_stats(centroid1, centroid2, time_minutes = 60):
  # Obtain individual latitude and longitude values from the centroids
  lon1, lat1 = centroid1
  lon2, lat2 = centroid2

  # Define default ellipsoid calculations
  geo = pyproj.Geod(ellps = 'WGS84')

  # Get the forward and backward azimuths plus distance (in m) via an inverse transformation
  fwd_azimuth, back_azimuth, distance_m = geo.inv(lon1, lat1, lon2, lat2)  

  # Convert distance to km
  distance_km = distance_m / 1000

  # Convert input minutes to seconds for calculating speed in m/s 
  time_seconds = time_minutes * 60

  # Calculate the speed in m/s
  speed_m_s = distance_m / time_seconds

  return distance_km, fwd_azimuth, speed_m_s

def get_max_ref(refs, labeled_cores):
  # Get the indices where the cores exist
  cores_i, cores_j = np.where(labeled_cores > 0)

  # Extract the reflectivity values
  ref_cores = refs[cores_i, cores_j]

  # Obtain the maximum value
  max_ref = max(ref_cores)

  return max_ref


# Make another function that calculates the average speed, max reflectivity, azimuth, etc.

