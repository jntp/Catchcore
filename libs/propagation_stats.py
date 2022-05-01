import math
import numpy as np
from scipy.ndimage.measurements import label
from skimage.measure import regionprops

def get_current_stats(refs, labeled_cores):
  # Get stats of each region
  # center_pt_lon
  # center_pt_lat 
  # max_ref
  # time_index
  labeled_image, num_features = label(1 * (labeled_cores > 0), np.ones((3, 3)))
  regions = regionprops(labeled_image, intensity_image = refs)
  for region in regions:
    print(region.label)
    print(region.centroid)

  print(refs.shape)

# Next steps... use math.floor on row, col of each centroid
# Use np.where to obtain lon, lat values... then store in list
# Get max reflectivity of each region

# Make another function that calculates the average speed, max reflectivity, azimuth, etc.

