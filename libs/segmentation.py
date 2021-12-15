import math
import numpy as np
from geopy.distance import geodesic 
from scipy.ndimage import binary_closing, binary_dilation
from scipy.ndimage.measurements import label
from skimage.morphology import disk, remove_small_objects
from skimage.measure import regionprops

def get_pixel_dimensions(max_lat, max_lon, min_lat, min_lon, ref_rows, ref_cols):
  # Create two "mock places" to find actual length/width
  place1 = (min_lat, min_lon)
  place2 = (max_lat, min_lon)

  ## Find pixel length (kms/pixel): actual_length / ref_rows
  # Find actual geodesic length of the two mock places
  actual_length = geodesic(place1, place2).km

  # Find the pixel length
  pixel_length = actual_length / ref_rows

  # Modify one of the "mock places" so latitudes are matching but not longitudes
  place2 = (min_lat, max_lon)

  ## Find pixel width (kms/pixel): actual_width / ref_cols
  # Find actual geodesic width of the two mock places
  actual_width = geodesic(place1, place2).km

  # Find the pixel width
  pixel_width = actual_width / ref_cols 

  return pixel_length, pixel_width 

def find_convective_cells(refs, min_ref = 45, min_size = 100):
  conv_refs = np.uint8(refs >= 45)  

  # Label features in the intense_refs array where intense pixels touch
  labeled_cells, num_feats_refs = label(conv_refs, np.ones((3, 3), dtype = int)) 

  # Remove small objects  
  labeled_cells = remove_small_objects(labeled_cells, min_size = min_size, connectivity = 2)

  return labeled_cells 

def close_holes(labeled_cells, conv_buffer):
  return binary_closing(labeled_cells > 0, structure = np.ones((3, 3)), iterations = conv_buffer)

def remove_wide_cells(refs, labeled_cells, max_width = 65): # 20 pixels is ~5.5 km
  labeled_image, num_feats = label(1 * (labeled_cells > 0), np.ones((3, 3))) 
  
  # Measure properties of the regions in the labeled image 
  regions = regionprops(labeled_image, intensity_image = refs)

  # Find the maximum width of each region
  for region in regions:
    coords = region.coords 

    # Extract unique y values (rows) of each region
    y_vals = np.unique(coords[:, 0])  

    # Create empty list to store widths of each "line" in a region
    widths = []

    # Find the width for each unique y value
    for y in y_vals: 
      # Extract x values that correspond to the y value
      y_array = np.where(coords == y)  
      y_indices = y_array[0] 
      y_coords = coords[y_indices]
      x_vals = y_coords[:, 1] # corresponds to cols   

      # Subtract different highest and lowest x value to find width
      width = max(x_vals) - min(x_vals)  
      widths.append(width) 
  
    # Check to see if region widths exceeds threshold
    if max(widths) > max_width: 
      # Set the region of the labeled image equal to zero if max width exceeds threshold
      ymin, xmin = np.min(region.coords[:, 0]), np.min(region.coords[:, 1])
      y, x = np.where(region.intensity_image > 0)
      labeled_image[ymin + y, xmin + x] = 0

  return labeled_image

def remove_adjacent_cells(refs, labeled_cells, max_dist = 100, min_slope = 1, y_thresh = 150):
  labeled_image, num_feats = label(1 * (labeled_cells > 0), np.ones((3, 3))) 
  regions = regionprops(labeled_image, refs)

  # Intialize list to store centroids
  centroids = [] 
  y_centroids = []
  x_centroids = []
  check_regions = [] # centroids of which to check "alignment" (along x axis) of other points

  for region in regions:
    centroids.append(region.centroid)
    y_centroids.append(region.centroid[0])
    x_centroids.append(region.centroid[1])
 
    for i, centroid in enumerate(centroids):
      # Don't calculate distance on first iteration, otherwise will prompt error
      if i == 0:
        continue

      # Check the distance of current centroid and all centroids
      if math.dist(region.centroid, centroids[i]) != 0 and math.dist(region.centroid, centroids[i]) <= max_dist:
        # Check if slope magnitude falls within acceptable threshold
        y1, x1 = region.centroid
        y2, x2 = centroids[i]
        slope = (y2 - y1) / (x2 - x1) # calculate slope

        if abs(slope) < min_slope:
          check_regions.append(region)

  # Convert to numpy array; used for array operations 
  y_centroids = np.array(y_centroids)

  # Check regions that meet the above criteria and see if they are aligned with other cells on y-axis
  for region in check_regions:
    lower_bound = region.centroid[1] - 50
    upper_bound = region.centroid[1] + 50
  
    aligned_points = np.where((x_centroids >= lower_bound) & (x_centroids <= upper_bound))[0]
    aligned_pts = np.array(aligned_points) 
    
    y_dist = abs(region.centroid[0] - y_centroids[aligned_pts]) 
    y_dist = y_dist[y_dist != 0] # Delete any zeros in array, which likely references same point 

    if len(aligned_pts > 1) and any(y_dist < y_thresh):
      continue
    else:
      # Set the region of the labeled image equal to zero if max width exceeds threshold
      ymin, xmin = np.min(region.coords[:, 0]), np.min(region.coords[:, 1])
      y, x = np.where(region.intensity_image > 0)
      labeled_image[ymin + y, xmin + x] = 0

  return labeled_image
      
# Merges cells within a specified search radius
def connect_cells(labeled_image, core_buffer):
  return binary_closing(labeled_image > 0, structure = disk(3), iterations = core_buffer)

# Connect cores if wtihin certain distance (gaps), checks to see if the axis falls within NCFR criteria
def check_axis(refs, labeled_cells, min_length = 250): # 80 pixels is ~20 km 
  # Based on find_lines in Haberlie and Ashley (2018)
  labeled_ncfr, num_feats = label(1 * (labeled_cells > 0), np.ones((3, 3)))
  regions = regionprops(labeled_ncfr, intensity_image = refs)

  for region in regions:
    # Check if axis length of each feature is lower than minimum length
    if region.major_axis_length < min_length:
      # Set all pixels within the feature equal to zero
      ymin, xmin = np.min(region.coords[:, 0]), np.min(region.coords[:, 1])
      y, x = np.where(region.intensity_image > 0)
      labeled_ncfr[ymin + y, xmin + x] = 0

  return labeled_ncfr



  

# Next step write extract_cores function
# Clean code and annotate properly (write function for the deletion of regions) 
