import math
import numpy as np
from geopy.distance import geodesic 
from scipy.ndimage.measurements import label
from skimage.morphology import remove_small_objects 
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

def find_convective_cells(refs, min_ref = 45, min_size = 10):
  intense_refs = refs[refs >= min_ref]

  # Label features in the intense_refs array where intense pixels touch
  labeled_cells, num_feats_refs = label(refs, np.ones((3, 3), dtype = int))

  # Remove small objects  
  labeled_cells = remove_small_objects(labeled_cells, min_size = min_size, connectivity = 2)

  return labeled_cells 

def remove_wide_cells(refs, labeled_cells, max_width = 20): # 20 pixels is ~5.25 km
  # Measure properties of the regions in the labeled image 
  regions = regionprops(labeled_cells, intensity_image = refs)

  # Create new matrix of zeros, which will be a labeled image replacing labeled_cells
  labeled_cells2 = np.zeros(labeled_cells.shape, dtype = int)

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
    if max(widths) < max_width: 
      # Update the labeled cells if width is less than the threshold
      labeled_cells2 += (labeled_cells == region.label) * refs

  return labeled_cells2 

# Next step Connect Gaps with Cores
# Will need to read literature on the general length of gaps
