import math
import numpy as np
from scipy.ndimage.measurements import label
from skimage.morphology import remove_small_objects 
from skimage.measure import regionprops

def find_convective_cells(refs, min_ref = 45, min_size = 10):
  intense_refs = refs[refs >= min_ref]

  # Label features in the intense_refs array where intense pixels touch
  labeled_cells, num_feats_refs = label(refs, np.ones((3, 3), dtype = int))

  # Remove small objects  
  labeled_cells = remove_small_objects(labeled_refs, min_size = min_size, connectivity = 2)

  return labeled_cells 
  # Measure properties of the regions in the labeled image
  # regions = regionprops(labeled_refs, intensity_image = refs)
  
  # Initialize another labeled image array
  # labeled_image = np.zeros(labeled_refs.shape, dtype = int)

def remove_wide_cells(refs, labeled_cells, max_width = 30): # still need to find conversion
  # General thinking: subtract longest x from centroid x to get width
  # Check if width > 5 km

  # Measure properties of the regions in the labeled image 
  regions = regionprops(labeled_cells, intensity_image = refs)
 
  # Test code of a particular region... delete later
  test_coord = regions[1].coords
  y_vals = np.unique(test_coord[:, 0]) 
  widths = []
  labeled_cells2 = np.zeros(labeled_cells.shape, dtype = int)

  for y in y_vals:
    test_array = np.where(test_coord == y) 
    test_indices = test_array[0]
    new_coord = test_coord[test_indices]
    x_vals = new_coord[:, 1]
    width = max(x_vals) - min(x_vals)
    widths.append(width)
 
  # if max(widths) < max_width:
    # Basically update the actual refs where the labels meet
    # labeled_cells2 += (labeled_cells == region.label) * refs
  
  for region in regions:
    # Rewrite test code here 
