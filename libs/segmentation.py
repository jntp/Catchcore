import numpy as np
from scipy.ndimage.measurements import label
from skimage.morphology import remove_small_objects 
from skimage.measure import regionprops

def find_convective_cells(refs, min_ref = 45, min_size = 10):
  intense_refs = refs[refs >= min_ref]

  # Label features in the intense_refs array where intense pixels touch
  labeled_refs, num_feats_refs = label(refs, np.ones((3, 3), dtype = int))

  # Remove small objects  
  remove_small_objects(labeled_refs, min_size = min_size, connectivity = 2, out = labeled_refs)

  return labeled_refs 
  # Measure properties of the regions in the labeled image
  # regions = regionprops(labeled_refs, intensity_image = refs)
  
  # Initialize another labeled image array
  # labeled_image = np.zeros(labeled_refs.shape, dtype = int)



