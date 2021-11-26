import numpy as np
from scipy.ndimage.measurements import label
from skimage.morphology import remove_small_objects 

def find_prelim_cores(refs, min_ref = 45):
  intense_refs = refs[refs >= 45]

  # Label features in the intense_refs array where intense pixels touch
  labeled_refs, num_feats_refs = label(refs, np.ones((3, 3), dtype = int))

  # Remove small objects  
