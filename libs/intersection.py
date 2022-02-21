import numpy as np
from skimage.measure import regionprops

def extract_core_centroids(refs, labeled_cores):
  regions = regionprops(labeled_cores, intensity_image = refs)
  centroids = []

  for region in regions:
    centroids.append(region.centroid) 

  return centroids

# See "Approximate and subdivide polygons" and find_contours in sci-kit-image manual
