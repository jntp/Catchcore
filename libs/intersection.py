import numpy as np
import geopandas as gpd
from skimage.measure import find_contours
from shapely.geometry import Polygon 

def get_core_contours(labeled_cores, lons, lats):
  # Find contours from labeled_cores
  contours = find_contours(labeled_cores, 0)

  # Create array to store contours of cores
  shapely_contours = []

  for contour in contours:
    # Extract indices of each contour; make accessible by converting to numpy array
    i = np.array(contour[:, 0], dtype = int)
    j = np.array(contour[:, 1], dtype = int)

    # Extract latitude and longitude coordinates corresponding NEXRAD crs
    lons_contour = lons[i, j]
    lats_contour = lats[i, j]

    # Arrange coordinates as tuple and store in array
    contour_points = []
    
    for k, lon in enumerate(lons_contour):
      contour_points.append((lons_contour[k], lats_contour[k]))

    # Make the shapely polygon from contour points
    shapely_contour = Polygon(contour_points) 

    # Store geometry in array
    shapely_contours.append(shapely_contour)

  return shapely_contours

def find_intersection:
  # Create Geoseries from polygon cores; will be needed for intersection 
  print("Test") 

