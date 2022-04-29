import numpy as np
import geopandas as gpd
from skimage.measure import find_contours
from shapely.geometry import Polygon
from shapely.ops import polygonize_full

## "Auxiliary functions" used to shorten code

def polygonize_watershed_boundary(watershed_boundary, code = 0):
  # Check if default or "work around" method would be used for polygonizing watershed
  if code == 0: # default
    # Polygonize the watershed boundary (currently a linearring or linestring)
    polygons, dangles, cuts, invalids = polygonize_full(watershed_boundary)

    # Return a set of polygons
    return polygons
  elif code == 1: # work around if default method does not work (i.e. SD watershed)
    # Get the coordinates of watershed boundary
    watershed_coords = watershed_boundary.coords

    # Make a polygon out of the boundaries
    polygon = Polygon(watershed_coords)

    # Return a SINGLE polygon
    return polygon

## Intersection Steps

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

    # Check if the polygon is valid
    if shapely_contour.is_valid:
      # Store geometry in array
      shapely_contours.append(shapely_contour)

  return shapely_contours

def find_intersection(polygon_cores, watershed_boundary, code = 0):
  # Create Geoseries from polygon cores; will be needed for intersection 
  cores_gs = gpd.GeoSeries(polygon_cores)

  # Polygonize the watershed boundary (currently a linearring or linestring)
  watershed_polygon = polygonize_watershed_boundary(watershed_boundary, code)

  # Get intersection
  intersections = cores_gs.intersection(watershed_polygon)

  return watershed_polygon, intersections

def check_area_intersections(intersections, watershed_polygon, threshold = 0.05):
  proportions = intersections.area / watershed_polygon.area
  is_sufficient = False

  # Check if sum of proportions meets the threshold
  if sum(proportions) >= threshold:
    is_sufficient = True

  return sum(proportions), is_sufficient

# Later Write an SD specific intersection function that seeks if the bounds of NCFR WILL intersect with
# bounds of watershed. But first write code for Propagation statistics
