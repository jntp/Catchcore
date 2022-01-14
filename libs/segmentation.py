import math
import numpy as np
from geopy.distance import geodesic 
from scipy.ndimage import binary_closing
from scipy.ndimage.measurements import label
from skimage.morphology import disk, remove_small_objects
from skimage.measure import regionprops 

## "Auxilliary functions" used to shorten code

def get_pixel_dimensions(max_lat, max_lon, min_lat, min_lon, ref_rows, ref_cols):
  """
  Finds the conversion between latitude/longitude and pixel width/length. 
  
  Parameters:
  max_lat - the greatest latitude
  max_lon - the greatest longitude
  min_lat - the smallest latitude
  min_lon - the smallest longitude
  ref_rows - total number of rows in 2D array
  ref_cols - total number of columns in 2D array

  Returns:
  pixel_length - length of a pixel (km/pixel)
  pixel_width - width of a pixel (km/pixel)
  """

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

def remove_region(region, labeled_image):
  """
  Removes a region by setting all of its x and y coordinates equal to zero in a labeled image
  
  Parameters:
  region - a single region or element in regions array; must be initialized using skimage.measure.regionprops()
  labeled_image - a labeled image, initialized using scipy.ndimage.measurements.label()

  Returns:
  labeled_image - the updated labeled image with the region removed
  """

  ymin, xmin = np.min(region.coords[:, 0]), np.min(region.coords[:, 1])
  y, x = np.where(region.intensity_image > 0)
  labeled_image[ymin + y, xmin + x] = 0

  return labeled_image 

## Segmentation steps

def find_convective_cells(refs, min_ref = 45, min_size = 100):
  """
  Extracts convective cells in a radar image based on a specified threshold 

  Parameters:
  refs - radar image (reflectivity)
  min_ref - minimum convective reflectivity threshold (default: 45 dbZ)
  min_size - minimum size of a convective cell (default: 100 pixels)

  Returns:
  labeled_cells - labeled image of cells that meet minimum reflectivity and size requirements
  """

  # Extract reflectivity based off threshold
  conv_refs = np.uint8(refs >= min_ref)  

  # Label features in the intense_refs array where intense pixels touch
  labeled_cells, num_feats_refs = label(conv_refs, np.ones((3, 3), dtype = int)) 

  # Remove small objects  
  labeled_cells = remove_small_objects(labeled_cells, min_size = min_size, connectivity = 2)

  return labeled_cells 

def close_holes(labeled_cells, conv_buffer):
  """
  Performs a binary closing (erosion then dilation) of an labeled image. Removes "small holes" on each label based
  on specified size.

  Parameters:
  labeled_cells - labeled image of convective cells
  conv_buffer - size of holes to close

  Returns:
  closing of the labeled image
  """

  return binary_closing(labeled_cells > 0, structure = np.ones((3, 3)), iterations = conv_buffer)

def remove_wide_cells(refs, labeled_cells, mode = 0, max_width = 80): 
  """
  Removes convective cells wider than a specified width.

  Parameters:
  refs - radar image (reflectivity)
  labeled_cells - labeled image of convective cells
  mode - 0 looks at greatest width of cell; 1 looks at width of centroid (default: 0)
  max_width - maximum width allowed for a convective cell (default: 65 pixels) 

  Returns:
  labeled_image - updated labeled image of convective cells with wide cells removed
  """

  # Create labeled image given labeled cells
  labeled_image, num_feats = label(1 * (labeled_cells > 0), np.ones((3, 3))) 
  
  # Measure properties of the regions in the labeled image 
  regions = regionprops(labeled_image, intensity_image = refs)

  # Find the maximum width of each region
  for region in regions:
    coords = region.coords 
    y_centroid = int(region.centroid[0]) # new code 

    # Extract unique y values (rows) of each region
    y_vals = np.unique(coords[:, 0])  

    if mode == 0:
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
        labeled_image = remove_region(region, labeled_image)
    elif mode == 1:
      # New Code!!!
      test_indices = np.where(coords == y_centroid)[0]
      test_coords = coords[test_indices]
      test_vals = test_coords[:, 1] 
    
      test_width = max(test_vals) - min(test_vals)
    
      if test_width > max_width:
        labeled_image = remove_region(region, labeled_image)  
     
  return labeled_image

def remove_adjacent_cells(refs, labeled_cells, min_size = 200, max_dist = 200, min_slope = 1, y_thresh = 25):
  """
  Removes cells that are adjacent (on x axis) to the "suspected" NCFR core. First checks if the centroids of 
  the cells fall within a certain distance. Then checks the slope to determine if the centroids of the cells 
  lie horizontally from each other. The last step involves checking whether each cells is "aligned" with 
  other cells along the vertical (y) axis. Removes the cell that is not aligned with other cells on y-axis.

  Parameters:
  refs - radar image (reflectivity)
  labeled_cells - labeled image of convective cells
  max_dist - maximum distance between centroids for cells to be considered "adjacent" to each other (default: 100px)
  min_slope - minimum slope for cells to be considered NOT "adjacent" to each other along the x-axis (default: 1)
  y_thresh - maximum distance along y-axis for cells to be considered associated with each other (default: 150px)

  Returns:
  labeled_image - updated labeled image with adjacent cell removed
  """

  # Create labeled image and regions
  labeled_image, num_feats = label(1 * (labeled_cells > 0), np.ones((3, 3))) 
  regions = regionprops(labeled_image, refs)

  # Remove small objects
  labeled_image = remove_small_objects(labeled_image, min_size = min_size, connectivity = 2)

  # Intialize list to store centroids and centroid information
  centroids = [] 
  y_centroids = []
  x_centroids = []
  check_regions = [] # centroids of which to check "alignment" (along x axis) of other points

  # Loop through regions and compare distance and slope of different regions/cells
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

        # If the slope magnitude lies below minimum value (suggesting adjacency), flag region for further investigation
        print(slope)       
        if abs(slope) < min_slope:
          check_regions.append(region)

  # Convert to numpy array; used for array operations 
  y_centroids = np.array(y_centroids)

  # Check regions that meet the above criteria and see if they are aligned with other cells on y-axis
  for region in check_regions:
    # Create lower and upper bound x values, will be used to determine if cells like on vertical (y) axis
    lower_bound = region.centroid[1] - 50
    upper_bound = region.centroid[1] + 50
  
    # Look for points that lie within the upper and lower x boundaries
    aligned_points = np.where((x_centroids >= lower_bound) & (x_centroids <= upper_bound))[0]
    aligned_pts = np.array(aligned_points) 
    
    # Obtains the vertical distance of points "aligned" on y-axis
    y_dist = abs(region.centroid[0] - y_centroids[aligned_pts]) 
    y_dist = y_dist[y_dist != 0] # Delete any zeros in array, which likely references same point 

    # Check if cell is aligned with at least once cell and the aligned cells are within threshold distance
    if len(aligned_pts > 1) and any(y_dist < y_thresh):
      # Skip to next iteration; region/cell passes test
      continue
    else:
      # Set the region of the labeled image equal to zero if max width exceeds threshold
      labeled_image = remove_region(region, labeled_image)

  # Remove small objects  
  labeled_image = remove_small_objects(labeled_image, min_size = min_size * 4, connectivity = 2)

  # Test
  labeled_image = check_axis(refs, labeled_image)

  return labeled_image
      
def connect_cells(labeled_image, core_buffer):
  """
  Merges cells within a specified search radius.

  Parameters:
  labeled_image - labeled image of convective cells
  core_buffer - maximum size of distance between cores that could qualify as NCFR (also known as "gaps")

  Returns:
  closing of labeled image that would result in the merging of cells meeting the criteria
  """

  return binary_closing(labeled_image > 0, structure = disk(3), iterations = core_buffer)

# Connect cores if wtihin certain distance (gaps), checks to see if the axis falls within NCFR criteria
def check_axis(refs, labeled_cells, min_length = 250): 
  """
  Checks if the major axis (or vertical axis of NCFR) of each convective region meets minimum length to be
  considered an NCFR. Eliminates all regions that don't meet the criteria.

  Parameters:
  refs - radar image (reflectivity)
  labeled_cells - labeled image of convective cells
  min_length - minimum length to be considered an NCFR

  Returns:
  labeled_ncfr - updated labeled image that shows the NCFR
  """

  # Create labeled image and regions
  labeled_ncfr, num_feats = label(1 * (labeled_cells > 0), np.ones((3, 3)))
  regions = regionprops(labeled_ncfr, intensity_image = refs)

  # Check length of each region
  for region in regions:
    # Check if axis length of each feature is lower than minimum length
    if region.major_axis_length < min_length:
      # Set all pixels within the feature equal to zero
      labeled_ncfr = remove_region(region, labeled_ncfr)

  return labeled_ncfr

def extract_cores(refs, labeled_ncfr, conv_buffer, min_ref = 45, min_size = 800):
  """
  Extracts the cores from a labeled NCFR image.

  Parameters:
  refs - radar image (reflectivity)
  labeled_ncfr - labeled image of an NCFR
  conv_buffer - size of holes to close 

  Returns:
  labeled_cores - updated labeled image of the NCFR cores
  """

  # Look for pixels that fail to meet convective criteria; remove them (or set equal to 0)
  y, x = np.where(refs < min_ref)
  labeled_ncfr[y, x] = 0
  labeled_cores = close_holes(labeled_ncfr, conv_buffer) # remove small holes from labeled regions

  # Remove small objects  
  # labeled_cores = remove_small_objects(labeled_cores, min_size = min_size, connectivity = 2)
  # PERHAPS MOVE THIS TO REMOVE_ADJACENT_CELLS!!!

  # Remove small holes
  labeled_cores = close_holes(labeled_cores, 3) 

  # Remove large objects
  # The thinking is remove large objects... then create a Step 8 which produces labeled NCFR
  # labeled_cores = remove_wide_cells(refs, labeled_cores, 1)

  return labeled_cores 
