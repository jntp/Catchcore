import numpy as np
from skimage.measure import regionprops

def extract_core_boundaries(refs, labeled_cores):
  regions = regionprops(labeled_cores, intensity_image = refs)

  for region in regions:
    print(region.label) 
    exterior_coords = []

    for coord in region.coords:
      i = coord[0]
      j = coord[1]

      # [i-1][j-1] [i-1][j] [i-1][j+1] [i][j-1] [i][j+1] [i+1][j-1] [i+1][j] [i+1][j+1]
      i_0s = []
      j_0s = []

      # Check for "corner cases"
      if i == 0 or j == 0:
        # "Label" the exterior boundary straightaway 
        exterior_coords.append([i, j]) 
        continue
      elif i == 0 and j == labeled_cores.shape[1] - 1:
        exterior_coords.append([i, j]) 
        continue
      elif i == labeled_cores.shape[0] - 1 and j == 0:
        exterior_coords.append([i, j]) 
        continue
      elif i == labeled_cores.shape[0] - 1 and j == labeled_cores.shape[1] - 1:
        exterior_coords.append([i, j]) 
        continue
      else: # Non-corner cases
        # Check for "edge cases"
        if i == 0:
          # [i][j-1] [i][j+1] [i+1][j-1] [i+1][j] [i+1][j+1]
          i_0s = [i, i, i+1, i+1, i+1] 
          j_0s = [j-1, j+1, j-1, j, j+1]
        elif j == 0:
          # [i-1][j] [i-1][j+1] [i][j+1] [i+1][j] [i+1][j+1]
          i_0s = [i-1, i-1, i, i+1, i+1]
          j_0s = [j, j+1, j+1, j, j+1] 
        elif i == labeled_cores.shape[0] - 1:
          # [i-1][j-1] [i-1][j] [i-1][j+1] [i][j-1] [i][j+1] 
          i_0s = [i-1, i-1, i-1, i, i]
          j_0s = [j-1, j, j+1, j-1, j+1] 
        elif j == labeled_cores.shape[1] - 1:
          # [i-1][j-1] [i-1][j] [i][j-1] [i+1][j-1] [i+1][j]
          i_0s = [i-1, i-1, i, i+1, i+1]
          j_0s = [j-1, j, j-1, j-1, j]
        else: # "Normal" cases
          i_0s = [i-1, i-1, i-1, i, i, i+1, i+1, i+1]
          j_0s = [j-1, j, j+1, j-1, j+1, j-1, j, j+1] 
      
      i_0s = np.array(i_0s)
      j_0s = np.array(j_0s)

      # print(i_0s)
      # print(j_0s)

      # print(labeled_cores[i_0s, j_0s])

      # Check if any of the elements in the array equal 0 or the empty spacei
      # Add to external boundary

  # -2147483648 is empty space 
  # Check if each coord stands on its own
  # Construct a new labeled array of the boundary
