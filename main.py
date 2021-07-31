import geopandas
import matplotlib.pyplot as plt
# from mpl_toolkits.basemap import Basemap # Don't use; deprecated!!!
from netCDF4 import Dataset
import numpy as np

## ** Main Function where everything happens **
def main():
  # Load NEXRAD data from netcdf4 file
  ncfile = '/media/jntp/D2BC15A1BC1580E1/NCFRs/20170217_18.nc'
  nexdata = Dataset(ncfile, mode = 'r')

  # Get data from netcdf file
  lons = nexdata['Longitude'][:][:]
  latd = nexdata['Latitude'][:][:]
  print(lons)

  # Get some parameters for a Stereographic Projection
  # lon_0 = lons.mean()
  # lat_0

if __name__ == '__main__':
  main() 

# Check cartopy documentation??? Maybe able to find more about obtaining x, y variables
