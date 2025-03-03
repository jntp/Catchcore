# Catchcore

## Introduction

Catchcore is a portmanteau between a hydrologic CATCHment and a CORE of a narrow cold frontal rainband (NCFR). An NCFR contains gaps and cores where cores are areas of very intense precipitation and gaps are areas of relatively light precipitation. When an NCFR core intersects with a catchment (or watershed) of interest, they become one--hence the name "Catchcore." This repository contains a Python segmentation algorithm meant to delineate the boundaries of (or "catch") the NCFR "cores" intersecting with CATCHments of interest.

![NCFR_Algorithm_Flow_3_whitebg](https://github.com/user-attachments/assets/cabea51d-a825-49e3-8dad-980a2319d5cb)

## Results

This algorithm follows a seven-step procedure to segment or delineate NCFR cores on Next Generation Weather Radar (NEXRAD) imagery, which are contained in NetCDF files. It can identify intersections between segmented NCFR cores and watersheds, in the form of Shapely polygons, and calculate the proportion of the watershed intersected. Results (propagation statistics) for each NCFR event are outputted via a csv file for every timestep--which displays: 
<ul>
  <li>The date and time of each timestep</li>
  <li>Whether an NCFR core intersects with the watershed</li>
  <li>The proportion of the watershed intersected (if above is true)</li>
  <li>The maximum reflectivity found in any of the segmented NCFR cores</li>
</ul>

For NCFR events where the algorithm correctly identifies the NCFR cores for multiple consecutive timesteps, you can obtain the csv file for the entire event--which displays for a centroid of an NCFR core:
<ul>
  <li>The start and end timestep</li>
  <li>The start and end position (lat, lon)</li>
  <li>The total distance traveled (km)</li>
  <li>The direction or azimuth of propagation (degrees)</li>
  <li>The average speed traveled (m/s)</li>
</ul>
