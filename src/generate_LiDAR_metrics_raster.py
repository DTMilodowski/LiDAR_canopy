#===============================================================================
# This script creates raster datasets for specified LiDAR canopy metrics.      #
# Initially focusing on PAI, quantified using the MacArthur-Horn method.       #
# The outline for the procedure is as follows:                                 #
# (1) Find the bounding box encompassing all the available LiDAR tiles This    #
#     will provide the raster extent                                           #
#     - alternatively one can always manually delimit a bounding box           #
# (2) Create the x and y coordinates of the grid for which the raster will be  #
#     produced, based on bounding box and desired resolution                   #
# (3) Loop through the las tiles in turn                                       #
#     (i)   load in with a 10 m buffer                                         #
#     (ii)  for bbox of tile                                                   #
#           - loop through x and y coord pairs within bounding box             #
#           - at each point sample point cloud and calculate metric            #
#           - (PAI, point density)                                             #
# (4) Write raster to geoTIFF                                                  #
#===============================================================================
import numpy as np
import sys
import os
import LiDAR_io as io
import LiDAR_tools as lidar
import auxilliary_functions as aux
import LiDAR_MacHorn_LAD_profiles as PAD
import structural_metrics as struct

#-------------------------------------------------------------------------------
# Input files
las_list = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/SAFE_las_files/las_list_full_path.txt' ## CHANGE AS REQUIRED
laz_files = False ## CHANGE AS REQUIRED

# Output files
PAI_raster = '/exports/csce/datastore/geos/users/dmilodow/BALI/LiDAR/Data/RasterData/SAFE_PAI.tif'
dens_raster = '/exports/csce/datastore/geos/users/dmilodow/BALI/LiDAR/Data/RasterData/SAFE_point_density.tif'


