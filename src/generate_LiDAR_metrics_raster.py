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

# Some parameters
min_PAD = 0.1
radius = 10.
max_height = 80.   
min_height = 2.     
layer_thickness = 1
heights = np.arange(0,max_height,layer_thickness)+layer_thickness
kappa = 0.7
raster_res = 5

#-------------------------------------------------------------------------------
# Phase one - get bounding box of all las tiles.
UR, LR, UL, LL = io.get_bbox_of_multiple_tiles(las_list)

# Phase two - construct host arrays
xmin = LL[0]
xmax = LR[0]
ymin = LL[1]
ymax = UL[1]

x_coords = np.arange(xmin+raster_res/2.,xmax,raster_res)
y_coords = np.arange(ymin+raster_res/2.,ymax,raster_res)

rows = y_coords.size
cols = x_coords.size
rows = np.arange(y_coords.size)
cols = np.arange(x_coords.size)

PAI = np.zeros((rows,cols))*np.nan
pt_density = np.zeros((rows,cols))

# Phase three - loop through las tiles and gradually fill the array
las_files = np.genfromtxt(file_list,delimiter=',',dtype='S256')
n_files = las_files.size
for i in range(0,n_files):
    lasFile = las.file.File(las_files[i],mode='r')
    max_xyz = lasFile.header.max
    min_xyz = lasFile.header.min
    lasFile.close()
