# This is a test function to develop scheme to transform tiled lidar products
# into a continuous raster

import numpy as np
import sys
import os
import laspy as las
import LiDAR_io as io
import LiDAR_tools as lidar
import auxilliary_functions as aux
import LiDAR_MacHorn_LAD_profiles as PAD
import structural_metrics as struct
import raster_io as raster
from matplotlib import pyplot as plt

#-------------------------------------------------------------------------------
# Input files
las_list = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/SAFE_las_files/las_list_full_path.txt' ## CHANGE AS REQUIRED
laz_files = False ## CHANGE AS REQUIRED

# Some parameters
radius = 10.
raster_res=25.

#-------------------------------------------------------------------------------
# Phase one - get bounding box of all las tiles.
UR, LR, UL, LL = io.get_bbox_of_multiple_tiles(las_list)

# Phase two - construct host arrays
xmin = LL[0]
xmax = LR[0]
ymin = LL[1]
ymax = UL[1]

LL=None; LR=None; UL=None; UR=None

x_coords = np.arange(xmin+raster_res/2.,xmax,raster_res)
y_coords = np.arange(ymin+raster_res/2.,ymax,raster_res)

rows = y_coords.size
cols = x_coords.size
rows_ii = np.arange(y_coords.size)
cols_jj = np.arange(x_coords.size)

test = np.zeros((rows,cols))*np.nan


# Phase three - loop through las tiles and gradually fill the array
las_files = np.genfromtxt(las_list,delimiter=',',dtype='S256')
n_files = las_files.size
for i in range(0,n_files):
    # get bbox of specific tile
    lasFile = las.file.File(las_files[i],mode='r')
    max_xyz = lasFile.header.max
    min_xyz = lasFile.header.min
    lasFile.close()

    # buffer this bounding box with the search radius
    E = max_xyz[0]+radius
    N = max_xyz[1]+radius
    W = min_xyz[0]-radius
    S = min_xyz[1]-radius

    # Read in LiDAR points for region of interest
    polygon = np.asarray([[W,N],[E,N],[E,S],[W,S]])
    
    # Get the positions and indices of grid points within the bounds of the tile
    row_mask = np.all((y_coords>=min_xyz[1],y_coords<max_xyz[1]),axis=0)
    col_mask = np.all((x_coords>=min_xyz[0],x_coords<max_xyz[0]),axis=0)
    x_coords_tile = x_coords[col_mask]
    y_coords_tile = y_coords[row_mask]
    rows_tile = rows_ii[row_mask]
    cols_tile = cols_jj[col_mask]

    # Loop through relevant portion of the grid, sampling at each grid square
    for ii in range(0,rows_tile.size):
        y_iter = y_coords_tile[ii]
        row_ii = rows_tile[ii]
        for jj in range(0,cols_tile.size):
            x_iter = x_coords_tile[jj]
            col_jj = cols_tile[jj]

            # mark point with iteration number   
            test[row_ii,col_jj] = i

plt.imshow(test,origin='lower');plt.show()
