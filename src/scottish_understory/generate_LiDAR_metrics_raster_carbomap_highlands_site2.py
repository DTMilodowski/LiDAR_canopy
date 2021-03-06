#===============================================================================
# This script creates raster datasets for specified LiDAR canopy metrics.      #
# Initially focusing on PAI, quantified using the MacArthur-Horn method.       #
# This is different to the original code as it uses the pixel boundaries to    #
# clip the point cloud rather than a pixel-centred circular neighbourhood      #
# The outline for the procedure is as follows:                                 #
# (1) Find the bounding box encompassing all the available LiDAR tiles This    #
#     will provide the raster extent                                           #
#     - alternatively one can always manually delimit a bounding box           #
# (2) Create the x and y coordinates of the grid for which the raster will be  #
#     produced, based on bounding box and desired resolution                   #
# (3) Loop through the las tiles in turn                                       #
#     (i)   load in with a 20 m buffer                                         #
#     (ii)  for bbox of tile                                                   #
#           - loop through x and y coord pairs within bounding box             #
#           - at each point sample point cloud and calculate metric            #
#             - sample with circular radius of and calculate metric            #
#           - (PAI, point density)                                             #
# (4) Write raster to geoTIFF                                                  #
#===============================================================================
import numpy as np
import sys
sys.path.append('../')
import os
import laspy as las
import LiDAR_io as io
import LiDAR_tools as lidar
import auxilliary_functions as aux
import LiDAR_MacHorn_LAD_profiles as pad
import structural_metrics as struct
import raster_io as raster

#-------------------------------------------------------------------------------
# Input files
las_list = '../las_lists/carbomap_site2_lastiles.txt'
savedir = '/home/dmilodow/DataStore_DTM/FOREST2020/LiDAR/carbomap_highlands/canopy_metrics/'
laz_files = False ## CHANGE AS REQUIRED - note that if laszip is installed, this
                  ## can safely be set to False for both las and laz files. If
                  ## True, this leads to a workaround for cases where laszip is
                  ## not available, in which las2las is called to write a temp
                  ## las file. This process is slow, so the recommendation is to
                  ## install laszip!

# Site ID
site = 'carbomap_site2'

# load in a raster to get coordinate system info for later
temp,temp_geoT,coord_sys = raster.load_GeoTIFF_band_and_georeferencing('/home/dmilodow/DataStore_DTM/FOREST2020/LiDAR/carbomap_highlands/carbomap_layers/Site2__Canopy_Height_Model_UTM.tif')

# Some parameters
radius = 10.#np.sqrt(10) # for buffer
max_height = 32.   # max tree height
min_height = 1.     # minimum height for inclusion in PAI
layer_thickness = 1 # thickness of vertical strata
heights = np.arange(0,max_height,layer_thickness)+layer_thickness
layers = heights.size
kappa = 1. # correction factor (clumping and leaf angles)
raster_res = 5. # resolution of output rasters

#-------------------------------------------------------------------------------
# Phase one - get bounding box of all las tiles.
# define bbox for region 3 (S)

N_ = 6207200
S_ = 6206000
W_ = 345000
E_ = 346400
"""
N_ = 6206500
S_ = 6206400
W_ = 345580
E_ = 346680
"""


bbox  = np.array([[W_,N_],[E_,N_],[E_,S_],[W_,S_]])

x_coords = np.arange(W_+raster_res/2.,E_+raster_res,raster_res)
y_coords = np.arange(S_+raster_res/2.,N_+raster_res,raster_res)

rows = y_coords.size
cols = x_coords.size
rows_ii = np.arange(y_coords.size)
cols_jj = np.arange(x_coords.size)

PAI = np.zeros((rows,cols))*np.nan
PAI_1_2m = np.zeros((rows,cols))*np.nan
PAI_2_5m = np.zeros((rows,cols))*np.nan
PAI_5m_up = np.zeros((rows,cols))*np.nan
can_height = np.zeros((rows,cols))*np.nan
PAD = np.zeros((rows,cols,layers))*np.nan
pulse_dens = np.zeros((rows,cols))
n_ground = np.zeros((rows,cols))

# Phase three - loop through las tiles and gradually fill the array
laz_files = io.find_las_files_by_polygon(las_list,bbox)
n_files = len(laz_files)
for i in range(0,n_files):
    print("Processing tile %i of %i" % (i+1,n_files))
    # get bbox of specific tile
    lasFile = las.file.File(laz_files[i],mode='r-')
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
    lidar_pts, starting_ids_for_trees, trees = io.load_lidar_data_by_polygon(las_list,
                                    polygon,laz_files=False,max_pts_per_tree = 5*10**5)
    N_trees = len(trees)

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

            # get pixel boundaries
            pixel_bbox = np.array([[x_iter-raster_res/2., y_iter-raster_res/2.], [x_iter-raster_res/2., y_iter+raster_res/2.],
                                   [x_iter+raster_res/2., y_iter+raster_res/2.], [x_iter+raster_res/2., y_iter-raster_res/2.],
                                   [x_iter-raster_res/2., y_iter-raster_res/2.]])

            # retrieve point clouds samples
            sample_pts = np.array([])
            for tt in range(0,N_trees):
                ids = trees[tt].query_ball_point([x_iter,y_iter], radius)
                if len(ids)>0:
                    if sample_pts.size==0:
                        sample_pts = lidar.filter_lidar_data_by_polygon(lidar_pts[np.asarray(ids)+starting_ids_for_trees[tt]],pixel_bbox)
                    else:
                        sample_iter = lidar.filter_lidar_data_by_polygon(lidar_pts[np.asarray(ids)+starting_ids_for_trees[tt]],pixel_bbox)
                        sample_pts = np.concatenate((sample_pts,sample_iter),axis=0)
                        sample_iter = None

            # If we have the returns, then calculate metric of interest - in
            # this case the PAI
            if sample_pts.size > 0:
                if np.sum(sample_pts[:,3]==1) > 0:
                    # calculate PAD profile
                    heights,first_return_profile,n_ground_returns = pad.bin_returns(sample_pts, max_height, layer_thickness)
                    PADprof = pad.estimate_LAD_MacArthurHorn(first_return_profile, n_ground_returns, layer_thickness, kappa,zero_nodata=False)

                    # vertically distributed PAI
                    PAD[row_ii,col_jj,:] = PADprof.copy()

                    # remove lowermost portion of profile
                    PAD_iter = PADprof.copy()
                    PAD_iter[heights<min_height]=0
                    PAI[row_ii,col_jj] = np.nansum(PAD_iter)

                    PAI_1_2m[row_ii,col_jj] = np.nansum(PAD_iter[np.all((heights>1,heights<=2),axis=0)])
                    PAI_2_5m[row_ii,col_jj] = np.nansum(PAD_iter[np.all((heights>2,heights<=5),axis=0)])
                    PAI_5m_up[row_ii,col_jj] = np.nansum(PAD_iter[heights>5])

                    can_height[row_ii,col_jj] = np.percentile(sample_pts[:,2],99)

                    # other metrics
                    pulse_dens[row_ii,col_jj] = np.sum(sample_pts[:,3]==1)/(raster_res**2.)
                    n_ground[row_ii,col_jj]= n_ground_returns

            sample_pts=None
    lidar_pts = None
    trees = None
    starting_ids_for_trees = None

np.savez('%s%s_metrics_5m' % (savedir,site),pulse_density=pulse_dens,pai=PAI,
                pad=PAD, n_ground=n_ground, pai_1_2m=PAI_1_2m, pai_2_5m=PAI_2_5m, pai_5m_up=PAI_5m_up,
                canopy_height=can_height)

metrics = np.load('%s%s_metrics_5m.npz' % (savedir,site))

# Now that the raster is filled, just need to write it to file
XMinimum = x_coords.min() - raster_res/2.
YMaximum = y_coords.max() + raster_res/2.
geoTransform = [ XMinimum, raster_res, 0, YMaximum, 0, -raster_res ]
metric_keys = list(metrics.keys())
for kk in range(0,len(metric_keys)):
    var = metric_keys[kk]
    print("\t\t\t Saving rasters: %s" % var)
    metrics[var][np.isnan(metrics[var])]=-9999
    raster.write_array_to_GeoTiff_with_coordinate_system(metrics[var], geoTransform, coord_sys,'%s%s_5m_%s.tif' % (savedir,site,var),north_up=False)
