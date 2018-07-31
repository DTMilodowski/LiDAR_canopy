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
import os
import laspy as las
import LiDAR_io as io
import LiDAR_tools as lidar
import auxilliary_functions as aux
import LiDAR_MacHorn_LAD_profiles as PAD
import structural_metrics as struct
import raster_io as raster

#-------------------------------------------------------------------------------
# Input files
laz_list = 'thetford2015_lazfiles.txt'
laz_files = True ## CHANGE AS REQUIRED

# Site ID
site = 'Thetford2015'


# Some parameters
min_PAD = 0.1
radius = np.sqrt(200)
max_height = 25.   
min_height = 2.     
layer_thickness = 1
heights = np.arange(0,max_height,layer_thickness)+layer_thickness
layers = heights.size
kappa = 1.
raster_res = 5

# Some georeferencing info
utm = 31

#-------------------------------------------------------------------------------
# Phase one - get bounding box of all las tiles.
UR, LR, UL, LL = io.get_bbox_of_multiple_tiles(laz_list,laz_files)

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

PAI = np.zeros((rows,cols))*np.nan
PAD = np.zeros((layers,rows,cols))*np.nan
Shannon = np.zeros((rows,cols))*np.nan
pt_dens = np.zeros((rows,cols))


# Phase three - loop through las tiles and gradually fill the array
laz_files = np.genfromtxt(laz_list,delimiter=',',dtype='S256')
n_files = laz_files.size
for i in range(0,n_files):
    print "Processing tile %i of %i" % (i+1,n_files)
    # get bbox of specific tile
    if laz_files:
        os.system("las2las %s temp.las" % laz_files[i])
        lasFile = las.file.File('temp.las',mode='r')
        max_xyz = lasFile.header.max
        min_xyz = lasFile.header.min
        lasFile.close()
        os.system("rm temp.las")
    else:
        lasFile = las.file.File(laz_files[i],mode='r')
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
    lidar_pts, starting_ids_for_trees, trees = io.load_lidar_data_by_polygon(las_list,polygon,max_pts_per_tree = 5*10**5, laz_files=laz_files)
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
                    heights,first_return_profile,n_ground_returns = PAD.bin_returns(sample_pts, max_height, layer_thickness)
                    PADprof = PAD.estimate_LAD_MacArthurHorn(first_return_profile, n_ground_returns, layer_thickness, kappa)
                    
                    # remove lowermost portion of profile
                    PAD_iter = PADprof.copy()
                    PAD_iter[heights<min_height]=0
                
                    PAI[row_ii,col_jj] = np.nansum(PAD_iter)

                    # vertically distributed PAI
                    PAD[:,row_ii,col_jj] = PAD_prof.copy()

                    # other metrics
                    pt_dens[row_ii,col_jj] = sample_pts.shape[0]/(np.pi*radius**2.)
                    Shannon[row_ii,col_jj] = struct.calculate_Shannon_index(PAD_iter[np.isnan(PAD_iter))
            sample_pts=None
    lidar_pts = None
    trees = None
    starting_ids_for_trees = None

np.savez('%s_metrics_10m' % site,point_density=pt_dens,pai=PAI,shannon=Shannon,shape=Shape,pai_02_10m=PAI10,pai_10_20m=PAI20,pai_20_30m=PAI30,pai_30_40m=PAI40,pai_40_50m=PAI50,pai_50_60m=PAI60,pai_60_70m=PAI70,pai_70_80m=PAI80,n_layers=layers, canopy_height = can_ht, mean=mean, std=std, skew=skew,kurt=kurt)

metrics = np.load('%s_metrics_10m.npz' % site)

# Now that the raster is filled, just need to write it to file
XMinimum = x_coords.min() - raster_res/2.
YMaximum = y_coords.max() + raster_res/2.
geoTransform = [ XMinimum, raster_res, 0, YMaximum, 0, -raster_res ]

for kk in range(0,len(metrics.keys())):
    var = metrics.keys()[kk]
    print "\t\t\t Saving rasters: %s" % var
    raster.write_raster_to_GeoTiff_UTM(metrics[var], geoTransform, ('%s_pointcloud_metrics_20m_%s' % (site,var)), utm)

"""
import sys
sys.path.append('/home/dmilodow/DataStore_DTM/FOREST2020/EOdata/EO_data_processing/src/plot_EO_data/colormap/')
import colormaps as cmaps
plt.register_cmap(name='inferno', cmap=cmaps.inferno)
plt.set_cmap(cmaps.inferno)
"""