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
#las_list = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/SAFE_las_files/las_list_full_path.txt' ## CHANGE AS REQUIRED
#las_list = 'maliau_las_list.txt' ## CHANGE AS REQUIRED
#las_list = 'danum_las_list.txt' ## CHANGE AS REQUIRED
#las_list = 'sepilok_laz_list.txt' ## CHANGE AS REQUIRED
las_list = 'gola_nomo_area_04_lazfiles.txt'
#las_list = 'gola_lalehun_lazfiles.txt'
#las_list = 'gola_transects_south_lazfiles.txt'
#las_list = 'gola_transects_central_lazfiles.txt'
#las_list = 'gola_central_area_01_02_lazfiles.txt'
#las_list = 'gola_mayengeima_area_03_lazfiles.txt'
laz_files = True ## CHANGE AS REQUIRED

# Site ID
#site = 'SAFE'
#site = 'maliau'
#site = 'danum'
#site = 'sepilok'
site = 'gola_nomo_area_04'
#site = 'gola_lalehun'
#site = 'gola_transects_south'
#site = 'gola_transects_central'
#site = 'gola_central_area_01_02'
#site = 'gola_mayengeima_area_03'


# Some parameters
min_PAD = 0.1
radius = np.sqrt(200)
max_height = 80.   
min_height = 2.     
layer_thickness = 1
heights = np.arange(0,max_height,layer_thickness)+layer_thickness
kappa = 0.7
raster_res = 20

# Some georeferencing info
#utm = 50
utm = 29
#EPSG = "32650" # WGS84 / UTM 50N

#-------------------------------------------------------------------------------
# Phase one - get bounding box of all las tiles.
UR, LR, UL, LL = io.get_bbox_of_multiple_tiles(las_list,laz_files)

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
PAI10 = np.zeros((rows,cols))*np.nan
PAI20 = np.zeros((rows,cols))*np.nan
PAI30 = np.zeros((rows,cols))*np.nan
PAI40 = np.zeros((rows,cols))*np.nan
PAI50 = np.zeros((rows,cols))*np.nan
PAI60 = np.zeros((rows,cols))*np.nan
PAI70 = np.zeros((rows,cols))*np.nan
PAI80 = np.zeros((rows,cols))*np.nan
Shannon = np.zeros((rows,cols))*np.nan
Shape = np.zeros((rows,cols))*np.nan
layers = np.zeros((rows,cols))*np.nan
pt_dens = np.zeros((rows,cols))
can_ht = np.zeros((rows,cols))*np.nan
mean = np.zeros((rows,cols))*np.nan
std = np.zeros((rows,cols))*np.nan
skew = np.zeros((rows,cols))*np.nan
kurt = np.zeros((rows,cols))*np.nan


# Phase three - loop through las tiles and gradually fill the array
las_files = np.genfromtxt(las_list,delimiter=',',dtype='S256')
n_files = las_files.size
for i in range(0,n_files):
    print "Processing tile %i of %i" % (i+1,n_files)
    # get bbox of specific tile
    if laz_files:
        os.system("las2las %s temp.las" % las_files[i])
        lasFile = las.file.File('temp.las',mode='r')
        max_xyz = lasFile.header.max
        min_xyz = lasFile.header.min
        lasFile.close()
        os.system("rm temp.las")
    else:
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
                        #sample_pts = lidar_pts[np.asarray(ids)+starting_ids_for_trees[tt]]
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
                
                    PAI[row_ii,col_jj] = np.sum(PAD_iter)

                    # vertically distributed PAI
                    PAI10[row_ii,col_jj] = np.sum(PAD_iter[np.all((heights>2,heights<=10),axis=0)])
                    PAI20[row_ii,col_jj] = np.sum(PAD_iter[np.all((heights>10,heights<=20),axis=0)])
                    PAI30[row_ii,col_jj] = np.sum(PAD_iter[np.all((heights>20,heights<=30),axis=0)])
                    PAI40[row_ii,col_jj] = np.sum(PAD_iter[np.all((heights>30,heights<=40),axis=0)])
                    PAI50[row_ii,col_jj] = np.sum(PAD_iter[np.all((heights>40,heights<=50),axis=0)])
                    PAI60[row_ii,col_jj] = np.sum(PAD_iter[np.all((heights>50,heights<=60),axis=0)])
                    PAI70[row_ii,col_jj] = np.sum(PAD_iter[np.all((heights>60,heights<=70),axis=0)])
                    PAI80[row_ii,col_jj] = np.sum(PAD_iter[np.all((heights>70,heights<=80),axis=0)])

                    # other metrics
                    pt_dens[row_ii,col_jj] = sample_pts.shape[0]/(np.pi*radius**2.)
                    Shannon[row_ii,col_jj] = struct.calculate_Shannon_index(PAD_iter)
                    Shape[row_ii,col_jj] = struct.calculate_Shannon_index(PAD_iter)
                    layers[row_ii,col_jj] = struct.calculate_number_of_contiguous_layers(heights,PAD_iter,min_PAD)
                    can_ht[row_ii,col_jj] = np.percentile(sample_pts[sample_pts[:,3]==1,2],99)
                    if np.sum(PAD_iter)>0:
                        mean[row_ii,col_jj],std[row_ii,col_jj],skew[row_ii,col_jj],kurt[row_ii,col_jj] = struct.calculate_moments_of_distribution(heights,PAD_iter) # note this function will kick out infs if there is only one layer occupied

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
