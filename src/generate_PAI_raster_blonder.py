#===============================================================================
# This is an adaptation of the generate_PAI_raster.py script to carry out some #
# analysis for Ben Blonder as part of his microclimate work. For each site:    #
# (1) Find the bounding box encompassing the plot (using his geoTIFF)          #
# (2) Create the x and y coordinates of the grid for which the raster will be  #
#     produced, based on this bounding box and desired resolution (1 m)        #
# (3) Loop through the las tiles in turn                                       #
#     (i)   load in LiDAR data with 10 m buffer                                #
#     (ii)  loop through x and y coord pairs within grid                       #
#     (iii) at each point sample point cloud and calculate PAI and point       #
#           density                                                            #
# (4) Write rasters to geoTIFF                                                 #
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
las_list = ['/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/SAFE_las_files/maliau_las_list_full_path.txt',
            '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/SAFE_las_files/las_list_full_path.txt',
            '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/SAFE_las_files/las_list_full_path.txt']

# Site ID
site = ['maliaubelian','safee','bso']

# Some parameters
min_PAD = 0.1
radius = 5.
max_height = 80.   
min_height = 2.     
layer_thickness = 1
heights = np.arange(0,max_height,layer_thickness)+layer_thickness
kappa = 0.7

# Some georeferencing info
utm = 50

for ss in range(0,3):
    #-------------------------------------------------------------------------------
    # Phase one - get georeferencing information from geotiff
    template_geotiff = '../Data/bblonder/dem_%s.tif' % site[ss]
    template, geoTransform, coord_sys = raster.load_GeoTIFF_band_and_georeferencing(template_geotiff)
    rows,cols = template.shape
    
    xmin = geoTransform[0]
    xmax = xmin + float(cols)*geoTransform[1]

    # check that y resolution is negative (standard formatting)
    if geoTransform[5]<0:
        ymax = geoTransform[3]
        ymin = ymax + float(rows)*geoTransform[5]
    else:
        ymin = geoTransform[3]
        ymax = ymin + float(rows)*geoTransform[5]
    
    raster_res = geoTransform[1] # assume x and y resolutions are the same
    x_coords = np.arange(xmin+raster_res/2.,xmax,raster_res)
    y_coords = np.arange(ymin+raster_res/2.,ymax,raster_res)

    # sanity check
    if rows != y_coords.size:
        print "PROBLEM - rows and y_coords have different size"
    if cols != x_coords.size:
        print "PROBLEM - cols and x_coords have different size"

    PAI = np.zeros((rows,cols))*np.nan
    pt_dens = np.zeros((rows,cols))


    # Phase three - loop through las tiles and gradually fill the array
    #------------------------------------------------------------------
    # - first find the relevant las tiles
    # define a bounding box around target points to load in point cloud around area of interest
    # due to the way the bbox is queried, need to buffer this at the moment. at some point I will
    # get round to fixing this, but at least it works!
    E = xmax+500
    N = ymax+500
    W = xmin-500
    S = ymin-500
    
    # Read in LiDAR points for region of interest
    polygon = np.asarray([[W,N],[E,N],[E,S],[W,S]])
    lidar_pts, starting_ids_for_trees, trees = io.load_lidar_data_by_polygon(las_list[ss],polygon,max_pts_per_tree = 5*10**5)

    N_trees = len(trees)
    # Loop through the grid, sampling at each grid square
    for ii in range(0,rows):
        for jj in range(0,cols):
            # retrieve point clouds samples        
            sample_pts = np.array([])
            for tt in range(0,N_trees):
                ids = trees[tt].query_ball_point([x_coords[jj],y_coords[ii]], radius)
                if len(ids)>0:
                    if sample_pts.size==0:
                        sample_pts = lidar_pts[np.asarray(ids)+starting_ids_for_trees[tt]]
                    else:
                        sample_pts = np.concatenate((sample_pts,lidar_pts[np.asarray(ids)+starting_ids_for_trees[tt]]),axis=0)
        
            # If we have the returns, then calculate metric of interest - in this case the PAI
            if sample_pts.size > 0:
                if np.sum(sample_pts[:,3]==1) > 0:
                    # calculate PAD profile
                    heights,first_return_profile,n_ground_returns = PAD.bin_returns(sample_pts, max_height, layer_thickness)
                    PADprof = PAD.estimate_LAD_MacArthurHorn(first_return_profile, n_ground_returns, layer_thickness, kappa)
                    
                    # remove lowermost portion of profile
                    PAD_iter = PADprof.copy()
                    PAD_iter[heights<min_height]=0
                
                    PAI[ii,jj] = np.sum(PAD_iter)
                    pt_dens[ii,jj] = sample_pts.shape[0]/(np.pi*radius**2.)

            sample_pts=None

    lidar_pts = None
    trees = None
    starting_ids_for_trees = None

    # Now that the raster is filled, just need to write it to file
    print "\t\t\t Saving rasters"
    raster.write_raster_to_GeoTiff_UTM(PAI, geoTransform, ('%s_pointcloud_metrics_05m_PAI' % (site[ss])), utm)
    raster.write_raster_to_GeoTiff_UTM(pt_dens, geoTransform, ('%s_pointcloud_metrics_05m_point_density' % (site[ss])), utm)
