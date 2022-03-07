###############################################################################################################
# This driver function analyses both LiDAR data and field inventory data to produce independent estimates of
# canopy structure.  These are compared against each other and their integrated LAD is compared against LAI
# estimates from hemispherical photographs.
###############################################################################################################
import sys
import numpy as np
sys.path.append('../')
import LiDAR_io as io
import LiDAR_tools as lidar
import auxilliary_functions as aux
import LiDAR_MacHorn_LAD_profiles as LAD1

#------------------------------------------------------------------------------------
# DIRECTORIES
# start by defining input files
las_file = '../Carbon_plot_point_cloud_buffer.las'
subplot_coordinate_file = '../BALI_subplot_coordinates_corrected.csv'

#------------------------------------------------------------------------------------
# PARAMETERS
# define important parameters for canopy profile estimation
Plots = ['E','Belian','B North']
N_plots = len(Plots)
max_height = 80
layer_thickness = 1.
n_layers = np.ceil(max_height/layer_thickness)
minimum_height = 2.
subplot_width=20.
#kappa = 0.70
kappa=0.50
heights = np.arange(0,max_height,layer_thickness)+layer_thickness

#------------------------------------------------------------------------------------
# DECLARATIONS
# define dictionaries to host the various canopy profiles and LAI estimates that will be produced
MacArthurHorn_LAD = {}

#------------------------------------------------------------------------------------
# LOADING DATA
# load coordinates and lidar points for target areas
subplot_polygons, subplot_labels = aux.load_boundaries(subplot_coordinate_file)
all_lidar_pts = io.load_lidar_data(las_file)
#------------------------------------------------------------------------------------
# MAIN ANALYSIS
# LiDAR PROFILES LOOP
# loop through all plots to be analysed
for pp in range(0,N_plots):
    print(Plots[pp])
    Plot_name=Plots[pp]
    n_subplots = subplot_polygons[Plot_name].shape[0]

    #------------------------------------------------------------------------------------
    # CLIP DATA TO PLOT
    # clip LiDAR point cloud to plot level (this makes subsequent processing much faster)
    n_coord_pairs = subplot_polygons[Plot_name].shape[0]*subplot_polygons[Plot_name].shape[1]
    coord_pairs = subplot_polygons[Plot_name].reshape(n_coord_pairs,2)
    bbox_polygon = aux.get_bounding_box(coord_pairs)
    plot_lidar_pts = lidar.filter_lidar_data_by_polygon(all_lidar_pts,bbox_polygon,filter_by_first_return_location=True)
    starting_ids, trees = io.create_KDTree(plot_lidar_pts) # build kd-tree for plot lidar points
    print("canopy height = ", np.percentile(plot_lidar_pts[plot_lidar_pts[:,3]==1,2],99), "m")

    #------------------------------------------------------------------------------------
    # SET UP ARRAYS TO HOST RESULTS
    # get some subplot-level information
    n_subplots = subplot_polygons[Plot_name].shape[0]

    # set up some arrays to host the MacArthur-Horn profiles
    LAD_MH = np.zeros((n_subplots, heights.size))

    #------------------------------------------------------------------------------------
    # LOOP THROUGH SUBPLOTS, CALCULATING CANOPY PROFILES
    pt_count = 0.
    # loop through subplots, calculating both return profiles and LAD distributions
    for i in range(0,n_subplots):
        subplot_index = int(subplot_labels[Plot_name][i]-1)
        # filter lidar points into subplot
        subplot_poly = subplot_polygons[Plot_name][i,:,:]
        sp_pts = lidar.filter_lidar_data_by_polygon(plot_lidar_pts,subplot_poly,
                                                    filter_by_first_return_location=True)
        pt_count += sp_pts.shape[0]

        # now get MacArthur-Horn profiles
        """
        heights,first_return_profile,n_ground_returns = LAD1.bin_returns(sp_pts, max_height, layer_thickness)
        LAD_MH[subplot_index,:] = LAD1.estimate_LAD_MacArthurHorn(first_return_profile,
                                                                    n_ground_returns,
                                                                    layer_thickness,
                                                                    kappa,
                                                                    zero_nodata=False)
        """
        heights,weighted_return_profile,weighted_n_ground_returns = LAD1.bin_returns_weighted_by_num_returns(sp_pts, max_height, layer_thickness)
        LAD_MH[subplot_index,:] = LAD1.estimate_LAD_MacArthurHorn(weighted_return_profile,
                                                                    weighted_n_ground_returns,
                                                                    layer_thickness,
                                                                    kappa,
                                                                    zero_nodata=False)
        
        # Check for columns for which no pulses hit ground without interception.
        # Below the depth at which the canopy is not penetrated by first returns
        # some methods are infinite. Thus we need to expand the search radius
        # iteratively, so that we can build up a full vertical profile. Note that
        # this potentially gives rise to coarsening effective resolution down the
        # profile, but this seems preferable to more crude gap-filling schemes.
        nodata_test = np.any(~np.isfinite(LAD_MH[subplot_index]))
        centre_x = np.mean(subplot_poly[0:4,0])
        centre_y = np.mean(subplot_poly[0:4,1])
        radius=subplot_width/2.

        while nodata_test:
            # expand neighbourhood for point cloud sample
            ids = trees[0].query_ball_point([centre_x,centre_y], radius)
            sp_pts_iter = plot_lidar_pts[ids]

            # get MacArthur-Horn profiles
            nodata_gaps = ~np.isfinite(LAD_MH[subplot_index])
            """
            heights,first_return_profile,n_ground_returns = LAD1.bin_returns(sp_pts_iter, max_height, layer_thickness)
            LAD_MH[subplot_index,nodata_gaps] = LAD1.estimate_LAD_MacArthurHorn(first_return_profile,
                                                                    n_ground_returns,
                                                                    layer_thickness,
                                                                    kappa,
                                                                    zero_nodata=False)[nodata_gaps]
            """
            heights,weighted_return_profile,weighted_n_ground_returns = LAD1.bin_returns_weighted_by_num_returns(sp_pts_iter, max_height, layer_thickness)
            LAD_MH[subplot_index,nodata_gaps] = LAD1.estimate_LAD_MacArthurHorn(weighted_return_profile,
                                                                    weighted_n_ground_returns,
                                                                    layer_thickness,
                                                                    kappa,
                                                                    zero_nodata=False)[nodata_gaps]

            # update check
            radius+=1.
            nodata_test = np.any(~np.isfinite(LAD_MH[subplot_index]))

    print("average point density = ", pt_count/10.**4, " pts/m^2")

    #------------------------------------------------------------------------------------
    # CLEANING AND STORING
    # now we have looped through and created the different profiles, need to account for any NaN's and apply minimum height
    # to the LAD distributions
    # - remove all profile values below minimum height prior to comparison
    mask = heights <= minimum_height
    LAD_MH[:,mask]=np.nan
    MacArthurHorn_LAD[Plot_name] = LAD_MH.copy()
#----------------------------------------------------------------------------
np.savez('lidar_canopy_profiles_adaptive_for_synthesis.npz',(MacArthurHorn_LAD))
