###############################################################################################################
# This driver function analyses both LiDAR data and field inventory data to produce independent estimates of
# canopy structure.  These are compared against each other and their integrated LAD is compared against LAI
# estimates from hemispherical photographs.
###############################################################################################################
import numpy as np
import LiDAR_io as io
import LiDAR_tools as lidar
import auxilliary_functions as aux
import LiDAR_radiative_transfer_LAD_profiles as PAD

#------------------------------------------------------------------------------------
# DIRECTORIES
# start by defining input files
las_file = 'Carbon_plot_point_cloud_buffer.las'
subplot_coordinate_file = 'BALI_subplot_coordinates_corrected.csv'
gps_pts_file = 'GPS_points_file_for_least_squares_fitting.csv'

# also define output directory (for saving data)
output_dir = '/home/dmilodow/DataStore_DTM/BALI/PAPERS/PaperDrafts/EstimatingCanopyStructureBALI/Profiles/'

#------------------------------------------------------------------------------------
# PARAMETERS
# define important parameters for canopy profile estimation
Plots = [b'Belian',b'E',b'B South']
#Plots = ['B North']
N_plots = len(Plots)
max_height = 80
max_return = 3
layer_thickness = 1.
minimum_height = 2.
plot_area = 10.**4
subplot_area = 20.*20.

heights_rad = np.arange(0,max_height+layer_thickness,layer_thickness)

#------------------------------------------------------------------------------------
# DECLARATIONS
# define dictionaries to host the various canopy profiles and LAI estimates that will be produced
planophile_PAD = {}
erectophile_PAD = {}
spherical_PAD = {}
planophile_PAD_mean = {}
erectophile_PAD_mean = {}
spherical_PAD_mean = {}
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

    #------------------------------------------------------------------------------------
    # SET UP ARRAYS TO HOST RESULTS
    # get some subplot-level information
    n_subplots = subplot_polygons[Plot_name].shape[0]
    for ss in range(0,n_subplots):
        subplot_index = int(subplot_labels[Plot_name][ss]-1)
        # filter lidar points into subplot
        sp_pts = lidar.filter_lidar_data_by_polygon(plot_lidar_pts,subplot_polygons[Plot_name][ss,:,:],filter_by_first_return_location=True)

    # set up some arrays to host the radiative transfer based profiles
    PAD_pla = np.zeros((n_subplots,heights_rad.size,max_return))
    PAD_ere = np.zeros((n_subplots,heights_rad.size,max_return))
    PAD_sph = np.zeros((n_subplots,heights_rad.size,max_return))

    #------------------------------------------------------------------------------------
    # LOOP THROUGH SUBPLOTS, CALCULATING CANOPY PROFILES
    pt_count = 0.
    # loop through subplots, calculating both return profiles and LAD distributions
    for i in range(0,n_subplots):
        subplot_index = int(subplot_labels[Plot_name][i]-1)
        # filter lidar points into subplot
        sp_pts = lidar.filter_lidar_data_by_polygon(plot_lidar_pts,subplot_polygons[Plot_name][i,:,:],filter_by_first_return_location=True)
        pt_count += sp_pts.shape[0]

        # first loop through the return numbers to calculate the radiative LAD profiles
        for rr in range(0,max_return):
            max_k=rr+1
            u,n,I,U = PAD.calculate_LAD_DTM(sp_pts,heights_rad,max_k,'spherical')
            PAD_sph[subplot_index,:,rr]=u.copy()
            u,n,I,U = PAD.calculate_LAD_DTM(sp_pts,heights_rad,max_k,'erectophile')
            PAD_ere[subplot_index,:,rr]=u.copy()
            u,n,I,U = PAD.calculate_LAD_DTM(sp_pts,heights_rad,max_k,'planophile')
            PAD_pla[subplot_index,:,rr]=u.copy()

    #------------------------------------------------------------------------------------
    # AVERAGE 1 ha PROFILES
    PAD_sph_mean = np.nansum(PAD_sph,axis=0)/(np.sum(np.isfinite(PAD_sph),axis=0)).astype('float')
    PAD_pla_mean = np.nansum(PAD_pla,axis=0)/(np.sum(np.isfinite(PAD_pla),axis=0)).astype('float')
    PAD_ere_mean = np.nansum(PAD_ere,axis=0)/(np.sum(np.isfinite(PAD_ere),axis=0)).astype('float')

    # CLEANING AND STORING
    # now we have looped through and created the different profiles, need to account for any NaN's and apply minimum height
    # to the LAD distributions
    # - remove all profile values below minimum height prior to comparison
    mask = np.max(heights_rad)-heights_rad<=minimum_height
    PAD_sph[:,mask]=np.nan
    PAD_ere[:,mask]=np.nan
    PAD_pla[:,mask]=np.nan

    # now the profiles are ready to be stored into their relevant dictionaries
    spherical_PAD[Plot_name] = PAD_sph.copy()
    planophile_PAD[Plot_name] = PAD_pla.copy()
    erectophile_PAD[Plot_name] = PAD_ere.copy()
    spherical_PAD_mean[Plot_name] = PAD_sph_mean.copy()
    planophile_PAD_mean[Plot_name] = PAD_pla_mean.copy()
    erectophile_PAD_mean[Plot_name] = PAD_ere_mean.copy()
#----------------------------------------------------------------------------
# Plot profiles
figure_name = 'figS1_leaf_angle_distributions_PAD.png'
figure_number = 111
csp.plot_leaf_angle_distribution_profile_comparison(figure_name,figure_number,heights_rad,
                        spherical_PAD_mean,erectophile_PAD_mean,planophile_PAD_mean)
