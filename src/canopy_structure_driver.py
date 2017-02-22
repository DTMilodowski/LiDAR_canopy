###############################################################################################################
# This driver function analyses both LiDAR data and field inventory data to produce independent estimates of
# canopy structure.  These are compared against each other and their integrated LAD is compared against LAI
# estimates from hemispherical photographs.
############################################################################################################### 
import numpy as np
import sys
from matplotlib import pyplot as plt
import LiDAR_tools as lidar
import auxilliary_functions as aux
import LiDAR_MacHorn_LAD_profiles as LAD1
import LiDAR_radiative_transfer_LAD_profiles as LAD2
import structural_metrics as structure

# start by defining input files
las_file = 'Carbon_plot_point_cloud_buffer.las'
subplot_coordinate_file = 'BALI_subplot_coordinates_corrected.csv'

# also define output directory (for saving figures)
output_dir = './Figures/'

# define important parameters for canopy profile estimation
Plots = ['LF','E','Belian','Seraya','B North','B South','DC1','DC2']
N_plots = len(Plots)
leaf_angle_dist = 'spherical'
max_height = 80
max_return = 3
layer_thickness = 1
n_layers = np.ceil(max_height/layer_thickness)
minimum_height = 2.

# define dictionaries to host the various canopy profiles and LAI estimates that will be produced
MacArthurHorn_LAD = {}
radiative_DTM_LAD = {}

# load coordinates and lidar points for target areas
subplot_polygons, subplot_labels = aux.load_boundaries(subplot_coordinate_file)
all_lidar_pts = lidar.load_lidar_data(las_file)

# loop through all plots to be analysed
for pp in range(0,N_plots):
    print Plots[pp]
    Plot_name=Plots[pp]
    # clip LiDAR point cloud to plot level (this makes subsequent processing much faster)
    n_coord_pairs = subplot_polygons[Plot_name].shape[0]*subplot_polygons[Plot_name].shape[1]
    coord_pairs = subplot_polygons[Plot_name].reshape(n_coord_pairs,2)
    bbox_polygon = aux.get_bounding_box(coord_pairs)
    plot_lidar_pts = lidar.filter_lidar_data_by_polygon(all_lidar_pts,bbox_polygon)
    
    # get some subplot-level information
    n_subplots = subplot_polygons[Plot_name].shape[0]

    # set up some arrays to host the radiative transfer based profiles
    heights_rad = np.arange(0,max_height+1)
    LAD_rad_DTM = np.zeros((n_subplots,heights_rad.size,max_return))

    # set up some arrays to host the MacArthur-Horn profiles
    heights = np.arange(0,max_height)+1
    LAD_MH = np.zeros((n_subplots, heights.size))

    # loop through subplots, calculating both return profiles and LAD distributions
    for i in range(0,n_subplots):
        print "Subplot: ", subplot_labels[Plot_name][i]
        # filter lidar points into subplot
        sp_pts = lidar.filter_lidar_data_by_polygon(plot_lidar_pts,subplot_polygons[Plot_name][i,:,:])

        # first of all, loop through the return numbers to calculate the radiative LAD profiles, accounting for imperfect penetration of LiDAR pulses into canopy
        for rr in range(0,max_return):
            max_k=rr+1
            u,n,I,U = LAD2.calculate_LAD_DTM(sp_pts,heights_rad,max_k,'spherical')
            LAD_rad_DTM[i,:,rr]=u.copy()

        # now get MacArthur-Horn profiles
        heights,first_return_profile,n_ground_returns = LAD1.bin_returns(sp_pts, max_height, layer_thickness)
        LAD_MH[i,:] = LAD1.estimate_LAD_MacArtherHorn(first_return_profile, n_ground_returns, layer_thickness, 1.)

    # now we have looped through and created the different profiles, need to account for any NaN's and apply minimum height
    # to the LAD distributions
    # - set NaN values to zero
    LAD_rad_DTM[np.isnan(LAD_rad_DTM)]=0
    LAD_MH[np.isnan(LAD_MH)]=0
    # - remove all profile values below minimum height prior to comparison
    mask = heights <= minimum_height
    LAD_MH[:,mask]=0
    mask = np.max(heights_rad)-heights_rad<=minimum_height
    LAD_rad_DTM[:,mask]=np.nan

    # Store arrays into respective dictionaries
    MacArthurHorn_LAD[Plot_name] = LAD_MH.copy()
    radiative_DTM_LAD[Plot_name] = LAD_rad_DTM.copy()


# next step is to take the canopy profiles and get the vertical and horizontal canopy structural metrics

# first up - create some dictionaries to host the structural metrics
frechet_dist_MH = {}
peak_heights_MH = {}
vertical_structural_variance_MH = {}
VSI_MH = {}
frechet_dist_rad = {}
peak_heights_rad = {}
vertical_structural_variance_rad = {}
VSI_rad = {}

for pp in range(0,N_plots):
    print Plots[pp]
    plot_name = Plots[pp]
    # get peaks
    print "\t- getting vertical metrics"
    peak_heights_MH[plot_name] = structure.retrieve_peaks(MacArthurHorn_LAD[Plot_name],heights)
    peak_heights_rad[plot_name] = structure.retrieve_peaks(radiative_DTM_LAD[Plot_name][:,:,-1],heights_rad)
    # get variance in layer heights
    vertical_structural_variance_MH[plot_name] = structure.calculate_vertical_structural_variance(peak_heights_MH[plot_name])
    vertical_structural_variance_rad[plot_name] = structure.calculate_vertical_structural_variance(peak_heights_rad[plot_name])
    # get VSI
    VSI_MH[plot_name] = structure.calculate_VSI(peak_heights_MH[plot_name])
    VSI_rad[plot_name] = structure.calculate_VSI(peak_heights_rad[plot_name])
    print "\t- getting horizontal metrics"
    # get mean Frechet distance
    frechet_dist_MH[plot_name] = structure.calculate_mean_Frechet_distance(MacArthurHorn_LAD[Plot_name],heights)
    frechet_dist_rad[plot_name] = structure.calculate_mean_Frechet_distance(radiative_DTM_LAD[Plot_name][:,:,-1],heights_rad)


