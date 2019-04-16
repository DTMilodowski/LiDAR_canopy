###############################################################################################################
# This driver function analyses both LiDAR data and field inventory data to produce independent estimates of
# canopy structure.  These are compared against each other and their integrated LAD is compared against LAI
# estimates from hemispherical photographs.
###############################################################################################################
import numpy as np
import sys
from matplotlib import pyplot as plt
import LiDAR_io as io
import LiDAR_tools as lidar
import auxilliary_functions as aux
import LiDAR_MacHorn_LAD_profiles as LAD1
import LiDAR_radiative_transfer_LAD_profiles as LAD2
import inventory_based_LAD_profiles as field
import least_squares_fitting as lstsq
from scipy import stats

from matplotlib import rcParams
from datetime import datetime

import load_field_data as cen

import canopy_structure_plots as csp

# Set up some basiic parameters for the plots
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = 8
rcParams['legend.numpoints'] = 1
axis_size = rcParams['font.size']+2

colour = ['#46E900','#1A2BCE','#E0007F']
rgb = [[70,233,0],[26,43,206],[224,0,127]]
labels = ['$1^{st}$', '$2^{nd}$', '$3^{rd}$', '$4^{th}$']

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
Plots = [b'LF',b'E',b'Belian',b'Seraya',b'B North',b'B South',b'DC1',b'DC2']
#Plots = ['B North']
N_plots = len(Plots)
leaf_angle_dist = 'spherical'
max_height = 80
max_return = 3
layer_thickness = 1.
n_layers = np.ceil(max_height/layer_thickness)
minimum_height = 2.
plot_area = 10.**4
subplot_area = 20.*20.
beta_min = 0.26 # 95% range reported for US tree species by Purves et al. (2007)
beta_max = 0.44 # 95% range reported for US tree species by Purves et al. (2007)
kappa = 0.70

heights = np.arange(0,max_height,layer_thickness)+layer_thickness
heights_rad = np.arange(0,max_height+layer_thickness,layer_thickness)

#------------------------------------------------------------------------------------
# DECLARATIONS
# define dictionaries to host the various canopy profiles and LAI estimates that will be produced
MacArthurHorn_LAD = {}
MacArthurHorn_LAD_mean = {}
radiative_LAD = {}
radiative_DTM_LAD = {}
radiative_DTM_LAD_old = {}
radiative_LAD_mean = {}
radiative_DTM_LAD_mean = {}
lidar_profiles ={}
lidar_profiles_adjusted ={}
penetration_limit ={}

MacArthurHorn_LAI = {}
radiative_LAI = {}
radiative_DTM_LAI = {}
MacArthurHorn_LAI_mean = {}
radiative_LAI_mean = {}
radiative_DTM_LAI_mean = {}

plot_point_cloud = {}
max_height_LiDAR={}
max_height_field={}
#------------------------------------------------------------------------------------
# LOADING DATA
# load coordinates and lidar points for target areas
subplot_polygons, subplot_labels = aux.load_boundaries(subplot_coordinate_file)
all_lidar_pts = io.load_lidar_data(las_file)

#------------------------------------------------------------------------------------
# MAIN ANALYSIS
# LiDAR PROFILES LOOP
# loop through all plots to be analysed
plot_point_cloud= np.load('plot_point_clouds.npz')['arr_0'][()]
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
    plot_point_cloud[Plots[pp]] = lidar.filter_lidar_data_by_polygon(all_lidar_pts,bbox_polygon,filter_by_first_return_location=True)
    plot_lidar_pts=plot_point_cloud[Plots[pp]]
    print("canopy height = ", np.percentile(plot_lidar_pts[plot_lidar_pts[:,3]==1,2],99), "m")

    #------------------------------------------------------------------------------------
    # SET UP ARRAYS TO HOST RESULTS
    # get some subplot-level information
    n_subplots = subplot_polygons[Plot_name].shape[0]
    max_height_LiDAR[Plot_name]=np.zeros(n_subplots)
    for ss in range(0,n_subplots):
        subplot_index = int(subplot_labels[Plot_name][ss]-1)
        # filter lidar points into subplot
        sp_pts = lidar.filter_lidar_data_by_polygon(plot_lidar_pts,subplot_polygons[Plot_name][ss,:,:],filter_by_first_return_location=True)
        max_height_LiDAR[Plot_name][ss]=np.max(sp_pts[:,2])

    # set up some arrays to host the radiative transfer based profiles
    LAD_rad = np.zeros((n_subplots,heights_rad.size,max_return))
    LAD_rad_DTM = np.zeros((n_subplots,heights_rad.size,max_return))
    LAD_rad_DTM_old = np.zeros((n_subplots,heights_rad.size,max_return))

    # set up some arrays to host the MacArthur-Horn profiles
    #heights = np.arange(0,max_height)+1
    LAD_MH = np.zeros((n_subplots, heights.size))

    # set up array to host the lidar return profiles
    lidar_return_profiles = np.zeros((n_subplots, heights_rad.size, max_return))
    lidar_return_profiles_adj = np.zeros((n_subplots, heights_rad.size, max_return))

    penetration_lim = np.zeros((n_subplots,heights.size))

    #------------------------------------------------------------------------------------
    # LOOP THROUGH SUBPLOTS, CALCULATING CANOPY PROFILES
    pt_count = 0.
    # loop through subplots, calculating both return profiles and LAD distributions
    for i in range(0,n_subplots):
        #print "Subplot: ", subplot_labels[Plot_name][i]
        subplot_index = int(subplot_labels[Plot_name][i]-1)
        # filter lidar points into subplot
        sp_pts = lidar.filter_lidar_data_by_polygon(plot_lidar_pts,subplot_polygons[Plot_name][i,:,:],filter_by_first_return_location=True)
        pt_count += sp_pts.shape[0]

        # first loop through the return numbers to calculate the radiative LAD profiles
        for rr in range(0,max_return):
            max_k=rr+1
            u,n,I,U = LAD2.calculate_LAD(sp_pts,heights_rad,max_k,'spherical')
            LAD_rad[subplot_index,:,rr]=u.copy()
        lidar_return_profiles[subplot_index,:,:n.shape[2]] = np.sum(n.copy(),axis=1)

        # now repeat but for adjusted profiles, accounting for imperfect penetration of LiDAR pulses into canopy
        for rr in range(0,max_return):
            max_k=rr+1
            u,n,I,U = LAD2.calculate_LAD_DTM(sp_pts,heights_rad,max_k,'spherical')
            LAD_rad_DTM[subplot_index,:,rr]=u.copy()
            u,n,I,U = LAD2.calculate_LAD_DTM_old(sp_pts,heights_rad,max_k,'spherical')
            LAD_rad_DTM_old[subplot_index,:,rr]=u.copy()
        lidar_return_profiles_adj[subplot_index,:,:n.shape[2]] = np.sum(n.copy(),axis=1)

        # now get MacArthur-Horn profiles
        heights,first_return_profile,n_ground_returns = LAD1.bin_returns(sp_pts, max_height, layer_thickness)
        LAD_MH[subplot_index,:] = LAD1.estimate_LAD_MacArthurHorn(first_return_profile, n_ground_returns, layer_thickness, kappa)
        # get the LiDAR penetration limits for the subplots
        penetration_lim[subplot_index,np.cumsum(first_return_profile)==0]=1.

    print("average point density = ", pt_count/10.**4, " pts/m^2")

    #------------------------------------------------------------------------------------
    # AVERAGE 1 ha PROFILES
    LAD_MH_mean = np.nansum(LAD_MH,axis=0)/(np.sum(np.isfinite(LAD_MH),axis=0)).astype('float')
    LAD_rad_mean = np.nansum(LAD_rad,axis=0)/(np.sum(np.isfinite(LAD_rad),axis=0)).astype('float')
    LAD_rad_DTM_mean = np.nansum(LAD_rad_DTM,axis=0)/(np.sum(np.isfinite(LAD_rad_DTM),axis=0)).astype('float')

    # CLEANING AND STORING
    # now we have looped through and created the different profiles, need to account for any NaN's and apply minimum height
    # to the LAD distributions
    # - remove all profile values below minimum height prior to comparison
    mask = heights <= minimum_height
    LAD_MH[:,mask]=np.nan
    mask = np.max(heights_rad)-heights_rad<=minimum_height
    LAD_rad[:,mask]=np.nan
    LAD_rad_DTM[:,mask]=np.nan
    LAD_rad_DTM_old[:,mask]=np.nan

    # now the profiles are ready to be stored into their relevant dictionaries
    lidar_profiles[Plot_name] = lidar_return_profiles.copy()
    lidar_profiles_adjusted[Plot_name] = lidar_return_profiles_adj.copy()
    penetration_limit[Plot_name] = penetration_lim.copy()

    MacArthurHorn_LAD[Plot_name] = LAD_MH.copy()
    radiative_LAD[Plot_name] = LAD_rad.copy()
    radiative_DTM_LAD[Plot_name] = LAD_rad_DTM.copy()
    radiative_DTM_LAD_old[Plot_name] = LAD_rad_DTM_old.copy()
    MacArthurHorn_LAD_mean[Plot_name] = LAD_MH_mean.copy()
    radiative_LAD_mean[Plot_name] = LAD_rad_mean.copy()
    radiative_DTM_LAD_mean[Plot_name] = LAD_rad_DTM_mean.copy()

    # Also ready to store their respective LAI profiles
    MacArthurHorn_LAI[Plot_name] = np.nansum(LAD_MH,axis=1)*layer_thickness
    radiative_LAI[Plot_name] = np.nansum(LAD_rad,axis=1)*layer_thickness
    radiative_DTM_LAI[Plot_name] = np.nansum(LAD_rad_DTM,axis=1)*layer_thickness
    MacArthurHorn_LAI_mean[Plot_name] = np.nansum(LAD_MH,axis=1)*layer_thickness
    radiative_LAI_mean[Plot_name] = np.nansum(LAD_rad,axis=1)*layer_thickness
    radiative_DTM_LAI_mean[Plot_name] = np.nansum(LAD_rad_DTM,axis=1)*layer_thickness
#----------------------------------------------------------------------------
np.savez('%splot_point_clouds.npz' % output_dir,plot_point_cloud)

np.savez('%slidar_canopy_profiles.npz' % output_dir,(MacArthurHorn_LAD,radiative_LAD,
        radiative_DTM_LAD,lidar_profiles,lidar_profiles_adjusted,penetration_limit))

np.savez('%slidar_PAI.npz' % output_dir,(MacArthurHorn_LAI,radiative_LAI,
        radiative_DTM_LAI,MacArthurHorn_LAI_mean,radiative_LAI_mean,
        radiative_DTM_LAI_mean))
