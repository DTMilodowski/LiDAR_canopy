###############################################################################################################
# This driver function analyses both LiDAR data and field inventory data to produce independent estimates of
# canopy structure.  These are compared against each other and their integrated LAD is compared against LAI
# estimates from hemispherical photographs.
###############################################################################################################
import numpy as np
from scipy import stats
import sys
import auxilliary_functions as aux
import canopy_structure_plots as csp
import inventory_based_LAD_profiles as field
import LiDAR_io as io

#------------------------------------------------------------------------------------
# DIRECTORIES
# start by defining input files
subplot_coordinate_file = 'BALI_subplot_coordinates_corrected.csv'
allometry_file = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/Regional/Allometry/Crown_depth_data_SEAsia.csv'
field_file = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/Local/SAFE_DANUM_carbonplots_FieldMapcensus2016.csv'
gps_pts_file = 'GPS_points_file_for_least_squares_fitting.csv'

# also define output directory (for saving figures)
data_dir = '/home/dmilodow/DataStore_DTM/BALI/PAPERS/PaperDrafts/EstimatingCanopyStructureBALI/Profiles/'
output_dir = '/home/dmilodow/DataStore_DTM/BALI/PAPERS/PaperDrafts/EstimatingCanopyStructureBALI/FiguresRevised/201910/'

#------------------------------------------------------------------------------------
# PARAMETERS
# define important parameters for canopy profile estimation
Plots = [b'LF',b'E',b'Belian',b'Seraya',b'B North',b'B South',b'DC1',b'DC2']
#Plots = ['B North']
N_plots = len(Plots)
n_subplots = 25
leaf_angle_dist = 'spherical'
max_height = 80
max_return = 3
layer_thickness = 1.
n_layers = np.ceil(max_height/layer_thickness)
minimum_height = 2.
plot_area = 10.**4
subplot_area = 20.*20.
kappa = 0.70

heights = np.arange(0,max_height,layer_thickness)+layer_thickness
heights_rad = np.arange(0,max_height+layer_thickness,layer_thickness)

#------------------------------------------------------------------------------------
# LOADING DATA
# load coordinates and lidar points for target areas
subplot_polygons, subplot_labels = aux.load_boundaries(subplot_coordinate_file)

# load field data and retrieve allometric relationships
field_data = field.load_crown_survey_data(field_file)

# Load LiDAR point clouds for the plots
plot_point_cloud= np.load('%splot_point_clouds.npz' % data_dir)['arr_0'][()]

# Load LiDAR canopy profiles
temp = np.load('%slidar_canopy_profiles_adaptive.npz' % data_dir)['arr_0'][()]
MacArthurHorn_PAD=temp[0]
MacArthurHorn_wt_PAD=temp[1]
radiative_PAD=temp[2]
radiative_DTM_PAD=temp[3]
lidar_profiles=temp[4]
lidar_profiles_adjusted=temp[5]
penetration_limit=temp[6]
temp=None

MacArthurHorn_PAD_mean = {}
MacArthurHorn_wt_PAD_mean = {}
radiative_PAD_mean= {}
radiative_DTM_PAD_mean= {}
for pp in range(0,N_plots):
    MacArthurHorn_PAD_mean[Plots[pp]] = np.nansum(MacArthurHorn_PAD[Plots[pp]],axis=0)/(np.sum(np.isfinite(MacArthurHorn_PAD[Plots[pp]]),axis=0)).astype('float')
    MacArthurHorn_wt_PAD_mean[Plots[pp]] = np.nansum(MacArthurHorn_wt_PAD[Plots[pp]],axis=0)/(np.sum(np.isfinite(MacArthurHorn_PAD[Plots[pp]]),axis=0)).astype('float')
    radiative_PAD_mean[Plots[pp]] = np.nansum(radiative_PAD[Plots[pp]],axis=0)/(np.sum(np.isfinite(radiative_PAD[Plots[pp]]),axis=0)).astype('float')
    radiative_DTM_PAD_mean[Plots[pp]] = np.nansum(radiative_DTM_PAD[Plots[pp]],axis=0)/(np.sum(np.isfinite(radiative_DTM_PAD[Plots[pp]]),axis=0)).astype('float')

# Load LiDAR PAI
temp = np.load('%slidar_PAI_adaptive.npz' % data_dir)['arr_0'][()]
MacArthurHorn_PAI=temp[0]
MacArthurHorn_wt_PAI=temp[1]
radiative_PAI=temp[2]
radiative_DTM_PAI=temp[3]
temp=None

# Load Inventory profiles
temp = np.load('%sinventory_canopy_profiles.npz' % data_dir)['arr_0'][()]
inventory_PAD=temp[0]
inventory_PAD_std=temp[1]
inventory_PAI=temp[2]
inventory_PAD_all=temp[3]
temp = None

#===============================================================================
# NOW MAKE PLOTS

#-------------------------------
# INTRODUCTION & METHODS
#-------------------------------
"""
# Figure 1 - Location map, with Hansen data and plot locations
"""
figure_name = output_dir+'Fig1_Location_map.png'
figure_number = 1
csp.plot_location_map(figure_name,figure_number)

"""
# Figure 2 sample point cloud - coloured by return number
"""
figure_name = output_dir+'Fig2_sample_point_cloud.png'
figure_number = 2
csp.plot_point_cloud(figure_name,figure_number,gps_pts_file,plot_point_cloud)


#-------------------------------
# RESULTS - STRUCTURAL CHANGES
#           ACROSS GRADIENT
#-------------------------------
"""
# Figure 3 - PAI plotted against basal area
"""
figure_name = output_dir + 'Fig3_PAI_vs_basal_area.png'
figure_number = 3

# Basal area (m^2 / ha) and standard errors
# data from table 1 of Riutta et al, GCB, 2018
BA = {}
BA['Belian']=(41.6,3.59)
BA['Seraya']=(34.7,2.74)
BA['LF']= (19.3,1.7)
BA['E']= (19.6,1.88)
BA['B North']=(11.1,1.81)
BA['B South']=(6.81,1.00)
BA['DC1']=(32.0,3.3)
BA['DC2']=(30.6,3.37)

colour = ['#46E900','#1A2BCE','#E0007F']
plot_colour = {}
plot_colour['Belian']=colour[0];plot_colour['Seraya']=colour[0];plot_colour['DC1']=colour[0]
plot_colour['DC2']=colour[0];plot_colour['LF']=colour[1];plot_colour['E']=colour[1]
plot_colour['B North']=colour[2];plot_colour['B South']=colour[2]

plot_marker = {}
plot_marker['Belian']='o';plot_marker['Seraya']='v';plot_marker['DC1']='^';plot_marker['DC2']='s'
plot_marker['LF']='o';plot_marker['E']='v';plot_marker['B North']='o';plot_marker['B South']='v'
plot_label = {}
plot_label['Belian']='MLA01';plot_label['Seraya']='MLA02';plot_label['DC1']='DAN04';plot_label['DC2']='DAN05'
plot_label['LF']='SAF04';plot_label['E']='SAF03';plot_label['B North']='SAF02';plot_label['B South']='SAF01'
csp.plot_LAI_vs_basal_area(figure_name,figure_number,MacArthurHorn_PAD,MacArthurHorn_PAD_mean,
                            MacArthurHorn_wt_PAD,MacArthurHorn_wt_PAD_mean,
                            radiative_DTM_PAD,radiative_DTM_PAD_mean,BA,plot_marker,plot_label,
                            plot_colour)

"""
# Figure 4 - Point clouds and profiles across degradation gradient
"""
figure_name = output_dir + 'Fig4_pointclouds_and_profiles.png'
figure_number = 4
gps_pts_file = 'GPS_points_file_for_least_squares_fitting.csv'
csp.plot_point_clouds_and_profiles(figure_name,figure_number, gps_pts_file,
                        plot_point_cloud,heights,heights_rad, lidar_profiles,
                        MacArthurHorn_PAD,MacArthurHorn_PAD_mean,
                        MacArthurHorn_wt_PAD,MacArthurHorn_wt_PAD_mean,
                        radiative_DTM_PAD,radiative_DTM_PAD_mean,
                        inventory_PAD,inventory_PAD_all)

"""
# Figure 5 - Sensitivity of PAI estimates to pulse density and resolution
# Figure 6 - Sensitivity analysis of unsampled voxels
# see sensitivity_analysis_figures_revised
"""
figure_name = output_dir + 'Fig5_PAI_sensitivity.png'

"""
# Figure 7 - Niche availability
# Plotting Shannon Index (panel a) and overstory PAD (panel b)
"""
figure_name = output_dir + 'Fig7_niche_availability.png'
figure_number = 7
csp.plot_niche_availability(figure_name,figure_number,MacArthurHorn_PAD,heights)

"""
#-------------------------------
# SUPPLEMENT
#-------------------------------
# Figure S1 - example crown model
"""
field_data = field.load_crown_survey_data(field_file)
a, b, CF, r_sq, p, H, D, H_i, PI_u, PI_l = field.retrieve_crown_allometry(allometry_file)
a_ht, b_ht, CF_ht, a_A, b_A, CF_A = field.calculate_allometric_equations_from_survey(field_data)
figure_number = 111
figure_name = output_dir+'figS1_crown_model_example'
Plot_name = b'Belian'
angle = 45.
csp.plot_canopy_model(figure_number,figure_name,Plot_name,field_data,angle,
a_ht, b_ht, CF_ht, a_A, b_A, CF_A, a, b, CF)


"""
# Figure S2 - Allometric models; include confidence intervals, and add vertical band
# illustrating the 10 cm DBH cutoff
"""
figure_name = output_dir + 'FigS2_allometric_relationships.png'
figure_number = 112
csp.plot_allometric_relationships(figure_name,figure_number,field_file,allometry_file)

"""
# Figure S3 - Leaf angle distributions
# see canopy_structure_paper_leaf_angle_dist_revised
"""

"""
# Figure S4 - "transmission ratio"
"""
las_file = 'Carbon_plot_point_cloud_buffer.las'
all_lidar_pts = io.load_lidar_data(las_file)
figure_number = 114
figure_name = output_dir+'figS4_transmittance_ratios.png'
csp.plot_transmittance_ratio(figure_number,figure_name,all_lidar_pts)


"""
# Figure S5- comparison of Detto vs. modified algorithm
"""
figure_number = 115
figure_name = output_dir+'figS5_LiDAR_profiles_comparison_test.png'
csp.plot_LiDAR_profiles_comparison(figure_name,figure_number,heights,heights_rad,
lidar_profiles,MacArthurHorn_PAD,MacArthurHorn_PAD_mean,
MacArthurHorn_wt_PAD,MacArthurHorn_wt_PAD_mean,
radiative_PAD,radiative_PAD_mean,
radiative_DTM_PAD,radiative_DTM_PAD_mean)

"""
# Figure S6 comparison of profiles for the two Danum sites
"""
figure_number = 116
figure_name = output_dir+'FigS6_pointclouds_and_profiles_Danum.png'
csp.plot_point_clouds_and_profiles_Danum(figure_name,figure_number, gps_pts_file,
                        plot_point_cloud,heights,heights_rad, lidar_profiles,
                        MacArthurHorn_PAD,MacArthurHorn_PAD_mean,
                        MacArthurHorn_wt_PAD,MacArthurHorn_wt_PAD_mean,
                        radiative_DTM_PAD,radiative_DTM_PAD_mean,
                        inventory_PAD,inventory_PAD_all)

"""
# Figure S7 - sensitivity of 1 ha profiles to pulse density
# Figure S8 - relative CIs across profiles for different pulse densities
# Figure S9 - sensitivity of 1 ha profiles to grid resolution
# Figure S10- relative CIs across profiles for different grid resolutions
# see sensitivity_analysis_figures_revised
"""
"""
#===============================================================================
# Summary statistics
"""
table_plots = [b'Belian',b'Seraya',b'DC1',b'DC2',b'E',b'LF',b'B North',b'B South']
print("Plot    \tMH\t+/-\tMHwt\t+/-\trad_2\t+/-\trad_3\t+/-\tcv\t+/-")
for pp,plot in enumerate(table_plots):
    mh = np.mean(MacArthurHorn_PAI[plot])
    mh_s = stats.sem(MacArthurHorn_PAI[plot])
    mh_wt = np.mean(MacArthurHorn_wt_PAI[plot])
    mh_wt_s = stats.sem(MacArthurHorn_wt_PAI[plot])
    r2 = np.mean(radiative_DTM_PAI[plot][:,1])
    r2_s = stats.sem(radiative_DTM_PAI[plot][:,1])
    r3 = np.mean(radiative_DTM_PAI[plot][:,2])
    r3_s = stats.sem(radiative_DTM_PAI[plot][:,2])
    cv = np.mean(inventory_PAI[plot])
    cv_s = stats.sem(inventory_PAI[plot])

    print('%s     \t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.2f\t' % (plot,mh,mh_s,
                            mh_wt,mh_wt_s,r2,r2_s,r3,r3_s,cv,cv_s))
