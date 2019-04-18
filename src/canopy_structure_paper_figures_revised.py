###############################################################################################################
# This driver function analyses both LiDAR data and field inventory data to produce independent estimates of
# canopy structure.  These are compared against each other and their integrated LAD is compared against LAI
# estimates from hemispherical photographs.
###############################################################################################################
import numpy as np
import sys
import auxilliary_functions as aux
import load_field_data as cen
import canopy_structure_plots as csp
import inventory_based_LAD_profiles as field

#------------------------------------------------------------------------------------
# DIRECTORIES
# start by defining input files
subplot_coordinate_file = 'BALI_subplot_coordinates_corrected.csv'
allometry_file = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/Regional/Allometry/Crown_depth_data_SEAsia.csv'
field_file = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/Local/SAFE_DANUM_carbonplots_FieldMapcensus2016.csv'
gps_pts_file = 'GPS_points_file_for_least_squares_fitting.csv'

# also define output directory (for saving figures)
data_dir = '/home/dmilodow/DataStore_DTM/BALI/PAPERS/PaperDrafts/EstimatingCanopyStructureBALI/Profiles/'
output_dir = '/home/dmilodow/DataStore_DTM/BALI/PAPERS/PaperDrafts/EstimatingCanopyStructureBALI/FiguresRevised/'

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
temp = np.load('%slidar_canopy_profiles.npz' % data_dir)['arr_0'][()]
MacArthurHorn_PAD=temp[0]
radiative_PAD=temp[1]
radiative_DTM_PAD=temp[2]
lidar_profiles=temp[3]
lidar_profiles_adjusted=temp[4]
penetration_limit=temp[5]
temp=None

MacArthurHorn_PAD_mean = {}
radiative_PAD_mean= {}
radiative_DTM_PAD_mean= {}
for pp in range(0,N_plots):
    MacArthurHorn_PAD_mean[Plots[pp]] = np.nansum(MacArthurHorn_PAD[Plots[pp]],axis=0)/(np.sum(np.isfinite(MacArthurHorn_PAD[Plots[pp]]),axis=0)).astype('float')
    radiative_PAD_mean[Plots[pp]] = np.nansum(radiative_PAD[Plots[pp]],axis=0)/(np.sum(np.isfinite(radiative_PAD[Plots[pp]]),axis=0)).astype('float')
    radiative_DTM_PAD_mean[Plots[pp]] = np.nansum(radiative_DTM_PAD[Plots[pp]],axis=0)/(np.sum(np.isfinite(radiative_DTM_PAD[Plots[pp]]),axis=0)).astype('float')


# Load LiDAR PAI
temp = np.load('%slidar_PAI.npz' % data_dir)['arr_0'][()]
MacArthurHorn_PAI=temp[0]
radiative_PAI=temp[1]
radiative_DTM_PAI=temp[2]
MacArthurHorn_PAI_mean=temp[3]
radiative_PAI_mean=temp[4]
radiative_DTM_PAI_mean=temp[5]
temp=None

# Load Inventory profiles
temp = np.load('%sinventory_canopy_profiles.npz' % data_dir)['arr_0'][()]
inventory_PAD=temp[0]
inventory_PAD_std=temp[1]
inventory_PAI=temp[2]
#inventory_PAI_std=temp['arr_3'][()]
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

"""
# Figure 3 - example crown model
"""
figure_number = 3
figure_name = output_dir+'fig3_crown_model_example'
Plot_name = b'Belian'
angle = 45.
csp.plot_canopy_model(figure_number,figure_name,Plot_name,field_data,angle,
                    a_ht, b_ht, CF_ht, a_A, b_A, CF_A, a, b, CF)

"""
# Figure 4 - Allometric models; include confidence intervals, and add vertical band
# illustrating the 10 cm DBH cutoff
"""
figure_name = output_dir + 'Fig4_allometric_relationships.png'
figure_number = 4
csp.plot_allometric_relationships(figure_name,figure_number,field_file,allometry_file)

#-------------------------------
# RESULTS - STRUCTURAL CHANGES
#           ACROSS GRADIENT
#-------------------------------
"""
# Figure 5 - Point clouds and profiles across degradation gradient
"""
figure_name = output_dir + 'Fig5_pointclouds_and_profiles.png'
figure_number = 5
gps_pts_file = 'GPS_points_file_for_least_squares_fitting.csv'
csp.plot_point_clouds_and_profiles(figure_name,figure_number, gps_pts_file,
                        plot_point_cloud,heights,heights_rad, lidar_profiles,
                        MacArthurHorn_PAD,MacArthurHorn_PAD_mean,radiative_DTM_PAD,
                        radiative_DTM_PAD_mean,inventory_PAD)

"""
# Figure 6 - Cross-plot canopy layers
"""
figure_name = output_dir + 'Fig6_crossplot_LiDAR_PAD_residual_profiles.png'
figure_number = 6
csp.plot_canopy_layer_residuals(figure_name,figure_number,heights,MacArthurHorn_PAD,
                MacArthurHorn_PAD_mean,radiative_DTM_PAD,radiative_DTM_PAD_mean)

"""
# Figure 7 - PAI plotted against basal area
"""
figure_name = output_dir + 'Fig7_PAI_vs_basal_area.png'
figure_number = 7

census_file = '/home/dmilodow/DataStore_DTM/BALI/BALI_Cplot_data/SAFE_CarbonPlots_TreeCensus.csv'
census = cen.collate_plot_level_census_data(census_file)

BA = {}
Plots_SAFE = ['Belian', 'Seraya', 'LF', 'E','B North', 'B South']
Plots_Danum = ['DC1', 'DC2']
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

for pp in range(0,len(Plots_SAFE)):
    temp_BA = np.zeros(n_subplots)
    for ss in range(0,n_subplots):
        # check basal area
        temp_BA[ss] = census[Plots_SAFE[pp]]['BasalArea'][ss,0]*n_subplots/100**2
    temp_BA[np.isnan(temp_BA)]=0
    BA[Plots_SAFE[pp]]=temp_BA.copy()

for pp in range(0,len(Plots_Danum)):
    temp_BA = np.zeros(n_subplots)
    for ss in range(0,n_subplots):
        indices = np.all((field_data['plot']==str.encode(Plots_Danum[pp]),
                                field_data['subplot']==ss+1,
                                np.isfinite(field_data['DBH'])),axis=0)
        temp_BA[ss] = np.sum((field_data['DBH'][indices]/2)**2*np.pi)*n_subplots/100.**2
    BA[Plots_Danum[pp]]=temp_BA.copy()

csp.plot_LAI_vs_basal_area(figure_name,figure_number,MacArthurHorn_PAD,MacArthurHorn_PAD_mean,
                            radiative_DTM_PAD,radiative_DTM_PAD_mean,BA,plot_marker,plot_label,
                            plot_colour)

#-------------------------------
# RESULTS - SENSITIVITY ANALYSIS
# see sensitivity analysis plots
#-------------------------------
# Figure 8 - Sensitivity analysis of vertical profiles to spatial resolution
# Comparison of OG vs Moderately Logged vs. Heavily Logged
# <see sensitivity analysis script>

# Figure 9 - Sensitivity analysis of unsampled voxels
# <see sensitivity analysis script>

# Figure 10 - Sensitivity analysis of vertical profiles to point density
# <see sensitivity analysis script>

# Figure 11 - Summary plots from sensitivity analysis for PAI vs. resolution
# and point density
# <see sensitivity analysis script>

#-------------------------------
# SUPPLEMENT
# METHODS
#-------------------------------
"""
# Figure S1 - "transmission ratio"
"""
figure_number = 111
figure_name = output_dir+'figS1_transmittance_ratios.png'
csp.plot_transmittance_ratio(figure_number,figure_name,all_lidar_pts)

"""
# Figure S2 - comparison of Detto vs. modified algorithm
"""
figure_number = 112
figure_name = output_dir+'figS2_LiDAR_profiles_comparison.png'
csp.plot_LiDAR_profiles_comparison(figure_name,figure_number,heights,heights_rad,
                        lidar_profiles,MacArthurHorn_PAD,MacArthurHorn_PAD_mean,
                        radiative_PAD,radiative_PAD_mean,
                        radiative_DTM_PAD,radiative_DTM_PAD_mean)

#-------------------------------
# SUPPLEMENT
# RESULTS
#-------------------------------
"""
# Figure S3 comparison of profiles for the two Danum sites
"""
figure_number = 113
figure_name = output_dir+'Fig3_pointclouds_and_profiles_Danum.png'
csp.plot_point_clouds_and_profiles_Danum(figure_name,figure_number, gps_pts_file,
                        plot_point_cloud,heights,heights_rad, lidar_profiles,
                        MacArthurHorn_PAD,MacArthurHorn_PAD_mean,radiative_DTM_PAD,
                        radiative_DTM_PAD_mean,inventory_PAD)

# Figure S4 - sensitivity analysis, confidence interval sensitivity to resolution

# Figure S5 - sensitivity analysis, confidence interval sensitivity to density
