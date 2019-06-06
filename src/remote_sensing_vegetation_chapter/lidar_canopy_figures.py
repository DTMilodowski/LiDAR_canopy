###############################################################################################################
# This driver function analyses both LiDAR data and field inventory data to produce independent estimates of
# canopy structure.  These are compared against each other and their integrated LAD is compared against LAI
# estimates from hemispherical photographs.
###############################################################################################################
import numpy as np
from scipy import stats
import xarray as xr
#import auxilliary_functions as aux
#import canopy_structure_plots as csp
from matplotlib import pyplot as plt
import seaborn as sns

#------------------------------------------------------------------------------------
# DIRECTORIES
# start by defining input files
las_file = '/home/dmilodow/DataStore_DTM/OTHER/RemoteSensingofVegetation_chapter/LiDAR/LNF_FR.las'
dem_file = '/home/dmilodow/DataStore_DTM/OTHER/RemoteSensingofVegetation_chapter/LiDAR/LNF_FR_dem.tif'
dsm_file = '/home/dmilodow/DataStore_DTM/OTHER/RemoteSensingofVegetation_chapter/LiDAR/LNF_FR_dsm.tif'
output_dir = '/home/dmilodow/DataStore_DTM/OTHER/RemoteSensingofVegetation_chapter/figures/'

#-----------------------------------------------------------------------------------
# Plot rasters (Figure 2 from LiDAR work)
dsm = xr.open_rasterio(dsm_file)
dsm.values[dsm.values>905]=np.nan
dsm.values[dsm.values < 0] = np.nan
x = dsm.coords['x'].values.copy()
y = dsm.coords['y'].values.copy()
origin = [x[0],y[0]]
dsm.coords['x']=x-x[0]
dsm.coords['y']=y-y[-1]

dem = xr.open_rasterio(dem_file)
dem.values[dem.values < 0] = np.nan
dem.coords['x']=x-x[0]
dem.coords['y']=y-y[-1]


chm = dsm-dem
mask = np.any((chm.values<0,chm.values>100),axis=0)
chm.values[mask] = np.nan

elev_range =[np.nanmin(dem.values),np.nanmax(dsm.values)]
height_range =[np.nanmin(chm.values),60]
fig,axes= plt.subplots(1,3,figsize = [9,5])
dsm.plot(vmin=elev_range[0],vmax=elev_range[1],cmap = 'viridis',ax=axes[0],
        add_colorbar=True,cbar_kwargs={'label': 'Elevation A.S.L. / m',
                    'orientation':'horizontal'})
dem.plot(vmin=elev_range[0],vmax=elev_range[1],cmap = 'viridis',ax=axes[1],
        add_colorbar=True,cbar_kwargs={'label': 'Elevation A.S.L. / m',
                    'orientation':'horizontal'})
chm.plot(vmin=height_range[0],vmax=height_range[1],cmap = 'viridis',ax=axes[2],
        add_colorbar=True,cbar_kwargs={'label': 'Height above ground / m',
                    'orientation':'horizontal'})

axes[0].set_title('Digital Surface Model',fontsize=12)
axes[1].set_title('Digital Elevation Model',fontsize=12)
axes[2].set_title('Canopy Height Model',fontsize=12)

yticklabels = axes[1].get_yticklabels() + axes[2].get_yticklabels()
axes[0].set_ylabel('Distance / m')
axes[1].set_ylabel('')
axes[2].set_ylabel('')
axes[1].set_xlabel('Distance / m')
axes[0].set_xlabel('')
axes[2].set_xlabel('')
plt.setp(yticklabels,visible=False)
for ax in axes:
    ax.yaxis.set_ticks_position('both')
    ax.set_aspect("equal")

fig.show()

#------------------------------------------------------------------------------------
# PARAMETERS
# define important parameters for canopy profile estimation
max_height = 80
layer_thickness = 1.
heights = np.arange(0,max_height,layer_thickness)+layer_thickness
n_layers = heights.size

"""plot_width = 100.
sample_res = np.array([2.,5.,10.,20.,50.])
keys = ['2m','5m','10m','20m','50m']
"""
kappa = 1.

# 2 ha sample to read in

N=6206000.
S=6205900.
W=344520.
E=344620.

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
inventory_PAD_all=temp[3]
temp = None

#===============================================================================
# Summary statistics
table_plots = [b'Belian',b'Seraya',b'DC1',b'DC2',b'E',b'LF',b'B North',b'B South']
print("Plot    \tMH\t+/-\trad_2\t+/-\trad_3\t+/-\tcv\t+/-")
for pp,plot in enumerate(table_plots):
    mh = np.mean(MacArthurHorn_PAI[plot])
    mh_s = stats.sem(MacArthurHorn_PAI[plot])
    r2 = np.mean(radiative_DTM_PAI[plot][:,1])
    r2_s = stats.sem(radiative_DTM_PAI[plot][:,1])
    r3 = np.mean(radiative_DTM_PAI[plot][:,2])
    r3_s = stats.sem(radiative_DTM_PAI[plot][:,2])
    cv = np.mean(inventory_PAI[plot])
    cv_s = stats.sem(inventory_PAI[plot])

    print('%s       \t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.2f\t' % (plot,mh,mh_s,
                            r2,r2_s,r3,r3_s,cv,cv_s))
    #print('%s\t    %.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f' % (plot,mh,mh_s,
    #                        r2,r2_s,r3,r3_s,cv))

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

"""
# Figure 8 - PAD distributions for distinct canopy layers
"""
figure_name = output_dir + 'Fig8_canopy_sublevel_PAD.png'
figure_number = 8
csp.plot_PAD_distributions_for_canopy_subdivisions(figure_name,figure_number,
                            heights,heights_rad,MacArthurHorn_PAD,
                            radiative_PAD,radiative_DTM_PAD)

"""
# Figure 11 - Cumulative PAD with Depth
"""
figure_name = output_dir + 'Fig11_cumulative_PAD_with_depth.png'
figure_number = 11
csp.plot_cumulative_PAD_vs_depth(figure_name,figure_number,MacArthurHorn_PAD,
                        radiative_DTM_PAD, method=0)
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
# Figure S2 - "transmission ratio"
"""
figure_number = 112
figure_name = output_dir+'figS2_transmittance_ratios.png'
csp.plot_transmittance_ratio(figure_number,figure_name,all_lidar_pts)

"""
# Figure S2 - comparison of Detto vs. modified algorithm
"""
figure_number = 113
figure_name = output_dir+'figS3_LiDAR_profiles_comparison.png'
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
