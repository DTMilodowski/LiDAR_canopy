###############################################################################################################
# This driver function analyses both LiDAR data and field inventory data to produce independent estimates of
# canopy structure.  These are compared against each other and their integrated LAD is compared against LAI
# estimates from hemispherical photographs.
###############################################################################################################
import numpy as np
import auxilliary_functions as aux
import inventory_based_LAD_profiles as field
import load_field_data as cen

#------------------------------------------------------------------------------------
# DIRECTORIES
# start by defining input files
allometry_file = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/Regional/Allometry/Crown_depth_data_SEAsia.csv'
field_file = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/Local/SAFE_DANUM_carbonplots_FieldMapcensus2016.csv'
stems_file = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/Local/SAFE_small_plot_tree_data.csv'

D1_stemcensus = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/Danum/carbon_plot_data_Danum_05042017_DTMformatting_DAN04.csv'
D2_stemcensus = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/Danum/carbon_plot_data_Danum_05042017_DTMformatting_DAN05.csv'
SAFE_stemcensus = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/SAFE_plots/SAFE_SmallTreeCensus_year1only.csv'

gps_pts_file = 'GPS_points_file_for_least_squares_fitting.csv'
subplot_coordinate_file = 'BALI_subplot_coordinates_corrected.csv'

# also define output directory (for saving figures)
output_dir = '/home/dmilodow/DataStore_DTM/BALI/PAPERS/PaperDrafts/EstimatingCanopyStructureBALI/Profiles/'

#------------------------------------------------------------------------------------
# PARAMETERS
# define important parameters for canopy profile estimation
Plots = [b'LF',b'E',b'Belian',b'Seraya',b'B North',b'B South',b'DC1',b'DC2']
#Plots = [b'DC1',b'DC2']
N_plots = len(Plots)
n_subplots=25
max_height = 80
layer_thickness = 1.
n_layers = np.ceil(max_height/layer_thickness)
minimum_height = 2.
plot_area = 10.**4
subplot_area = 20.*20.

heights = np.arange(0,max_height,layer_thickness)+layer_thickness

#------------------------------------------------------------------------------------
# DECLARATIONS
# define dictionaries to host the various canopy profiles and LAI estimates that will be produced

inventory_LAD = {}
inventory_LAD_std = {}
inventory_LAI = {}
inventory_LAD_all = {}

#------------------------------------------------------------------------------------
# LOADING DATA
subplot_polygons, subplot_labels = aux.load_boundaries(subplot_coordinate_file)
# load field data and retrieve allometric relationships
field_data = field.load_crown_survey_data(field_file)
DBH_BAAD, H_BAAD, D_BAAD = field.load_BAAD_crown_allometry_data(allometry_file)

# load stem data
#stem_data = field.load_SAFE_small_plot_data(stems_file)
DC1_stem_data = field.load_Danum_stem_census(D1_stemcensus) # all subplots
DC2_stem_data = field.load_Danum_stem_census(D2_stemcensus) # all subplots
SAFE_stem_data = field.load_SAFE_small_stem_census(SAFE_stemcensus) # subset of subplots

#------------------------------------------------------------------------------------
# MAIN ANALYSIS
# FIELD INVENTORY ANALYSIS - MONTE-CARLO UNCERTAINTY PROPAGATION
# First get allometric relationships required
# 1) height-depth allometry
a, b, CF, r_sq, p, temp, PI_u, PI_l = field.log_log_linear_regression(H_BAAD,D_BAAD)
# 2) DBH-Crown area allometry
a_A, b_A, CF_A, r_sq_A, p_A, temp, PI_u_A, PI_l_A = field.log_log_linear_regression(field_data['DBH_field'],field_data['CrownArea'])
# 3) DBH-Height allometry
a_ht, b_ht, CF_ht, r_sq_ht, p_ht, temp, PI_u_ht, PI_l_ht = field.log_log_linear_regression(field_data['DBH_field'],field_data['Height_field'])

for pp in range(0,N_plots):
    for ss in range(0,n_subplots):# mask out dead and broken trees
        dead_mask = np.all((field_data['dead_flag1']==-1,field_data['dead_flag2']==-1,field_data['dead_flag3']==-1),axis=0)
        brokenlive_mask = field_data['brokenlive_flag']==-1
        mask = np.all((field_data['subplot']==ss-1,field_data['plot']==Plots[pp],np.isfinite(field_data['DBH_field']),np.isfinite(field_data['Xfield']),dead_mask,brokenlive_mask),axis=0)


BAAD_data={}
BAAD_data['D']=D_BAAD
BAAD_data['Ht']=H_BAAD
BAAD_data['DBH']=DBH_BAAD
error = {}
error['Ht']=[.1,.1]
error['DBH']=[0.,.02]
error['Area']=[0.,.05]

n_iter = 100
for pp in range(0,N_plots):
    print(Plots[pp])
    Plot_name=Plots[pp]
    # set up array to host inventory profiles

    # mask out dead and broken trees
    dead_mask = np.all((field_data['dead_flag1']==-1,field_data['dead_flag2']==-1,
                                        field_data['dead_flag3']==-1),axis=0)
    brokenlive_mask = (field_data['brokenlive_flag']==-1)
    mask = np.all((field_data['plot']==Plot_name,np.isfinite(field_data['DBH_field']),
                                            dead_mask,brokenlive_mask),axis=0)

    Ht = field_data['Height_field'][mask]
    DBH = field_data['DBH_field'][mask]
    Area = field_data['CrownArea'][mask]
    x0 = field_data['Xfield'][mask]
    y0 = field_data['Yfield'][mask]

    # now building the canopy model
    buff = 20
    xmin = np.nanmin(field_data['Xfield'][mask])
    xmax = np.nanmax(field_data['Xfield'][mask])
    ymin = np.nanmin(field_data['Yfield'][mask])
    ymax = np.nanmax(field_data['Yfield'][mask])
    x = np.arange(xmin-buff,xmax+buff,1.)+0.5
    y = np.arange(ymin-buff,ymax+buff,1.)+0.5
    z = np.arange(0,80.,1.)+0.5

    # ITERATE MONTE-CARLO PROCEDURE
    #------------------------------------------------------------------------------------
    field_profiles = field.calculate_crown_volume_profiles_mc(x,y,z,x0,y0,Ht,DBH,Area,
                                        a_ht,b_ht,a_A,b_A,a,b,
                                        field_data,BAAD_data,n_iter=n_iter)
    #field_profiles = calculate_crown_volume_profiles_mc_with_measurement_error(x,y,z,
    #                                    x0,y0,Ht,DBH,Area,a_ht,b_ht,a_A,b_A,a_D,b_D,error,
    #                                    field_data,BAAD_data,n_iter=n_iter)

    # add small stem contributions
    smallstem_profiles = np.zeros((n_subplots,heights.size))
    for ss in range(0,n_subplots):
        subplot_index = int(subplot_labels[Plot_name][ss]-1)
        if Plot_name == b'DC1':
            Ht,Area,Depth,StemDensity = field.calculate_crown_dimensions_for_stem_distributions(DC1_stem_data['dbh'],
                                                DC1_stem_data['stem_density'][:,subplot_index],
                                                a_ht, b_ht, CF_ht, a_A, b_A, CF_A, a, b, CF)
        elif Plot_name == b'DC2':
            Ht,Area,Depth,StemDensity = field.calculate_crown_dimensions_for_stem_distributions(DC2_stem_data['dbh'],
                                                DC2_stem_data['stem_density'][:,subplot_index],
                                                a_ht, b_ht, CF_ht, a_A, b_A, CF_A, a, b, CF)
        else:
            Ht,Area,Depth,StemDensity = field.calculate_crown_dimensions_for_stem_distributions(SAFE_stem_data[Plot_name]['dbh'],
                                                SAFE_stem_data[Plot_name]['stem_density'][:,subplot_index],
                                                a_ht, b_ht, CF_ht, a_A, b_A, CF_A, a, b, CF)

        smallstem_profiles[subplot_index,:] = field.calculate_LAD_profiles_ellipsoid_from_stem_size_distributions(heights,
                                            Area, Depth, Ht, StemDensity)

    field_profiles[:]+=np.mean(smallstem_profiles,axis=0)

    inventory_LAD[Plot_name] = np.mean(field_profiles,axis=0)
    inventory_LAD_std[Plot_name] = np.std(field_profiles,axis=0)
    inventory_LAI[Plot_name] = np.sum(field_profiles,axis=1)*layer_thickness
    inventory_LAD_all[Plot_name] = field_profiles.copy()

np.savez('%sinventory_canopy_profiles.npz' % output_dir,(inventory_LAD,
                            inventory_LAD_std,inventory_LAI,inventory_LAD_all))
