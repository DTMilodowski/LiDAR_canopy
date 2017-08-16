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
import inventory_based_LAD_profiles as field
import least_squares_fitting as lstsq
from scipy import stats

from matplotlib import rcParams
from datetime import datetime

sys.path.append('/home/dmilodow/DataStore_DTM/BALI/SPA_BALI_data_and_analysis/scripts/field_data/')
import load_field_data as cen

# Set up some basiic parameters for the plots
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = 8
rcParams['legend.numpoints'] = 1
axis_size = rcParams['font.size']+2

# start by defining input files
las_file = 'Carbon_plot_point_cloud_buffer.las'
subplot_coordinate_file = 'BALI_subplot_coordinates_corrected.csv'
allometry_file = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/Regional/Allometry/Crown_depth_data_SEAsia.csv'
field_file = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/Local/SAFE_DANUM_carbonplots_FieldMapcensus2016.csv'
LAI_file = '/home/dmilodow/DataStore_DTM/BALI/BALI_Cplot_data/SAFE_CarbonPlots_LAI_fromHemisphericalPhotos.csv'
stems_file = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/Local/SAFE_small_plot_tree_data.csv'

D1_stemcensus = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/Danum/carbon_plot_data_Danum_05042017_DTMformatting_DAN04.csv'
D2_stemcensus = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/Danum/carbon_plot_data_Danum_05042017_DTMformatting_DAN05.csv'
SAFE_stemcensus = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/SAFE_plots/SAFE_SmallTreeCensus_year1only.csv'

# also define output directory (for saving figures)
output_dir = '/home/dmilodow/DataStore_DTM/BALI/PAPERS/PaperDrafts/EstimatingCanopyStructureBALI/Figures/'

# define important parameters for canopy profile estimation
Plots = ['LF','E','Belian','Seraya','B North','B South','DC1','DC2']
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
beta = 0.6 # morphological parameter in inventory-based canopy model
kappa = 0.76

# define dictionaries to host the various canopy profiles and LAI estimates that will be produced
MacArthurHorn_LAD = {}
radiative_LAD = {}
radiative_DTM_LAD = {}
inventory_LAD = {}
lidar_profiles ={}
lidar_profiles_adjusted ={}

MacArthurHorn_LAI = {}
radiative_LAI = {}
radiative_DTM_LAI = {}
inventory_LAI = {}
Hemisfer_LAI = {}

plot_point_cloud = {}

# load coordinates and lidar points for target areas
subplot_polygons, subplot_labels = aux.load_boundaries(subplot_coordinate_file)
all_lidar_pts = lidar.load_lidar_data(las_file)

# load field data and retrieve allometric relationships
field_data = field.load_crown_survey_data(field_file)
a, b, CF, r_sq, p, H, D = field.retrieve_crown_allometry(allometry_file)
a_ht, b_ht, CF_ht, a_A, b_A, CF_A = field.calculate_allometric_equations_from_survey(field_data)

# load stem data
#stem_data = field.load_SAFE_small_plot_data(stems_file)
DC1_stem_data = field.load_Danum_stem_census(D1_stemcensus) # all subplots
DC2_stem_data = field.load_Danum_stem_census(D2_stemcensus) # all subplots
SAFE_stem_data = field.load_SAFE_small_stem_census(SAFE_stemcensus) # subset of subplots

# load LAI estimates from hemiphotos
field_LAI = aux.load_field_LAI(LAI_file)

# loop through all plots to be analysed
for pp in range(0,N_plots):
    print Plots[pp]
    Plot_name=Plots[pp]
    # clip LiDAR point cloud to plot level (this makes subsequent processing much faster)
    n_coord_pairs = subplot_polygons[Plot_name].shape[0]*subplot_polygons[Plot_name].shape[1]
    coord_pairs = subplot_polygons[Plot_name].reshape(n_coord_pairs,2)
    bbox_polygon = aux.get_bounding_box(coord_pairs)
    plot_lidar_pts = lidar.filter_lidar_data_by_polygon(all_lidar_pts,bbox_polygon)    
    plot_point_cloud[Plots[pp]]=plot_lidar_pts.copy()
    print "canopy height = ", np.percentile(plot_lidar_pts[plot_lidar_pts[:,3]==1,2],99), "m"

    # get some subplot-level information
    n_subplots = subplot_polygons[Plot_name].shape[0]
    # set up some arrays to host the radiative transfer based profiles
    heights_rad = np.arange(0,max_height+layer_thickness,layer_thickness)
    LAD_rad = np.zeros((n_subplots,heights_rad.size,max_return))
    LAD_rad_DTM = np.zeros((n_subplots,heights_rad.size,max_return))

    # set up some arrays to host the MacArthur-Horn profiles
    #heights = np.arange(0,max_height)+1
    heights = np.arange(0,max_height,layer_thickness)+layer_thickness
    LAD_MH = np.zeros((n_subplots, heights.size))
    
    # set up array to host the lidar return profiles
    lidar_return_profiles = np.zeros((n_subplots, heights_rad.size, max_return))
    lidar_return_profiles_adj = np.zeros((n_subplots, heights_rad.size, max_return))

    # set up array to host inventory profiles
    field_LAD_profiles = np.zeros((n_subplots,heights.size))

    # set up array to hold the hemisfer LAI estimate
    LAI_hemisfer = np.zeros(n_subplots)


    # Also get estimate of lower canopy volume from small SAFE plots
    """
    if Plot_name in ['Belian','Seraya','DC1','DC2']:
        block = 'VJR'
    elif Plot_name in ['B North', 'B South']:
        block = 'B'
    elif Plot_name == 'E':
        block = 'E'
    elif Plot_name == 'LF':
        block = 'LF'
    else:
        block = '?'

    Ht,Area,Depth,StemDensity = field.calculate_crown_dimensions_small_plots(stem_data[block]['dbh'],stem_data[block]['height'],stem_data[block]['stem_density'],a_ht, b_ht, CF_ht, a_A, b_A, CF_A, a, b, CF)
    smallstem_LAD_profile = field.calculate_LAD_profiles_from_stem_size_distributions(heights, Area, Depth, Ht, StemDensity, beta)
    """

    pt_count = 0.
    # loop through subplots, calculating both return profiles and LAD distributions
    for i in range(0,n_subplots):
        #print "Subplot: ", subplot_labels[Plot_name][i]
        subplot_index = subplot_labels[Plot_name][i]-1
        # filter lidar points into subplot
        sp_pts = lidar.filter_lidar_data_by_polygon(plot_lidar_pts,subplot_polygons[Plot_name][i,:,:])
        pt_count += sp_pts.shape[0]
        # first of all, loop through the return numbers to calculate the radiative LAD profiles
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
        lidar_return_profiles_adj[subplot_index,:,:n.shape[2]] = np.sum(n.copy(),axis=1)

        # now get MacArthur-Horn profiles
        heights,first_return_profile,n_ground_returns = LAD1.bin_returns(sp_pts, max_height, layer_thickness)
        LAD_MH[subplot_index,:] = LAD1.estimate_LAD_MacArthurHorn(first_return_profile, n_ground_returns, layer_thickness, kappa)

        # now get field inventory estimate
        mask = np.all((field_data['plot']==Plot_name,field_data['subplot']==subplot_labels[Plot_name][i]),axis=0)
        Ht,Area,Depth = field.calculate_crown_dimensions(field_data['DBH_field'][mask],field_data['Height_field'][mask],field_data['CrownArea'][mask], a_ht, b_ht, CF_ht, a_A, b_A, CF_A, a, b, CF)
        field_LAD_profiles[subplot_index,:], CanopyV = field.calculate_LAD_profiles_generic(heights, Area, Depth, Ht, beta, subplot_area)


        if Plot_name == 'DC1':
            Ht,Area,Depth,StemDensity = field.calculate_crown_dimensions_for_stem_distributions(DC1_stem_data['dbh'],DC1_stem_data['stem_density'][:,subplot_index],a_ht, b_ht, CF_ht, a_A, b_A, CF_A, a, b, CF)
        elif Plot_name == 'DC2':
            Ht,Area,Depth,StemDensity = field.calculate_crown_dimensions_for_stem_distributions(DC2_stem_data['dbh'],DC2_stem_data['stem_density'][:,subplot_index],a_ht, b_ht, CF_ht, a_A, b_A, CF_A, a, b, CF)
        else:
            Ht,Area,Depth,StemDensity = field.calculate_crown_dimensions_for_stem_distributions(SAFE_stem_data[Plot_name]['dbh'],SAFE_stem_data[Plot_name]['stem_density'][:,subplot_index],a_ht, b_ht, CF_ht, a_A, b_A, CF_A, a, b, CF)

        smallstem_LAD_profile = field.calculate_LAD_profiles_from_stem_size_distributions(heights, Area, Depth, Ht, StemDensity, beta)
        field_LAD_profiles[subplot_index,:]+=smallstem_LAD_profile
        # now load in the LAI estimates from the hemispherical photographs
        Hemisfer_mask = np.all((field_LAI['Subplot']==subplot_labels[Plot_name][i],field_LAI['Plot']==Plot_name),axis=0)
        LAI_hemisfer[subplot_index] = field_LAI['LAI'][Hemisfer_mask]
    
    print "average point density = ", pt_count/10.**4, " pts/m^2"
    # now we have looped through and created the different profiles, need to account for any NaN's and apply minimum height
    # to the LAD distributions
    # - set NaN values to zero
    LAD_rad[np.isnan(LAD_rad)]=0
    LAD_rad_DTM[np.isnan(LAD_rad_DTM)]=0
    LAD_MH[np.isnan(LAD_MH)]=0
    # - remove all profile values below minimum height prior to comparison
    mask = heights <= minimum_height
    LAD_MH[:,mask]=0
    mask = np.max(heights_rad)-heights_rad<=minimum_height
    LAD_rad[:,mask]=0
    LAD_rad_DTM[:,mask]=0

    # now the profiles are ready to be stored into their relevant dictionaries
    lidar_profiles[Plot_name] = lidar_return_profiles.copy()
    lidar_profiles_adjusted[Plot_name] = lidar_return_profiles_adj.copy()

    MacArthurHorn_LAD[Plot_name] = LAD_MH.copy()
    radiative_LAD[Plot_name] = LAD_rad.copy()
    radiative_DTM_LAD[Plot_name] = LAD_rad_DTM.copy()
    inventory_LAD[Plot_name] = field_LAD_profiles.copy()

    # Also ready to store their respective LAI profiles
    MacArthurHorn_LAI[Plot_name] = np.sum(LAD_MH,axis=1)*layer_thickness
    radiative_LAI[Plot_name] = np.sum(LAD_rad,axis=1)*layer_thickness
    radiative_DTM_LAI[Plot_name] = np.sum(LAD_rad_DTM,axis=1)*layer_thickness
    inventory_LAI[Plot_name] = np.sum(field_LAD_profiles,axis=1)*layer_thickness
    Hemisfer_LAI[Plot_name] = LAI_hemisfer.copy()

#########################################################################################
# Now plot some figures
#----------------------------------------------------------------------------------------
# Some useful arrays for getting stats (dictionaries not so good for this)
hemisfer = np.zeros(N_plots)
mh = np.zeros(N_plots)
rad2 = np.zeros(N_plots)
rad3 = np.zeros(N_plots)
radDTM2 = np.zeros(N_plots)
radDTM3 = np.zeros(N_plots)
vol = np.zeros(N_plots)
for i in range(0,N_plots):
    hemisfer[i]=np.mean(Hemisfer_LAI[Plots[i]])
    mh[i]=np.mean(MacArthurHorn_LAI[Plots[i]])
    rad2[i] = np.mean(radiative_LAI[Plots[i]][:,-2])
    radDTM2[i] = np.mean(radiative_DTM_LAI[Plots[i]][:,-2])
    rad3[i] = np.mean(radiative_LAI[Plots[i]][:,-1])
    radDTM3[i] = np.mean(radiative_DTM_LAI[Plots[i]][:,-1])
    vol[i]=np.mean(inventory_LAI[Plots[i]])

#----------------------------------------------------------------------------------------
# Figure 1: Sample point clouds for three sites: Belian, LF, B North

# first up, going to need to find affine transformation to rotate the point cloud for
# easy plotting

# The GPS coordinates for the plot locations can be used to find this transformation matrix
gps_pts_file = 'GPS_points_file_for_least_squares_fitting.csv'
datatype = {'names': ('plot', 'x', 'y', 'x_prime', 'y_prime'), 'formats': ('S32','f16','f16','f16','f16')}
plot_coords = np.genfromtxt(gps_pts_file, delimiter = ',',dtype=datatype)
plot_coords['plot'][plot_coords['plot']=='danum_1'] = 'DC1'
plot_coords['plot'][plot_coords['plot']=='danum_2'] = 'DC2'
plot_coords['plot'][plot_coords['plot']=='maliau_belian'] = 'Belian'
plot_coords['plot'][plot_coords['plot']=='maliau_seraya'] = 'Seraya'
plot_coords['plot'][plot_coords['plot']=='B_north'] = 'B North'
plot_coords['plot'][plot_coords['plot']=='B_south'] = 'B South'
plot_coords['plot'][plot_coords['plot']=='LFE'] = 'LF'

affine = {}
for pp in range(0, N_plots):
    # first get points for a given plot and build matrices - note that I've reversed xy and xy_prime in this case to reverse the rotation-translation
    mask = plot_coords['plot']==Plots[pp]
    x = plot_coords['x_prime'][mask]
    y = plot_coords['y_prime'][mask]
    x_prime = plot_coords['x'][mask]
    y_prime = plot_coords['y'][mask]
    affine[Plots[pp]]=lstsq.least_squares_affine_matrix(x,y,x_prime,y_prime)

    Xi = np.asarray([plot_point_cloud[Plots[pp]][:,0],plot_point_cloud[Plots[pp]][:,1],np.ones(plot_point_cloud[Plots[pp]].shape[0])])
    Xi_prime = np.dot(affine[Plots[pp]],Xi)

    plot_point_cloud[Plots[pp]][:,0]=Xi_prime[0]
    plot_point_cloud[Plots[pp]][:,1]=Xi_prime[1]

plt.figure(1, facecolor='White',figsize=[8,12])

colour = ['#46E900','#1A2BCE','#E0007F']
rgb = [[70,233,0],[26,43,206],[224,0,127]]
labels = ['$1^{st}$', '$2^{nd}$', '$3^{rd}$', '$4^{th}$']

# Belian
ax1a = plt.subplot2grid((3,1),(0,0))
ax1a.annotate('a - Maliau Reserve, MAO01', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax1a.set_ylabel('Height / m',fontsize=axis_size)
plt.gca().set_aspect('equal', adjustable='box-forced')

# LF
ax1b = plt.subplot2grid((3,1),(1,0),sharey=ax1a,sharex=ax1a)
ax1b.annotate('b - SAFE landscape, SAF04', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax1b.set_ylabel('Height / m',fontsize=axis_size)
plt.gca().set_aspect('equal', adjustable='box-forced')

# B North
ax1c = plt.subplot2grid((3,1),(2,0),sharey=ax1a,sharex=ax1a)
ax1c.annotate('c - SAFE landscape, SAF02', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax1c.set_ylabel('Height / m',fontsize=axis_size)
ax1c.set_xlabel('Horizontal distance / m',fontsize=axis_size)
plt.gca().set_aspect('equal', adjustable='box-forced')

fig1_plots = ['Belian', 'LF', 'B North']
axes = [ax1a, ax1b, ax1c]
for pp in range(0,3):
    plot_lidar_pts = plot_point_cloud[fig1_plots[pp]]
    for k in range(0,max_return):
    
        mask = np.all((plot_lidar_pts[:,0]>=0,plot_lidar_pts[:,0]<=100,plot_lidar_pts[:,1]>=0,plot_lidar_pts[:,0]<=100,plot_lidar_pts[:,3]==k+1),axis=0)
        points_x = 100-plot_lidar_pts[mask][:,0]#-np.min(plot_lidar_pts[:,0])
        points_z = plot_lidar_pts[mask][:,2]
        points_y = plot_lidar_pts[mask][:,1]#-np.min(plot_lidar_pts[:,1])
        
        alpha_max = 0.2
        colours = np.zeros((points_x.size,4))
        colours[:,0]=rgb[k][0]/255.
        colours[:,1]=rgb[k][1]/255.
        colours[:,2]=rgb[k][2]/255.
        colours[:,3]=alpha_max*(1-points_x/(points_x.max()+1))
        axes[pp].scatter(points_y,points_z,marker='o',c=colours,edgecolors='none',s=2)#,label=labels[k])
        axes[pp].scatter(0,0,marker='o',c=colours[0,0:3],edgecolors='none',s=2,label=labels[k])

ax1a.set_ylim(0,80)
ax1a.set_xlim(0,100)
ax1a.legend(loc=1,fontsize=axis_size)
plt.tight_layout()
plt.savefig(output_dir+'fig1_plot_pointclouds.png')


#-----------------------------------------------------------------------------------------
# Figure 2: Transmission ratios
# This figure plots the number of vegetation returns that propagate through further into
# the canopy, defined as the transmittance ratio.
 
# A figure illustrating transmittance ratio between successive returns 
plt.figure(2, facecolor='White',figsize=[3,3])
ax2 = plt.subplot2grid((1,1),(0,0))
ax2.set_xlabel('return number')
ax2.set_ylabel('transmittance ratio')
for i in range(0,4):
    if i==0:
        ax2.plot(1,1,'o',color='blue')
    else:
        N_veg = float(np.all((all_lidar_pts[:,3]==i,all_lidar_pts[:,4]==1),axis=0).sum())
        N_i = float((all_lidar_pts[:,3]==i+1).sum())
        ax2.plot(i+1,N_i/N_veg,'o',color='#1A2BCE')

ax2.set_ylim(0,1.1)
ax2.set_xlim(0.2,4.8)
plt.tight_layout()
plt.savefig(output_dir+'fig2_transmittance_ratios.png')

#-----------------------------------------------------------------------------------------
# Figure 3- Allometric relationships used to construct the inventory-based LAD profiles.
plt.figure(3,facecolor='White',figsize=[9,3])
ax3a = plt.subplot2grid((1,3),(0,0))
ax3a.set_ylabel('Crown Depth / m')
ax3a.set_xlabel('Height / m')
ax3a.annotate('a - Height-Crown Depth', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax3a.plot(H,D,'.',color='#1A2BCE',alpha=0.1)
H_mod = np.linspace(0.,H.max(),1000)
D_mod = CF*a*H_mod**b
ax3a.plot(H_mod,D_mod,'-',color='black')
eq = '$D=%.3fH^{%.3f}$' % (CF*a, b)
ax3a.annotate(eq, xy=(0.95,0.05), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='bottom', fontsize=10)

ax3b = plt.subplot2grid((1,3),(0,1))
ax3b.annotate('b - DBH-Crown Area', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax3b.set_xlabel('DBH / cm')
ax3b.set_ylabel('Crown Area / m$^2$')
ax3b.plot(field_data['DBH_field'],field_data['CrownArea'],'.',color='#1A2BCE',alpha=0.1)
DBH_mod = np.linspace(0.,np.nanmax(field_data['DBH_field']),1000)
CA_mod = CF_A*a_A*DBH_mod**b_A
ax3b.plot(DBH_mod,CA_mod,'-',color='black')
eq = '$A=%.3fDBH^{%.3f}$' % (CF_A*a_A, b_A)
ax3b.annotate(eq, xy=(0.95,0.05), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='bottom', fontsize=10)

ax3c = plt.subplot2grid((1,3),(0,2),sharex=ax3b)
ax3c.annotate('c - DBH-Height', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax3c.set_xlabel('DBH / cm')
ax3c.set_ylabel('Height / m')
ax3c.plot(field_data['DBH_field'],field_data['Height_field'],'.',color='#1A2BCE',alpha=0.1)
Ht_mod = CF_ht*a_ht*DBH_mod**b_ht
ax3c.plot(DBH_mod,Ht_mod,'-',color='black')
eq = '$H=%.3fDBH^{%.3f}$' % (CF_ht*a_ht, b_ht)
ax3c.annotate(eq, xy=(0.95,0.05), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='bottom', fontsize=10)


plt.tight_layout()
plt.savefig(output_dir+'fig4_allometric_relationships.png')


#-----------------------------------------------------------------------------------------
# Figure 4: Canopy profile comparisons across the gradient
plt.figure(4, facecolor='White',figsize=[8,12])
# Belian
# - returns
ax4a = plt.subplot2grid((3,5),(0,0))
ax4a.annotate('a', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
#ax4a.annotate('MAO01, OG', xy=(0.95,0.89), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)
ax4a.annotate('LiDAR returns', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=9)
ax4a.set_ylabel('Height / m',fontsize=axis_size)
# - MacHorn
ax4b = plt.subplot2grid((3,5),(0,1),sharey=ax4a)
ax4b.annotate('b', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax4b.annotate('MacArthur-Horn', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=9)
# - Detto
ax4c = plt.subplot2grid((3,5),(0,2),sharey=ax4a,sharex=ax4b)
ax4c.annotate('c', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax4c.annotate('rad. trans.\n(Detto)', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=9)
ax4c.set_title('Old-growth forest, MAO01', fontsize=10)
# - Corrected rad trans
ax4d = plt.subplot2grid((3,5),(0,3),sharey=ax4a,sharex=ax4b)
ax4d.annotate('d', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax4d.annotate('rad. trans.\n(new)', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=9)
# - Inventory
ax4e = plt.subplot2grid((3,5),(0,4),sharey=ax4a,sharex=ax4b)
ax4e.annotate('e', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax4e.annotate('crown volume', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=9)

# LF
# - returns
ax4f = plt.subplot2grid((3,5),(1,0), sharex = ax4a, sharey = ax4a)
ax4f.annotate('f', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
#ax4f.annotate('SAF04, ML', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)
ax4f.set_ylabel('Height / m',fontsize=axis_size)
# - MacHorn
ax4g = plt.subplot2grid((3,5),(1,1),sharey=ax4a, sharex=ax4b)
ax4g.annotate('g', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
# - Detto
ax4h = plt.subplot2grid((3,5),(1,2),sharey=ax4a,sharex=ax4b)
ax4h.annotate('h', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax4h.set_title('Moderately logged forest, SAF04', fontsize=10)
# - Corrected rad trans
ax4i = plt.subplot2grid((3,5),(1,3),sharey=ax4a,sharex=ax4b)
ax4i.annotate('i', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
# - Inventory
ax4j = plt.subplot2grid((3,5),(1,4),sharey=ax4a,sharex=ax4b)
ax4j.annotate('j', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)


# B North
# - returns
ax4k = plt.subplot2grid((3,5),(2,0), sharex = ax4a, sharey = ax4a)
ax4k.annotate('k', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
#ax4k.annotate('SAF02, HL', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)
ax4k.set_ylabel('Height / m',fontsize=axis_size)
ax4k.set_xlabel('Number of returns\n(x1000)',fontsize=axis_size,horizontalalignment='center')
# - MacHorn
ax4l = plt.subplot2grid((3,5),(2,1),sharey=ax4a, sharex=ax4b)
ax4l.annotate('l', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax4l.set_xlabel('PAD\n(m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size,horizontalalignment='center')

# - Detto
ax4m = plt.subplot2grid((3,5),(2,2),sharey=ax4a,sharex=ax4b)
ax4m.annotate('m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax4m.set_xlabel('PAD\n(m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size,horizontalalignment='center')
ax4m.set_title('Heavily logged forest, SAF02', fontsize=10)
# - Corrected rad trans
ax4n = plt.subplot2grid((3,5),(2,3),sharey=ax4a,sharex=ax4b)
ax4n.annotate('n', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax4n.set_xlabel('PAD\n(m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size,horizontalalignment='center')
# - Inventory
ax4o = plt.subplot2grid((3,5),(2,4),sharey=ax4a,sharex=ax4b)
ax4o.annotate('o', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax4o.set_xlabel('Crown Volume\n(m$^3$m$^{-2}$m$^{-1}$)',fontsize=axis_size,horizontalalignment='center')

fig2_plots = ['Belian', 'LF', 'B North']
axes1 = [ax4a, ax4f, ax4k]
axes2 = [ax4b, ax4g, ax4l]
axes3 = [ax4c, ax4h, ax4m]
axes4 = [ax4d, ax4i, ax4n]
axes5 = [ax4e, ax4j, ax4o]
for pp in range(0,3):
    Plot_name = fig2_plots[pp]
    
    # plot lidar profile
    return_dist_adj = np.sum(lidar_profiles_adjusted[Plot_name],axis=0)
    return_dist     = np.sum(lidar_profiles[Plot_name],axis=0)
    labels = ['$1^{st}$', '$2^{nd}$', '$3^{rd}$', '$4^{th}$']
    for k in range(0,max_return):
        axes1[pp].plot(return_dist[:,k]/1000.,np.max(heights_rad)-heights_rad,'-',c=colour[k],linewidth=1,label=labels[k])

    # plot macarthur horn profile
    for i in range(0,n_subplots):
        axes2[pp].fill_betweenx(heights[2:],0,MacArthurHorn_LAD[Plot_name][i,2:],color=colour[0],alpha=0.05)
    axes2[pp].plot(np.mean(MacArthurHorn_LAD[Plot_name],axis=0)[2:],heights[2:],'-',c=colour[0],linewidth=2)

    # plot detto profile
    for i in range(0,n_subplots):
        axes3[pp].fill_betweenx(heights_rad[3:],0,radiative_LAD[Plot_name][i,:-3,-1][::-1],color=colour[1],alpha=0.05)
    axes3[pp].plot(np.mean(radiative_LAD[Plot_name][:,:-3,1],axis=0)[::-1],heights_rad[3:],'-',c=colour[1],linewidth=2)

    # plot corrective radiative transfer profile
    for i in range(0,n_subplots):
        axes4[pp].fill_betweenx(heights_rad[3:],0,radiative_DTM_LAD[Plot_name][i,:-3,-1][::-1],color=colour[1],alpha=0.05)
    axes4[pp].plot(np.mean(radiative_DTM_LAD[Plot_name][:,:-3,1],axis=0)[::-1],heights_rad[3:],'-',c=colour[1],linewidth=2)

    # field inventory
    for i in range(0,n_subplots):
        axes5[pp].fill_betweenx(heights[2:],0,inventory_LAD[Plot_name][i,2:],color=colour[2],alpha=0.05)
    axes5[pp].plot(np.mean(inventory_LAD[Plot_name],axis=0)[2:],heights[2:],'-',c=colour[2],linewidth=2)



ax4a.set_ylim(0,80)
ax4a.set_xlim(0,29)
ax4a.legend(loc='lower right')
ax4b.set_xlim(xmin=0,xmax=0.7)

ax4l.locator_params(axis='x',nbins=5)
ax4m.locator_params(axis='x',nbins=5)
ax4n.locator_params(axis='x',nbins=5)
ax4o.locator_params(axis='x',nbins=5)
ax4k.locator_params(axis='x',nbins=5)


yticklabels = ax4b.get_yticklabels() + ax4c.get_yticklabels() + ax4d.get_yticklabels() + ax4e.get_yticklabels() + ax4g.get_yticklabels() + ax4h.get_yticklabels() + ax4i.get_yticklabels() + ax4j.get_yticklabels() + ax4l.get_yticklabels() + ax4m.get_yticklabels() + ax4n.get_yticklabels() + ax4o.get_yticklabels()
xticklabels = ax4a.get_xticklabels() + ax4b.get_xticklabels() + ax4c.get_xticklabels() + ax4d.get_xticklabels() + ax4e.get_xticklabels() + ax4f.get_xticklabels() + ax4g.get_xticklabels() + ax4h.get_xticklabels() + ax4i.get_xticklabels() + ax4j.get_xticklabels() 

plt.setp(yticklabels,visible=False)
plt.setp(xticklabels,visible=False)
plt.subplots_adjust(hspace=0.2, wspace = 0.1)

#plt.tight_layout()
plt.savefig(output_dir+'fig5_plot_LAD_profiles.png')
#plt.show()

#--------------------------------------------------------------------------------------
# Figure 5: comparison of PAI
# annotate with stats
r_sq_a = [aux.get_rsquared_annotation(mh,rad2), aux.get_rsquared_annotation(mh,rad3)]
r_sq_b = [aux.get_rsquared_annotation(mh,radDTM2), aux.get_rsquared_annotation(mh,radDTM3)]

plt.figure(5, facecolor='White',figsize=[7,4])
ax5a = plt.subplot2grid((1,2),(0,0))
ax5a.set_xlabel('PAI$_{MacArthur-Horn}$',fontsize=axis_size)
ax5a.set_ylabel('PAI$_{rad}$',fontsize=axis_size)
ax5a.annotate('a - radiative transfer (Detto)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax5a.annotate('$k_{max}=2$; ' + r_sq_a[0] + '\n$k_{max}=3$; ' + r_sq_a[1], xy=(0.95,0.90), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)

ax5a.plot([0,20],[0,20],'--',color='black',alpha=0.3)
for i in range(0,N_plots):
    for k in range(1,max_return):
        ax5a.plot(MacArthurHorn_LAI[Plots[i]],radiative_LAI[Plots[i]][:,k],'.',color=colour[k],alpha=0.5)
    
for i in range(0,N_plots):
    for k in range(1,max_return):
        x_err=np.std(MacArthurHorn_LAI[Plots[i]])/np.sqrt(n_subplots)
        y_err=np.std(radiative_LAI[Plots[i]][:,k])/np.sqrt(n_subplots)
        ax5a.errorbar(np.mean(MacArthurHorn_LAI[Plots[i]]),np.mean(radiative_LAI[Plots[i]][:,k]),xerr=x_err,yerr=y_err,marker='o',color=colour[k])

ax5b = plt.subplot2grid((1,2),(0,1), sharex = ax5a, sharey = ax5a)
ax5b.set_xlabel('PAI$_{MacArthur-Horn}$',fontsize=axis_size)
ax5b.set_ylabel('PAI$_{rad}$',fontsize=axis_size)
ax5b.annotate('b - radiative transfer (modified)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax5b.annotate('$k_{max}=2$; ' + r_sq_b[0] + '\n$k_{max}=3$; ' + r_sq_b[1], xy=(0.95,0.90), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)

ax5b.plot([0,20],[0,20],'--',color='black',alpha=0.3)
for i in range(0,N_plots):
    for k in range(1,max_return):
        ax5b.plot(MacArthurHorn_LAI[Plots[i]],radiative_DTM_LAI[Plots[i]][:,k],'.',color=colour[k],alpha=0.5)
    
for i in range(0,N_plots):
    for k in range(1,max_return):
        x_err=np.std(MacArthurHorn_LAI[Plots[i]])/np.sqrt(n_subplots)
        y_err=np.std(radiative_DTM_LAI[Plots[i]][:,k])/np.sqrt(n_subplots)
        if i == 0:
            leg_label = '$k_{max}=$' + str(k+1) 
            ax5b.errorbar(np.mean(MacArthurHorn_LAI[Plots[i]]),np.mean(radiative_DTM_LAI[Plots[i]][:,k]),xerr=x_err,yerr=y_err,marker='o',color=colour[k],label = leg_label)
        else:
            ax5b.errorbar(np.mean(MacArthurHorn_LAI[Plots[i]]),np.mean(radiative_DTM_LAI[Plots[i]][:,k]),xerr=x_err,yerr=y_err,marker='o',color=colour[k])

# configs
ax5a.legend(loc=4)
ax5b.legend(loc=4)
ax5a.set_xlim((0,20))
ax5a.set_ylim((0,20))

plt.tight_layout()
plt.savefig(output_dir+'fig6_lidar_LAI_comparison.png')
#plt.show()

#--------------------------------------------------------------------------------------
# Figure 6: LAI vs. hemiphotos

# annotate with stats
r_sq_a = aux.get_rsquared_annotation(hemisfer,mh)
r_sq_b1 =  aux.get_rsquared_annotation(hemisfer,rad3)
r_sq_b2 = aux.get_rsquared_annotation(hemisfer,radDTM3)
r_sq_c = aux.get_rsquared_annotation(hemisfer,vol)

plt.figure(6, facecolor='White',figsize=[9,4])
ax6a = plt.subplot2grid((1,3),(0,0))
ax6a.set_xlabel('PAI$_{Hemisfer}$',fontsize=axis_size)
ax6a.set_ylabel('PAI$_{MacArthur-Horn}$',fontsize=axis_size)
ax6a.annotate('a - MacArthur-Horn', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax6a.annotate(r_sq_a, xy=(0.95,0.90), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)
ax6a.plot([0,20],[0,20],'--',color='black',alpha=0.3)

for i in range(0,N_plots):
    ax6a.plot(Hemisfer_LAI[Plots[i]],MacArthurHorn_LAI[Plots[i]],'.',color=colour[0],alpha=0.5)

for i in range(0,N_plots):
    x_err=np.std(Hemisfer_LAI[Plots[i]])/np.sqrt(n_subplots)
    y_err=np.std(MacArthurHorn_LAI[Plots[i]])/np.sqrt(n_subplots)
    ax6a.errorbar(np.mean(Hemisfer_LAI[Plots[i]]),np.mean(MacArthurHorn_LAI[Plots[i]]),x_err,y_err,'o',color='black')
#ax6a.plot(LAI_hemi_mod, LAI_MH_mod, '-', color = 'k')


ax6b = plt.subplot2grid((1,3),(0,1), sharex=ax6a, sharey=ax6a)
ax6b.annotate('b - radiative transfer', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax6b.annotate(r_sq_b2 + ' (modified)\n' + r_sq_b1 + ' (Detto)', xy=(0.95,0.90), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)
ax6b.set_xlabel('PAI$_{Hemisfer}$',fontsize=axis_size)
ax6b.set_ylabel('PAI$_{rad}$',fontsize=axis_size)
ax6b.plot([0,20],[0,20],'--',color='black',alpha=0.3)
for i in range(0,N_plots):
    ax6b.plot(Hemisfer_LAI[Plots[i]],radiative_LAI[Plots[i]][:,1],'.',color='0.5',alpha=0.5)
    ax6b.plot(Hemisfer_LAI[Plots[i]],radiative_DTM_LAI[Plots[i]][:,1],'.',color=colour[1],alpha=0.5)

for i in range(0,N_plots):
    x_err=np.std(Hemisfer_LAI[Plots[i]])/np.sqrt(n_subplots)
    y_err1=np.std(radiative_LAI[Plots[i]][:,1])/np.sqrt(n_subplots)
    y_err2=np.std(radiative_DTM_LAI[Plots[i]][:,1])/np.sqrt(n_subplots)
    ax6b.errorbar(np.mean(Hemisfer_LAI[Plots[i]]),np.mean(radiative_LAI[Plots[i]][:,1]),xerr=x_err,yerr=y_err1,marker='o',color='black',mfc='white')
    ax6b.errorbar(np.mean(Hemisfer_LAI[Plots[i]]),np.mean(radiative_DTM_LAI[Plots[i]][:,1]),xerr=x_err,yerr=y_err2,marker='o',color='black')

#ax6b.plot(LAI_hemi_mod, LAI_rad_mod, '-', color = 'k')


ax6c = plt.subplot2grid((1,3),(0,2), sharex=ax6a)
ax6c.annotate('c - field inventory', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax6c.annotate(r_sq_c, xy=(0.95,0.90), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)
ax6c.set_xlabel('PAI$_{Hemisfer}$',fontsize=axis_size)
ax6c.set_ylabel('crown volume / m$^3$m$^{-2}$',fontsize=axis_size)

for i in range(0,N_plots):
    ax6c.plot(Hemisfer_LAI[Plots[i]],inventory_LAI[Plots[i]],'.',color=colour[2],alpha=0.5)

for i in range(0,N_plots):
    x_err=np.std(Hemisfer_LAI[Plots[i]])/np.sqrt(n_subplots)
    y_err=np.std(inventory_LAI[Plots[i]])/np.sqrt(n_subplots)
    ax6c.errorbar(np.mean(Hemisfer_LAI[Plots[i]]),np.mean(inventory_LAI[Plots[i]]),xerr=x_err,yerr=y_err,marker='o',color='black')


ax6a.set_xlim((0,10))
ax6a.set_ylim((0,20))
ax6c.set_ylim(ymin=0)
plt.tight_layout()
plt.savefig(output_dir+'fig11_LAI_hemiphoto_comparison.png')
#plt.show()



#--------------------------------------------------------------------------------------
# Figure 7: LAI vs. canopy volume

# annotate with stats
r_sq_a = aux.get_rsquared_annotation(vol,mh)
r_sq_b=  [aux.get_rsquared_annotation(vol,rad2),aux.get_rsquared_annotation(vol,rad3)]
r_sq_c = [aux.get_rsquared_annotation(vol,radDTM2),aux.get_rsquared_annotation(vol,radDTM3)]

plt.figure(7, facecolor='White',figsize=[9,4])
ax7a = plt.subplot2grid((1,3),(0,0))
ax7a.set_xlabel('crown volume / m$^3$m$^{-2}$',fontsize=axis_size)
ax7a.set_ylabel('PAI',fontsize=axis_size)
ax7a.annotate('a - MacArthur-Horn', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax7a.annotate(r_sq_a, xy=(0.95,0.90), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)
for i in range(0,N_plots):
    ax7a.plot(inventory_LAI[Plots[i]],MacArthurHorn_LAI[Plots[i]],'.',color=colour[0],alpha=0.5)

for i in range(0,N_plots):
    x_err=np.std(inventory_LAI[Plots[i]])/np.sqrt(n_subplots)
    y_err=np.std(MacArthurHorn_LAI[Plots[i]])/np.sqrt(n_subplots)
    ax7a.errorbar(np.mean(inventory_LAI[Plots[i]]),np.mean(MacArthurHorn_LAI[Plots[i]]),xerr=x_err,yerr=y_err,marker='o',color='black')


ax7b = plt.subplot2grid((1,3),(0,1), sharex=ax7a, sharey=ax7a)
ax7b.annotate('b - radiative transfer (Detto)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax7b.annotate('$k_{max}=2$; ' + r_sq_b[0] + '\n$k_{max}=3$; ' + r_sq_b[1], xy=(0.95,0.90), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)
ax7b.set_xlabel('crown volume / m$^3$m$^{-2}$',fontsize=axis_size)
ax7b.set_ylabel('PAI',fontsize=axis_size)
for k in range(1,3):
    for i in range(0,N_plots):
        ax7b.plot(inventory_LAI[Plots[i]],radiative_LAI[Plots[i]][:,k],'.',color=colour[k],alpha=0.5)

for k in range(1,3):
    for i in range(0,N_plots):
        x_err=np.std(inventory_LAI[Plots[i]])/np.sqrt(n_subplots)
        y_err=np.std(radiative_LAI[Plots[i]][:,k])/np.sqrt(n_subplots)
        if i==0:
            leg_label = '$k_{max}=$' + str(k+1) 
            ax7b.errorbar(np.mean(inventory_LAI[Plots[i]]),np.mean(radiative_LAI[Plots[i]][:,k]),xerr=x_err,yerr=y_err,marker='o',color=colour[k],label=leg_label)
        else:
            ax7b.errorbar(np.mean(inventory_LAI[Plots[i]]),np.mean(radiative_LAI[Plots[i]][:,k]),xerr=x_err,yerr=y_err,marker='o',color=colour[k])

ax7c = plt.subplot2grid((1,3),(0,2), sharex=ax7a, sharey=ax7a)
ax7c.annotate('c - radiative transfer (modified)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax7c.annotate('$k_{max}=2$; ' + r_sq_c[0] + '\n$k_{max}=3$; ' + r_sq_c[1], xy=(0.95,0.90), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)
ax7c.set_ylabel('PAI',fontsize=axis_size)
ax7c.set_xlabel('crown volume / m$^3$m$^{-2}$',fontsize=axis_size)
for k in range(1,3):
    for i in range(0,N_plots):
        ax7c.plot(inventory_LAI[Plots[i]],radiative_DTM_LAI[Plots[i]][:,k],'.',color=colour[k],alpha=0.5)

for k in range(1,3):
    for i in range(0,N_plots):
        x_err=np.std(inventory_LAI[Plots[i]])/np.sqrt(n_subplots)
        y_err=np.std(radiative_DTM_LAI[Plots[i]][:,k])/np.sqrt(n_subplots)
        if i==0:
            leg_label = '$k_{max}=$' + str(k+1) 
            ax7c.errorbar(np.mean(inventory_LAI[Plots[i]]),np.mean(radiative_DTM_LAI[Plots[i]][:,k]),xerr=x_err,yerr=y_err,marker='o',color=colour[k],label=leg_label)
        else:
            ax7c.errorbar(np.mean(inventory_LAI[Plots[i]]),np.mean(radiative_DTM_LAI[Plots[i]][:,k]),xerr=x_err,yerr=y_err,marker='o',color=colour[k])
            
ax7c.legend(loc=4)

ax7a.set_ylim((0,20))
plt.tight_layout()
plt.savefig(output_dir+'fig10_LiDAR_LAI_canopy_volume_comparison.png')
#plt.show()
#--------------------------------------------------------------------------------------
#Figure 8:  Plot total basal area against LAI from the two preferred methods

census_file = '/home/dmilodow/DataStore_DTM/BALI/BALI_Cplot_data/SAFE_CarbonPlots_TreeCensus.csv'
census = cen.collate_plot_level_census_data(census_file)

BA = {}
Plots_SAFE = ['Belian', 'Seraya', 'LF', 'E','B North', 'B South']
Plots_Danum = ['DC1', 'DC2']
plot_colour = {}
plot_colour['Belian']=colour[0]
plot_colour['Seraya']=colour[0]
plot_colour['DC1']=colour[0]
plot_colour['DC2']=colour[0]
plot_colour['LF']=colour[1]
plot_colour['E']=colour[1]
plot_colour['B North']=colour[2]
plot_colour['B South']=colour[2]

plot_marker = {}
plot_marker['Belian']='o'
plot_marker['Seraya']='v'
plot_marker['DC1']='^'
plot_marker['DC2']='s'
plot_marker['LF']='o'
plot_marker['E']='v'
plot_marker['B North']='o'
plot_marker['B South']='v'
plot_label = {}
plot_label['Belian']='MAO01'
plot_label['Seraya']='MAO02'
plot_label['DC1']='DAN04'
plot_label['DC2']='DAN05'
plot_label['LF']='SAF04'
plot_label['E']='SAF03'
plot_label['B North']='SAF02'
plot_label['B South']='SAF01'

for pp in range(0,len(Plots_SAFE)):
    temp_BA = np.zeros(n_subplots)
    for ss in range(0,n_subplots):
        # check basal area
        temp_BA[ss] = census[Plots_SAFE[pp]]['BasalArea'][ss,0]*25/100**2
    temp_BA[np.isnan(temp_BA)]=0
    BA[Plots_SAFE[pp]]=temp_BA.copy()

for pp in range(0,len(Plots_Danum)):
    temp_BA = np.zeros(n_subplots)
    for ss in range(0,n_subplots):
        indices = np.all((field_data['plot']==Plots_Danum[pp],field_data['subplot']==ss+1,np.isfinite(field_data['DBH'])),axis=0)
        temp_BA[ss] = np.sum((field_data['DBH'][indices]/2)**2*np.pi)*25./100**2
    BA[Plots_Danum[pp]]=temp_BA.copy()

plt.figure(8, facecolor='White',figsize=[7,4])
ax8a = plt.subplot2grid((1,2),(0,0))
ax8a.set_xlabel('basal area / m$^2$ha$^{-1}$',fontsize=axis_size)
ax8a.set_ylabel('PAI',fontsize=axis_size)
ax8a.annotate('a - MacArthur-Horn', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

for i in range(0,N_plots):
    ax8a.plot(BA[Plots[i]],MacArthurHorn_LAI[Plots[i]],'.',color=plot_colour[Plots[i]],alpha=0.5)

for i in range(0,N_plots):
    x_err=np.std(BA[Plots[i]])/np.sqrt(n_subplots)
    y_err=np.std(MacArthurHorn_LAI[Plots[i]])/np.sqrt(n_subplots)
    ax8a.errorbar(np.mean(BA[Plots[i]]),np.mean(MacArthurHorn_LAI[Plots[i]]),xerr=x_err,yerr=y_err,marker=plot_marker[Plots[i]],markerfacecolor=plot_colour[Plots[i]], markeredgecolor='black',ecolor='black',linestyle='None')


ax8b = plt.subplot2grid((1,2),(0,1), sharex=ax8a, sharey=ax8a)
ax8b.annotate('b - radiative transfer (modified)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax8b.set_ylabel('PAI',fontsize=axis_size)
ax8b.set_xlabel('basal area / m$^2$ha$^{-1}$',fontsize=axis_size)

k = 1
for i in range(0,N_plots):
    ax8b.plot(BA[Plots[i]],radiative_DTM_LAI[Plots[i]][:,k],'.',color=plot_colour[Plots[i]],alpha=0.5)

for i in range(0,N_plots):
    x_err=np.std(BA[Plots[i]])/np.sqrt(n_subplots)
    y_err=np.std(radiative_DTM_LAI[Plots[i]][:,k])/np.sqrt(n_subplots)
    ax8b.errorbar(np.mean(BA[Plots[i]]),np.mean(radiative_DTM_LAI[Plots[i]][:,k]),xerr=x_err,yerr=y_err,marker=plot_marker[Plots[i]],markerfacecolor=plot_colour[Plots[i]], markeredgecolor='black',ecolor='black',linestyle='None',label = plot_label[Plots[i]])

ax8a.set_ylim((0,20))
ax8b.legend(loc='lower right',fontsize=axis_size-2)

# get R-sq for 1 ha plots
basalarea = np.zeros(N_plots)
MH_LAI = np.zeros(N_plots)
radDTM_LAI = np.zeros(N_plots)
for i in range(0,N_plots):
    basalarea[i]=np.mean(BA[Plots[i]])
    MH_LAI[i] = np.mean(MacArthurHorn_LAI[Plots[i]])
    radDTM_LAI[i] = np.mean(radiative_DTM_LAI[Plots[i]][:,-1])

# annotate with stats
r_sq_a = aux.get_rsquared_annotation(basalarea,MH_LAI)
r_sq_b =  aux.get_rsquared_annotation(basalarea,radDTM_LAI)
ax8a.annotate(r_sq_a, xy=(0.95,0.90), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)
ax8b.annotate(r_sq_b, xy=(0.95,0.90), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)

plt.tight_layout()
plt.savefig(output_dir+'fig12_LiDAR_LAI_BasalArea_comparison.png')
#plt.show()
#--------------------------------------------------------------------------------------
# Now put together a table that has all the LAI estimates
f = open(output_dir+'BALI_LAI_table.csv',"w") #opens file
f.write("Plot, BasalArea, error, CrownVolume, error, PAI_MH, error, PAI_Detto, error,  PAI_Detto, error, PAI_new, error, PAI_new, error, PAI_hemisfer\n")
for i in range(0,N_plots):
    f.write(Plots[i]+", ")
    f.write(str(np.mean(BA[Plots[i]]))+", "+str(np.std(BA[Plots[i]])/np.sqrt(n_subplots)) + ", ")
    f.write(str(np.mean(inventory_LAI[Plots[i]]))+", "+str(np.std(inventory_LAI[Plots[i]])/np.sqrt(n_subplots)) + ", ")
    f.write(str(np.mean(MacArthurHorn_LAI[Plots[i]]))+", "+str(np.std(MacArthurHorn_LAI[Plots[i]])/np.sqrt(n_subplots)) + ", ")
    f.write(str(np.mean(radiative_LAI[Plots[i]][:,1]))+", "+str(np.std(radiative_LAI[Plots[i]][:,1])/np.sqrt(n_subplots)) + ", ")
    f.write(str(np.mean(radiative_LAI[Plots[i]][:,2]))+", "+str(np.std(radiative_LAI[Plots[i]][:,2])/np.sqrt(n_subplots)) + ", ")
    f.write(str(np.mean(radiative_DTM_LAI[Plots[i]][:,1]))+", "+str(np.std(radiative_DTM_LAI[Plots[i]][:,1])/np.sqrt(n_subplots)) + ", ")
    f.write(str(np.mean(radiative_DTM_LAI[Plots[i]][:,2]))+", "+str(np.std(radiative_DTM_LAI[Plots[i]][:,2])/np.sqrt(n_subplots)) + ", ")
    f.write(str(np.mean(Hemisfer_LAI[Plots[i]]))+", "+str(np.std(Hemisfer_LAI[Plots[i]])/np.sqrt(n_subplots)) + "\n")
f.close()



#--------------------------------------------------------------------------------------
# Supplementary figures
figS_plots = ['Seraya', 'DC1', 'DC2','E', 'B South']
figS_codes = ['MAO02', 'DAN04', 'DAN05', 'SAF03', 'SAF01']
for pp in range(0,5):
    Plot_name = figS_plots[pp]
    plt.figure(9, facecolor='White',figsize=[8,12])
    axSa = plt.subplot2grid((3,5),(0,0),rowspan=2, colspan=5)
    axSa.annotate('a - ' + figS_codes[pp], xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    axSa.set_ylabel('Height / m',fontsize=axis_size)
    axSa.set_xlabel('Horizontal distance / m',fontsize=axis_size)
    plt.gca().set_aspect('equal', adjustable='box-forced')
    plot_lidar_pts = plot_point_cloud[figS_plots[pp]]
    for k in range(0,max_return):
    
        mask = np.all((plot_lidar_pts[:,0]>=0,plot_lidar_pts[:,0]<=100,plot_lidar_pts[:,1]>=0,plot_lidar_pts[:,0]<=100,plot_lidar_pts[:,3]==k+1),axis=0)
        points_x = 100-plot_lidar_pts[mask][:,0]
        points_z = plot_lidar_pts[mask][:,2]
        points_y = plot_lidar_pts[mask][:,1]
        
        alpha_max = 0.2
        colours = np.zeros((points_x.size,4))
        colours[:,0]=rgb[k][0]/255.
        colours[:,1]=rgb[k][1]/255.
        colours[:,2]=rgb[k][2]/255.
        colours[:,3]=alpha_max*(1-points_x/(points_x.max()+1))
        axSa.scatter(points_y,points_z,marker='o',c=colours,edgecolors='none',s=2)
        axSa.scatter(0,0,marker='o',c=colours[0,0:3],edgecolors='none',s=2,label=labels[k])

    axSa.set_ylim(0,80)
    axSa.set_xlim(0,100)
    axSa.legend(loc=1,fontsize=axis_size)

    axSb = plt.subplot2grid((3,5),(2,0), sharey = axSa)
    axSb.annotate('b', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    axSb.set_ylabel('Height / m',fontsize=axis_size)
    axSb.set_xlabel('Number of returns\n(x1000)',fontsize=axis_size,horizontalalignment='center')
    # - MacHorn
    axSc = plt.subplot2grid((3,5),(2,1),sharey=axSa)
    axSc.annotate('c', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    axSc.set_xlabel('PAD\n(m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size,horizontalalignment='center')

    # - Detto
    axSd = plt.subplot2grid((3,5),(2,2),sharey=axSa,sharex=axSc)
    axSd.annotate('d', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    axSd.set_xlabel('PAD\n(m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size,horizontalalignment='center')
    # - Corrected rad trans
    axSe = plt.subplot2grid((3,5),(2,3),sharey=axSa,sharex=axSc)
    axSe.annotate('e', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    axSe.set_xlabel('PAD\n(m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size,horizontalalignment='center')
    # - Inventory
    axSf = plt.subplot2grid((3,5),(2,4),sharey=axSa,sharex=axSc)
    axSf.annotate('f', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    axSf.set_xlabel('Crown Volume\n(m$^3$m$^{-2}$m$^{-1}$)',fontsize=axis_size,horizontalalignment='center')
    
    # plot lidar profile
    return_dist_adj = np.sum(lidar_profiles_adjusted[Plot_name],axis=0)
    return_dist     = np.sum(lidar_profiles[Plot_name],axis=0)
    labels = ['$1^{st}$', '$2^{nd}$', '$3^{rd}$', '$4^{th}$']
    for k in range(0,max_return):
        axSb.plot(return_dist[:,k]/1000.,np.max(heights_rad)-heights_rad,'-',c=colour[k],linewidth=1,label=labels[k])

    # plot macarthur horn profile
    for i in range(0,n_subplots):
        axSc.fill_betweenx(heights[2:],0,MacArthurHorn_LAD[Plot_name][i,2:],color=colour[0],alpha=0.05)
    axSc.plot(np.mean(MacArthurHorn_LAD[Plot_name],axis=0)[2:],heights[2:],'-',c=colour[0],linewidth=2)

    # plot detto profile
    for i in range(0,n_subplots):
        axSd.fill_betweenx(heights_rad[3:],0,radiative_LAD[Plot_name][i,:-3,-1][::-1],color=colour[1],alpha=0.05)
    axSd.plot(np.mean(radiative_LAD[Plot_name][:,:-3,1],axis=0)[::-1],heights_rad[3:],'-',c=colour[1],linewidth=2)

    # plot corrective radiative transfer profile
    for i in range(0,n_subplots):
        axSe.fill_betweenx(heights_rad[3:],0,radiative_DTM_LAD[Plot_name][i,:-3,-1][::-1],color=colour[1],alpha=0.05)
    axSe.plot(np.mean(radiative_DTM_LAD[Plot_name][:,:-3,1],axis=0)[::-1],heights_rad[3:],'-',c=colour[1],linewidth=2)

    # field inventory
    for i in range(0,n_subplots):
        axSf.fill_betweenx(heights[2:],0,inventory_LAD[Plot_name][i,:][2:],color=colour[2],alpha=0.05)
    axSf.plot(np.mean(inventory_LAD[Plot_name],axis=0)[2:],heights[2:],'-',c=colour[2],linewidth=2)



    axSa.set_ylim(0,80)
    axSb.set_xlim(0,29)
    axSc.set_xlim(xmin=0,xmax=0.7)
    #axSf.set_xlim(xmin=0,xmax=2.0)

    axSb.locator_params(axis='x',nbins=5)
    axSc.locator_params(axis='x',nbins=5)
    axSd.locator_params(axis='x',nbins=5)
    axSe.locator_params(axis='x',nbins=5)
    axSf.locator_params(axis='x',nbins=5)


    yticklabels = axSc.get_yticklabels() + axSd.get_yticklabels() + axSe.get_yticklabels() + axSf.get_yticklabels()

    plt.setp(yticklabels,visible=False)
    plt.subplots_adjust(wspace = 0.1)

    plt.savefig(output_dir+'figS'+str(pp+1)+'_plot_pointclouds.png')


