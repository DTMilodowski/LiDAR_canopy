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

sys.path.append('/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/UtilityTools/')
import statistics_tools as stats

# start by defining input files
las_file = 'Carbon_plot_point_cloud_buffer.las'
subplot_coordinate_file = 'BALI_subplot_coordinates_corrected.csv'
allometry_file = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/Regional/Allometry/Crown_depth_data_SEAsia.csv'
field_file = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/Data/Local/SAFE_carbonplots_FieldMapcensus2016.csv'
LAI_file = '/home/dmilodow/DataStore_DTM/BALI/BALI_Cplot_data/SAFE_CarbonPlots_LAI_fromHemisphericalPhotos.csv'

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
plot_area = 10.**4
subplot_area = 20.*20.
beta = 0.6 # morphological parameter in inventory-based canopy model

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

# load coordinates and lidar points for target areas
subplot_polygons, subplot_labels = aux.load_boundaries(subplot_coordinate_file)
all_lidar_pts = lidar.load_lidar_data(las_file)

# load field data and retrieve allometric relationships
field_data = field.load_crown_survey_data(field_file)
a, b, CF, r_sq, p, H, D = field.retrieve_crown_allometry(allometry_file)
a_ht, b_ht, CF_ht, a_A, b_A, CF_A = field.calculate_allometric_equations_from_survey(field_data)

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
    
    # get some subplot-level information
    n_subplots = subplot_polygons[Plot_name].shape[0]

    # set up some arrays to host the radiative transfer based profiles
    heights_rad = np.arange(0,max_height+1)
    LAD_rad = np.zeros((n_subplots,heights_rad.size,max_return))
    LAD_rad_DTM = np.zeros((n_subplots,heights_rad.size,max_return))

    # set up some arrays to host the MacArthur-Horn profiles
    heights = np.arange(0,max_height)+1
    LAD_MH = np.zeros((n_subplots, heights.size))
    
    # set up array to host the lidar return profiles
    lidar_return_profiles = np.zeros((n_subplots, heights_rad.size, max_return))
    lidar_return_profiles_adj = np.zeros((n_subplots, heights_rad.size, max_return))

    # set up array to host inventory profiles
    field_LAD_profiles = np.zeros((n_subplots,heights.size))

    # set up array to hold the hemisfer LAI estimate
    LAI_hemisfer = np.zeros(n_subplots)

    # loop through subplots, calculating both return profiles and LAD distributions
    for i in range(0,n_subplots):
        print "Subplot: ", subplot_labels[Plot_name][i]
        # filter lidar points into subplot
        sp_pts = lidar.filter_lidar_data_by_polygon(plot_lidar_pts,subplot_polygons[Plot_name][i,:,:])

        # first of all, loop through the return numbers to calculate the radiative LAD profiles
        for rr in range(0,max_return):
            max_k=rr+1
            u,n,I,U = LAD2.calculate_LAD(sp_pts,heights_rad,max_k,'spherical')
            LAD_rad[i,:,rr]=u.copy()
        lidar_return_profiles[i,:,:] = np.sum(n.copy(),axis=1)

        # now repeat but for adjusted profiles, accounting for imperfect penetration of LiDAR pulses into canopy
        for rr in range(0,max_return):
            max_k=rr+1
            u,n,I,U = LAD2.calculate_LAD_DTM(sp_pts,heights_rad,max_k,'spherical')
            LAD_rad_DTM[i,:,rr]=u.copy()
        lidar_return_profiles_adj[i,:,:] = np.sum(n.copy(),axis=1)

        # now get MacArthur-Horn profiles
        heights,first_return_profile,n_ground_returns = LAD1.bin_returns(sp_pts, max_height, layer_thickness)
        LAD_MH[i,:] = LAD1.estimate_LAD_MacArtherHorn(first_return_profile, n_ground_returns, layer_thickness, 1.)

        # now get field inventory estimate
        mask = np.all((field_data['plot']==Plot_name,field_data['subplot']==subplot_labels[Plot_name][i]),axis=0)
        Ht,Area,Depth = field.calculate_crown_dimensions(field_data['DBH_field'][mask],field_data['Height_field'][mask],field_data['CrownArea'][mask], a_ht, b_ht, CF_ht, a_A, b_A, CF_A, a, b, CF)
        field_LAD_profiles[i,:], CanopyV = field.calculate_LAD_profiles_generic(heights, Area, Depth, Ht, beta, subplot_area)
        # now load in the LAI estimates from the hemispherical photographs
        Hemisfer_mask = np.all((field_LAI['Subplot']==subplot_labels[Plot_name][i],field_LAI['Plot']==Plot_name),axis=0)
        LAI_hemisfer[i] = field_LAI['LAI'][Hemisfer_mask]

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
    MacArthurHorn_LAI[Plot_name] = np.sum(LAD_MH,axis=1)
    radiative_LAI[Plot_name] = np.sum(LAD_rad,axis=1)
    radiative_DTM_LAI[Plot_name] = np.sum(LAD_rad_DTM,axis=1)
    inventory_LAI[Plot_name] = np.sum(field_LAD_profiles,axis=1)
    Hemisfer_LAI[Plot_name] = LAI_hemisfer.copy()

#########################################################################################
# Now plot some figures
#----------------------------------------------------------------------------------------
# Figure 1 - Transmittance ratios for the point cloud.
# This figure plots the number of vegetation returns that propagate through further into
# the canopy, defined as the transmittance ratio.
 
# A figure illustrating transmittance ratio between successive returns 
plt.figure(1, facecolor='White',figsize=[4,4])
ax1 = plt.subplot2grid((1,1),(0,0))
ax1.set_xlabel('return number')
ax1.set_ylabel('transmittance ratio')
for i in range(0,4):
    if i==0:
        ax1.plot(1,1,'o',color='blue')
    else:
        N_veg = float(np.all((all_lidar_pts[:,3]==i,all_lidar_pts[:,4]==1),axis=0).sum())
        N_i = float((all_lidar_pts[:,3]==i+1).sum())
        ax1.plot(i+1,N_i/N_veg,'o',color='blue')
ax1.set_ylim(0,1.1)
ax1.set_xlim(0,5)
plt.tight_layout()
plt.savefig(output_dir+'transmittance_ratios.png')
plt.show()

#----------------------------------------------------------------------------------------
# Figure 2 - Allometric relationships used to construct the inventory-based LAD profiles.
plt.figure(2, facecolor='White',figsize=[9,3])
ax2a = plt.subplot2grid((1,3),(0,0))
ax2a.set_ylabel('Crown Depth / m')
ax2a.set_xlabel('Height / m')
ax2a.annotate('a - Height-Crown Depth', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax2a.plot(H,D,'.',color='blue',alpha=0.1)
H_mod = np.linspace(0.,H.max(),1000)
D_mod = CF*a*H_mod**b
ax2a.plot(H_mod,D_mod,'-',color='black')
eq = '$D=%.3fH^{%.3f}$' % (CF*a, b)
ax2a.annotate(eq, xy=(0.95,0.05), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='bottom', fontsize=10)

ax2b = plt.subplot2grid((1,3),(0,1))
ax2b.annotate('b - DBH-Crown Area', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax2b.set_xlabel('DBH / cm')
ax2b.set_ylabel('Crown Area / m$^2$')
ax2b.plot(field_data['DBH_field'],field_data['CrownArea'],'.',color='red',alpha=0.1)
DBH_mod = np.linspace(0.,np.nanmax(field_data['DBH_field']),1000)
CA_mod = CF_A*a_A*DBH_mod**b_A
ax2b.plot(DBH_mod,CA_mod,'-',color='black')
eq = '$A=%.3fDBH^{%.3f}$' % (CF_A*a_A, b_A)
ax2b.annotate(eq, xy=(0.95,0.05), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='bottom', fontsize=10)

ax2c = plt.subplot2grid((1,3),(0,2),sharex=ax2b)
ax2c.annotate('c - DBH-Height', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax2c.set_xlabel('DBH / cm')
ax2c.set_ylabel('Height / m')
ax2c.plot(field_data['DBH_field'],field_data['Height_field'],'.',color='red',alpha=0.1)
Ht_mod = CF_ht*a_ht*DBH_mod**b_ht
ax2c.plot(DBH_mod,Ht_mod,'-',color='black')
eq = '$H=%.3fDBH^{%.3f}$' % (CF_ht*a_ht, b_ht)
ax2c.annotate(eq, xy=(0.95,0.05), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='bottom', fontsize=10)


plt.tight_layout()
plt.savefig(output_dir+'BALI_allometries.png')
plt.show()

#-------------------------------------------------------------------------------------------------
# Figure 3 - Comparison of LiDAR return profiles and the radiative transfer profiles with/without
# compensation for "lost" returns

for pp in range(0,N_plots):
    Plot_name=Plots[pp]
    return_dist = np.mean(lidar_profiles[Plot_name],axis=0)
    return_dist_adj = np.mean(lidar_profiles_adjusted[Plot_name],axis=0)
    LAD_r = np.mean(radiative_LAD[Plot_name],axis=0)
    LAD_r_adj = np.mean(radiative_DTM_LAD[Plot_name],axis=0)


    plt.figure(3, facecolor='White',figsize=[10,6.5])
    ##plt.title(Plot_name)
    # lidar
    ax3a = plt.subplot2grid((1,3),(0,0))
    colour = ['black','blue','red','orange']
    labels = ['$1^{st}$', '$2^{nd}$', '$3^{rd}$', '$4^{th}$']
    ax3a.annotate('a', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    for i in range(0,max_return):
        ax3a.plot(return_dist[:,i]/1000.,np.max(heights_rad)-heights_rad,'-',c=colour[i],linewidth=1,label=labels[i])
        ax3a.plot(return_dist_adj[:,i]/1000.,np.max(heights_rad)-heights_rad,'--',c=colour[i],linewidth=1)
    ax3a.set_ylim(0,80)
    ax3a.set_ylabel('Height / m')
    ax3a.set_xlabel('Number of returns (x1000)')
    ax3a.legend(loc=1)
    ax3a.set_xlim(xmax=1.2*return_dist.max()/1000.)
    
    #Radiative Transfer initial
    ax3b = plt.subplot2grid((1,3),(0,1))
    ax3b.annotate('b', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    for i in range(0,max_return):
        ax3b.plot(LAD_r[:,i],np.max(heights_rad)-heights_rad,'-',c=colour[i],linewidth=1)
    ax3b.set_ylim(0,80)
    ax3b.set_xlabel('LAD$_{rad}$ / m$^2$m$^{-1}$')

    #Radiative Transfer adjusted
    ax3c = plt.subplot2grid((1,3),(0,2),sharex=ax3b)
    ax3c.annotate(Plot_name, xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)
    ax3c.annotate('c', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    for i in range(0,max_return):
        ax3c.plot(LAD_r_adj[:,i],np.max(heights_rad)-heights_rad,'-',c=colour[i],linewidth=1)
    ax3c.set_ylim(0,80)
    ax3c.set_xlabel('LAD$_{rad}$ / m$^2$m$^{-1}$')
    

    ax3c.set_xlim(xmax=0.7)
  
    ax3c.locator_params(axis='x',nbins=5)

    plt.tight_layout()
    plt.savefig(output_dir+Plot_name+'_LAD_radiative_kmax_'+str(max_return)+'.png')
    plt.show()

#----------------------------------------------------------------------------------------------------------------
# Figure 4 - Comparison of subplot and mean LiDAR return profiles and the different methods for estimating LAD
for pp in range(0,N_plots):
    Plot_name=Plots[pp]
    return_dist_adj = np.mean(lidar_profiles_adjusted[Plot_name],axis=0)

    LAD_r = np.mean(radiative_LAD[Plot_name],axis=0)
    LAD_r_adj = np.mean(radiative_DTM_LAD[Plot_name],axis=0)


    plt.figure(4, facecolor='White',figsize=[10,6.5])
    ##plt.title(Plot_name)
    # lidar
    ax4a = plt.subplot2grid((1,4),(0,0))
    colour = ['black','blue','red','orange']
    labels = ['$1^{st}$', '$2^{nd}$', '$3^{rd}$', '$4^{th}$']
    ax4a.annotate('a - ' + Plot_name, xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    for k in range(0,max_return):
        ax4a.plot(return_dist_adj[:,k]/1000.,np.max(heights_rad)-heights_rad,'-',c=colour[i],linewidth=1)
    ax4a.set_ylim(0,80)
    ax4a.set_ylabel('Height / m')
    ax4a.set_xlabel('Number of returns (x1000)')
    ax4a.legend(loc=1)
    ax4a.set_xlim(xmax=1.2*return_dist.max()/1000.)
    
    # MacArthur-Horn method
    ax4b = plt.subplot2grid((1,4),(0,1))
    ax4b.annotate('b', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax4b.annotate('MacArthur-Horn', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)
    for i in range(0,n_subplots):
        ax4b.plot(MacArthurHorn_LAD[Plot_name][i,:],heights,'-',c='blue',linewidth=0.5,alpha=0.5)
    ax4b.plot(np.mean(MacArthurHorn_LAD[Plot_name],axis=0),heights,'-',c='blue',linewidth=1)
    ax4b.set_ylim(0,80)
    ax4b.set_xlabel('LAD$_{MacArthur-Horn}$ / m$^2$m$^{-1}$')

    #Radiative Transfer adjusted
    ax4c = plt.subplot2grid((1,4),(0,2),sharex=ax3b)
    ax4c.annotate('c', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax4c.annotate('Radiative transfer', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)
    for i in range(0,max_return):
        ax4c.plot(radiative_DTM_LAD[Plot_name][i,:,-1],np.max(heights_rad)-heights_rad,'-',c='red',linewidth=0.5,alpha=0.5)
    ax4c.plot(np.mean(radiative_DTM_LAD[Plot_name][:,:,-1],axis=0),np.max(heights_rad)-heights_rad,'-',c='red',linewidth=1)
    ax4c.set_ylim(0,80)
    ax4c.set_xlabel('LAD$_{rad}$ / m$^2$m$^{-1}$')

    #Field Inventory
    ax4d = plt.subplot2grid((1,4),(0,3))
    ax4d.annotate('Field inventory', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)
    ax4d.annotate('d', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    for i in range(0,max_return):
        ax4d.plot(inventory_LAD[Plot_name][i,:],heights,'-',c='green',alpha=0.5,linewidth=0.5)
    ax4d.plot(np.mean(inventory_LAD_r_adj[Plot_name],axis=0),heights,'-',c='green',linewidth=1)
    ax4d.set_ylim(0,80)
    ax4d.set_xlabel('Canopy Volume / m$^3$m$^{-2}$')
    

    ax4c.set_xlim(xmax=0.7)
  
    ax4c.locator_params(axis='x',nbins=5)
    ax4d.locator_params(axis='x',nbins=5)

    plt.tight_layout()
    plt.savefig(output_dir+Plot_name+'_LAD_comparison_kmax_'+str(max_return)+'.png')
    plt.show()



#----------------------------------------------------------------------------------------------------------------
# Figure 5 - Comparing LAI from hemiphotos against integrated LAD from the different methods
# Plot up the subplot-level estimates, in addition to the 
plt.figure(4, facecolor='White',figsize=[9,3])
ax7 = plt.subplot2grid((1,3),(0,0))
ax7.set_xlabel('LAI$_{Hemisfer}$')
ax7.set_ylabel('LAI$_{MacArthur-Horn}$')
ax7.annotate('a - MacArthur-Horn', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax7.plot([0,20],[0,20],'--',color='black',alpha=0.3)
for pp in range(0,N_plots):
    ax7.plot(Hemisfer_LAI[Plots[i]],MacArthurHorn_LAI[Plots[i]],'.',color='blue',alpha=0.5)

for pp in range(0,N_plots):
    x_err=np.std(Hemisfer_LAI[Plots[i]])/np.sqrt(N_plots)
    y_err=np.std(MacArthurHorn_LAI[Plots[i]])/np.sqrt(N_plots)
    ax7.errorbar(np.mean(Hemisfer_LAI[Plots[i]]),np.mean(MacArthurHorn_LAI[Plots[i]]),xerr=x_err,yerr=y_err,'o',color='black')


ax8 = plt.subplot2grid((1,3),(0,1), sharex=ax7, sharey=ax7)
ax8.annotate('b - Radiative transfer', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax8.set_xlabel('LAI$_{Hemisfer}$')
ax8.set_ylabel('LAI$_{rad}$')
ax8.plot([0,20],[0,20],'--',color='black',alpha=0.3)
for pp in range(0,N_plots):
    ax7.plot(Hemisfer_LAI[Plots[i]],radiative_LAI[Plots[i]],'.',color='0.5',alpha=0.5)
    ax7.plot(Hemisfer_LAI[Plots[i]],radiative_DTM_LAI[Plots[i]],'.',color='blue',alpha=0.5)

for pp in range(0,N_plots):
    x_err=np.std(Hemisfer_LAI[Plots[i]])/np.sqrt(N_plots)
    y_err1=np.std(radiative_LAI[Plots[i]])/np.sqrt(N_plots)
    y_err2=np.std(radiative_DTM_LAI[Plots[i]])/np.sqrt(N_plots)
    ax7.errorbar(np.mean(Hemisfer_LAI[Plots[i]]),np.mean(radiative_LAI[Plots[i]]),xerr=x_err,yerr=y_err1,'o',color='black',mfc='white')
    ax7.errorbar(np.mean(Hemisfer_LAI[Plots[i]]),np.mean(radiative_DTM_LAI[Plots[i]]),xerr=x_err,yerr=y_err2,'o',color='black')


ax9 = plt.subplot2grid((1,3),(0,2), sharex=ax7)
ax9.annotate('c - Field inventory', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax9.set_xlabel('LAI$_{Hemisfer}$')
ax9.set_ylabel('Canopy Volume / $m^3m^{-2}$')

#ax9.plot([0,20],[0,20],'--',color='black',alpha=0.3)


for i in range(0,len(Plots)):
    ax7.plot(Hemisfer[Plots[i]],MacArthurHorn_native[Plots[i]],'.',color='blue',alpha=0.5)
    ax8.plot(Hemisfer[Plots[i]],radiative_spherical[Plots[i]],'.',color='red', alpha=0.5)
    ax8.plot(Hemisfer[Plots[i]],radiative_spherical_1st[Plots[i]],'.',color='black', alpha=0.5)
    ax9.plot(MacArthurHorn_native[Plots[i]],radiative_spherical[Plots[i]],'.',color='red', alpha=0.5)
    ax9.plot(MacArthurHorn_native[Plots[i]],radiative_spherical_1st[Plots[i]],'.',color='black', alpha=0.5)
    

ax7.set_xlim((0,20))
ax7.set_ylim((0,20))
plt.tight_layout()
plt.savefig('GEM_subplot_compilation.png')
plt.show()


