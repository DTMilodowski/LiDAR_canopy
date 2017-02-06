import numpy as np
import sys
from matplotlib import pyplot as plt
import LiDAR_tools as lidar
import auxilliary_functions as aux
import LiDAR_MacHorn_LAD_profiles as LAD1
import LiDAR_radiative_transfer_LAD_profiles as LAD2

las_file = 'Carbon_plot_point_cloud_buffer.las'
#subplot_coordinate_file = 'BALI_subplot_coordinates_corrected.csv'
coordinate_file = 'BALI_plot_coordinates.csv'
Plots = ['LF','E','Belian','Seraya','B North','B South','DC1','DC2']
N_plots = len(Plots)
leaf_angle_dist = 'spherical'
max_height = 80
max_return = 4
layer_thickness = 1
n_layers = np.ceil(max_height/layer_thickness)
minimum_height = 2.

# store profiles in dictionaries
radiative_spherical_LAD = {}
radiative_spherical_adjusted_LAD = {}
lidar_profiles ={}
lidar_profiles_adjusted = {}

#subplot_polygons, subplot_labels = aux.load_boundaries(subplot_coordinate_file)
subplot_polygons = aux.load_generic_boundaries(coordinate_file)
all_pts = lidar.load_lidar_data(las_file)

for pp in range(0,N_plots):
    Plot_name=Plots[pp]
    bbox_polygon = subplot_polygons[Plot_name] 
    print Plot_name
    lidar_pts = lidar.filter_lidar_data_by_polygon(all_pts,bbox_polygon)

    heights_rad = np.arange(0,max_height+1)

    u,n,I,U = LAD2.calculate_LAD(lidar_pts,heights_rad,max_return,'spherical')
    LAD_profiles_spherical=u.copy()
    lidar_profiles[Plot_name] = np.sum(n.copy(),axis=1)
    #derive correction factor for number of second returns
    #N_1veg = float(np.all((sp_pts[:,3]==1,sp_pts[:,4]==1),axis=0).sum())
    #N_2 = float((sp_pts[:,3]==2).sum())
    #n[:,:,1]*=N_1veg/N_2
    u,n,I,U = LAD2.calculate_LAD_DTM(lidar_pts,heights_rad,max_return,'spherical')
    LAD_profiles_spherical_adjusted=u.copy()
    lidar_profiles_adjusted[Plot_name] = np.sum(n.copy(),axis=1)

    LAD_profiles_spherical[np.isnan(LAD_profiles_spherical)]=0
    LAD_profiles_spherical_adjusted[np.isnan(LAD_profiles_spherical_adjusted)]=0

    # remove all profile values below minimum height prior to comparison
    #mask = heights <= minimum_height
    mask = np.max(heights_rad)-heights_rad<=minimum_height
    LAD_profiles_spherical[mask]=0
    LAD_profiles_spherical_adjusted[mask]=0

    # store profiles in dictionaries
    radiative_spherical_LAD[Plot_name] = LAD_profiles_spherical[:-1]
    radiative_spherical_adjusted_LAD[Plot_name] = LAD_profiles_spherical_adjusted[:-1]

heights = np.arange(1,81)
for pp in range(0,N_plots):
    Plot_name=Plots[pp]
    plt.figure(1, facecolor='White',figsize=[10,6.5])
    ##plt.title(Plot_name)
    # lidar
    ax11 = plt.subplot2grid((1,3),(0,0))
    colour = ['black','blue','red','orange']
    labels = ['$1^{st}$', '$2^{nd}$', '$3^{rd}$', '$4^{th}$']
    ax11.annotate('a', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    for i in range(0,max_return):
        ax11.plot(lidar_profiles[Plot_name][:,i]/1000.,np.max(heights_rad)-heights_rad,'-',c=colour[i],linewidth=1,label=labels[i])
        ax11.plot(lidar_profiles_adjusted[Plot_name][:,i]/1000.,np.max(heights_rad)-heights_rad,'--',c=colour[i],linewidth=1)
    ax11.set_ylim(0,80)
    ax11.set_ylabel('Height / m')
    ax11.set_xlabel('Number of returns (x1000)')
    ax11.legend(loc=1)
    
    #Radiative Transfer initial
    ax12 = plt.subplot2grid((1,3),(0,1))
    ax12.annotate('b', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax12.plot(radiative_spherical_LAD[Plot_name],np.max(heights)-heights+1,'-',c='k',linewidth=1)
    ax12.set_ylim(0,80)
    ax12.set_xlabel('LAD$_{rad}$ / m$^2$m$^{-1}$')

    #Radiative Transfer adjusted
    ax13 = plt.subplot2grid((1,3),(0,2),sharex=ax12)
    ax13.annotate(Plot_name, xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)
    ax13.annotate('c', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax13.plot(radiative_spherical_adjusted_LAD[Plot_name],np.max(heights)-heights+1,'-',c='k',linewidth=1)
    ax13.set_ylim(0,80)
    ax13.set_xlabel('LAD$_{rad}$ / m$^2$m$^{-1}$')
    

    ax12.set_xlim(xmax=0.7)
    
    ax11.locator_params(axis='x',nbins=5)
    ax12.locator_params(axis='x',nbins=5)

    print Plot_name, radiative_spherical_LAD[Plot_name].sum(),  radiative_spherical_adjusted_LAD[Plot_name].sum()
    plt.tight_layout()
    plt.savefig(Plot_name+'_LAD_radiative_comparison_full_plot_inversion.png')
    plt.show()

# A figure illustrating transmittance ratio between successive returns 
plt.figure(2, facecolor='White',figsize=[4,4])
ax21.set_xlabel('return number')
ax21.set_ylabel('transmittance ratio')
for i in range(0,4):
    if i=1:
        ax21.plot(1,1,'o',color='blue')
    else:
        N_veg = float(np.all((all_pts[:,3]==i+1,all_pts[:,4]==1),axis=0).sum())
        N_i = float((all_pts[:,3]==i+1).sum())
        ax21.plot(i+1,N_i/N_veg,'o',color='blue')
ax21.set_ylim(0,1.1)
ax21.set_xlim(0,5)
