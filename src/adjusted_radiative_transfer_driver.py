import numpy as np
import sys
from matplotlib import pyplot as plt
import LiDAR_tools as lidar
import auxilliary_functions as aux
import LiDAR_MacHorn_LAD_profiles as LAD1
import LiDAR_radiative_transfer_LAD_profiles as LAD2

las_file = 'Carbon_plot_point_cloud_buffer.las'
subplot_coordinate_file = 'BALI_subplot_coordinates_corrected.csv'
Plots = ['LF','E','Belian','Seraya','B North','B South','DC1','DC2']
N_plots = len(Plots)
leaf_angle_dist = 'spherical'
max_height = 80
max_return = 2
layer_thickness = 1
n_layers = np.ceil(max_height/layer_thickness)
minimum_height = 2.

# store profiles in dictionaries
radiative_spherical_LAD = {}
radiative_spherical_adjusted_LAD = {}
lidar_profiles ={}

subplot_polygons, subplot_labels = aux.load_boundaries(subplot_coordinate_file)

all_pts = lidar.load_lidar_data(las_file)

for pp in range(0,N_plots):
    Plot_name=Plots[pp]
    n_coord_pairs = subplot_polygons[Plot_name].shape[0]*subplot_polygons[Plot_name].shape[1]
    coord_pairs = subplot_polygons[Plot_name].reshape(n_coord_pairs,2)
    bbox_polygon = aux.get_bounding_box(coord_pairs)

    lidar_pts = lidar.filter_lidar_data_by_polygon(all_pts,bbox_polygon)

    n_subplots = subplot_polygons[Plot_name].shape[0]
    subplot_lidar_profiles = np.zeros((n_subplots,n_layers))
    subplot_LAD_profiles_spherical = np.zeros((n_subplots,n_layers+1))
    subplot_LAD_profiles_spherical_adjusted = np.zeros((n_subplots,n_layers+1))
    n_ground_returns = np.zeros(n_subplots)
    subplot_LAI = np.zeros(n_subplots)

    heights_rad = np.arange(0,max_height+1)

    for i in range(0,n_subplots):
        print "Subplot: ", subplot_labels[Plot_name][i]
        sp_pts = lidar.filter_lidar_data_by_polygon(lidar_pts,subplot_polygons[Plot_name][i,:,:])
        #print sp_pts[sp_pts[:,3]==1].shape[0]/20./20.
        heights,subplot_lidar_profiles[i,:],n_ground_returns[i] = LAD1.bin_returns(sp_pts, max_height, layer_thickness)

        u,n,I,U = LAD2.calculate_LAD(sp_pts,heights_rad,max_return,'spherical')
        subplot_LAD_profiles_spherical[i,:]=u.copy()

        #derive correction factor for number of second returns
        #N_1veg = float(np.all((sp_pts[:,3]==1,sp_pts[:,4]==1),axis=0).sum())
        #N_2 = float((sp_pts[:,3]==2).sum())
        #n[:,:,1]*=N_1veg/N_2
        u,n,I,U = LAD2.calculate_LAD_DTM(sp_pts,heights_rad,max_return,'spherical')
        u,n,I,U = LAD2.calculate_LAD(sp_pts,heights_rad,max_return,'spherical',n)
        subplot_LAD_profiles_spherical_adjusted[i,:]=u.copy()

    subplot_LAD_profiles_spherical[np.isnan(subplot_LAD_profiles_spherical)]=0
    subplot_LAD_profiles_spherical_adjusted[np.isnan(subplot_LAD_profiles_spherical_adjusted)]=0

    # remove all profile values below minimum height prior to comparison
    #mask = heights <= minimum_height
    mask = np.max(heights_rad)-heights_rad<=minimum_height
    subplot_LAD_profiles_spherical[:,mask]=0
    subplot_LAD_profiles_spherical_adjusted[:,mask]=0

    # store profiles in dictionaries
    lidar_profiles[Plot_name] = subplot_lidar_profiles
    radiative_spherical_LAD[Plot_name] = subplot_LAD_profiles_spherical[:,:-1]
    radiative_spherical_adjusted_LAD[Plot_name] = subplot_LAD_profiles_spherical_adjusted[:,:-1]

for pp in range(0,N_plots):
    Plot_name=Plots[pp]
    plt.figure(3, facecolor='White',figsize=[10,6.5])
    #plt.title(Plot_name)
    # lidar
    ax31 = plt.subplot2grid((1,3),(0,0))
    ax31.annotate('a', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    for i in range(0,n_subplots):
        ax31.plot(lidar_profiles[Plot_name][i,:]/1000.,heights,'-',c='k',alpha=0.2,linewidth=1)
    ax31.plot(np.mean(lidar_profiles[Plot_name],axis=0)/1000.,heights,'-',c='k',linewidth=1)
    ax31.set_ylim(0,80)
    ax31.set_ylabel('Height / m')
    ax31.set_xlabel('Number of returns (x1000)')
    
    #Radiative Transfer initial
    ax32 = plt.subplot2grid((1,3),(0,1))
    ax32.annotate('b', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    for i in range(0,n_subplots):
        ax32.plot(radiative_spherical_LAD[Plot_name][i,:],np.max(heights)-heights+1,'-',c='r',alpha=0.2,linewidth=1)
    ax32.plot(np.mean(radiative_spherical_LAD[Plot_name],axis=0),np.max(heights)-heights+1,'-',c='k',linewidth=1)
    ax32.set_ylim(0,80)
    ax32.set_xlabel('LAD$_{rad}$ / m$^2$m$^{-1}$')

    #Radiative Transfer adjusted
    ax33 = plt.subplot2grid((1,3),(0,2),sharex=ax32)
    ax33.annotate('c', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    for i in range(0,n_subplots):
        ax33.plot(radiative_spherical_adjusted_LAD[Plot_name][i,:],np.max(heights)-heights+1,'-',c='r',alpha=0.2,linewidth=1)
    ax33.plot(np.mean(radiative_spherical_adjusted_LAD[Plot_name],axis=0),np.max(heights)-heights+1,'-',c='k',linewidth=1)
    ax33.set_ylim(0,80)
    ax33.set_xlabel('LAD$_{rad}$ / m$^2$m$^{-1}$')
    

    ax32.set_xlim(xmax=0.7)
    ax31.set_xlim(xmax=3*np.mean(lidar_profiles[Plot_name],axis=0).max()/1000)

    ax31.locator_params(axis='x',nbins=5)
    ax32.locator_params(axis='x',nbins=5)

    print np.mean(radiative_spherical_adjusted_LAD[Plot_name],axis=0).sum()
    plt.tight_layout()
    plt.savefig(Plot_name+'_LAD_radiative_comparison.png')
    plt.show()
