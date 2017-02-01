import numpy as np
import sys
from matplotlib import pyplot as plt
import LiDAR_tools as lidar
import auxilliary_functions as aux
import LiDAR_MacHorn_LAD_profiles as LAD1
import LiDAR_radiative_transfer_LAD_profiles as LAD2

sys.path.append('/home/dmilodow/DataStore_DTM/BALI/MetDataProcessing/UtilityTools/')
import statistics_tools as stats

LAI_file = '/home/dmilodow/DataStore_DTM/BALI/BALI_Cplot_data/SAFE_CarbonPlots_LAI_fromHemisphericalPhotos.csv'

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
MacArthurHorn_native = {}
radiative_spherical = {}
radiative_spherical_1st = {}
radiative_planophile = {}
radiative_erectophile = {}
Hemisfer = {}

subplot_polygons, subplot_labels = aux.load_boundaries(subplot_coordinate_file)

field_LAI = aux.load_field_LAI(LAI_file)
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
    subplot_LAD_profiles_planophile = np.zeros((n_subplots,n_layers+1))
    subplot_LAD_profiles_erectophile = np.zeros((n_subplots,n_layers+1))
    subplot_LAD_profiles_native = np.zeros((n_subplots,n_layers))
    subplot_LAD_profiles_spherical_1stOnly=np.zeros((n_subplots,n_layers+1))
    n_ground_returns = np.zeros(n_subplots)
    subplot_LAI = np.zeros(n_subplots)

    heights_rad = np.arange(0,max_height+1)

    for i in range(0,n_subplots):
        print "Subplot: ", subplot_labels[Plot_name][i]
        sp_pts = lidar.filter_lidar_data_by_polygon(lidar_pts,subplot_polygons[Plot_name][i,:,:])
        #print sp_pts[sp_pts[:,3]==1].shape[0]/20./20.
        heights,subplot_lidar_profiles[i,:],n_ground_returns[i] = LAD1.bin_returns(sp_pts, max_height, layer_thickness)
        subplot_LAI[i] = field_LAI['LAI'][np.all((field_LAI['Subplot']==subplot_labels[Plot_name][i],field_LAI['Plot']==Plot_name),axis=0)]

        subplot_LAD_profiles_native[i,:] = LAD1.estimate_LAD_MacArtherHorn(subplot_lidar_profiles[i,:],n_ground_returns[i],layer_thickness,1.)

        u,n,I,U = LAD2.calculate_LAD(sp_pts,heights_rad,1,'spherical')
        subplot_LAD_profiles_spherical_1stOnly[i,:]=u.copy()

        u,n,I,U = LAD2.calculate_LAD(sp_pts,heights_rad,max_return,'spherical')
        subplot_LAD_profiles_spherical[i,:]=u.copy()

        u,n,I,U = LAD2.calculate_LAD(sp_pts,heights_rad,max_return,'planophile',n)
        subplot_LAD_profiles_planophile[i,:]=u.copy()

        u,n,I,U = LAD2.calculate_LAD(sp_pts,heights_rad,max_return,'erectophile',n)
        subplot_LAD_profiles_erectophile[i,:]=u.copy()
        
    subplot_LAD_profiles_spherical[np.isnan(subplot_LAD_profiles_spherical)]=0
    subplot_LAD_profiles_spherical_1stOnly[np.isnan(subplot_LAD_profiles_spherical_1stOnly)]=0
    subplot_LAD_profiles_planophile[np.isnan(subplot_LAD_profiles_planophile)]=0
    subplot_LAD_profiles_erectophile[np.isnan(subplot_LAD_profiles_erectophile)]=0
    kmin = 0.20
    kmax = 5.
    kinc = 0.005
    misfit, ks, best_k_LAD_profiles, best_k = LAD1.minimise_misfit_for_k(kmin,kmax,kinc,subplot_LAI,subplot_lidar_profiles,n_ground_returns,layer_thickness, minimum_height)
    # option here would be to monte-carlo the position estimates, then create misfit vs. ks "envelopes", but runtime for one iteration is ~1 minute, so this is likely to be prohibitive for more than ~500.


    # remove all profile values below minimum height prior to comparison
    mask = heights <= minimum_height
    best_k_LAD_profiles[:,mask]=0
    subplot_LAD_profiles_native[:,mask]=0
    mask = np.max(heights_rad)-heights_rad<=minimum_height
    subplot_LAD_profiles_spherical[:,mask]=0
    subplot_LAD_profiles_spherical_1stOnly[:,mask]=0
    subplot_LAD_profiles_planophile[:,mask]=0
    subplot_LAD_profiles_erectophile[:,mask]=0


    # store profiles in dictionaries
    MacArthurHorn_native[Plot_name] = np.sum(subplot_LAD_profiles_native,axis=1)
    radiative_spherical[Plot_name] = np.sum(subplot_LAD_profiles_spherical[:,:-1],axis=1)
    radiative_spherical_1st[Plot_name] = np.sum(subplot_LAD_profiles_spherical_1stOnly[:,:-1],axis=1)
    radiative_planophile[Plot_name] = np.sum(subplot_LAD_profiles_planophile[:,:-1],axis=1)
    radiative_erectophile[Plot_name] = np.sum(subplot_LAD_profiles_erectophile[:,:-1],axis=1)
    Hemisfer[Plot_name] = subplot_LAI

    # print some stats to screen
    Rsq_1 = stats.calculate_r_squared(subplot_LAI,np.sum(subplot_LAD_profiles_native,axis=1))
    Rsq_2 = stats.calculate_r_squared(subplot_LAI,np.sum(best_k_LAD_profiles,axis=1))
    Rsq_3 = stats.calculate_r_squared(subplot_LAI,np.sum(subplot_LAD_profiles_spherical[:,:-1],axis=1))
    Rsq_4 = stats.calculate_r_squared(subplot_LAI,np.sum(subplot_LAD_profiles_planophile[:,:-1],axis=1))
    Rsq_5 = stats.calculate_r_squared(subplot_LAI,np.sum(subplot_LAD_profiles_erectophile[:,:-1],axis=1))
    print "================================================================="
    print Plot_name
    print "\tHemisfer"
    print "\t\tLAI = ", np.mean(subplot_LAI), "+/-", np.std(subplot_LAI)
    print "Mac-Horn native"
    print "\t\tLAI = ", np.mean(np.sum(subplot_LAD_profiles_native[:,1:],axis=1)), "+/-", np.std(np.sum(subplot_LAD_profiles_native[:,1:],axis=1)), "; R^2 = ", Rsq_1
    print "Mac-Horn calibrated"
    print "\t\tLAI = ", np.mean(np.sum(best_k_LAD_profiles[:,1:],axis=1)), "+/-", np.std(np.sum(best_k_LAD_profiles[:,1:],axis=1)), "; R^2 = ", Rsq_2
    print "Radiative transfer - spherical"
    print "\t\tLAI = ", np.mean(np.sum(subplot_LAD_profiles_spherical[:,:-2],axis=1)), "+/-", np.std(np.sum(subplot_LAD_profiles_spherical[:,:-2],axis=1)), "; R^2 = ", Rsq_3
    print "Radiative transfer - planophile"
    print "\t\tLAI = ", np.mean(np.sum(subplot_LAD_profiles_planophile[:,:-2],axis=1)), "+/-", np.std(np.sum(subplot_LAD_profiles_planophile[:,:-2],axis=1)), "; R^2 = ", Rsq_4
    print "Radiative transfer - erectophile"
    print "\t\tLAI = ", np.mean(np.sum(subplot_LAD_profiles_erectophile[:,:-2],axis=1)), "+/-", np.std(np.sum(subplot_LAD_profiles_erectophile[:,:-2],axis=1)), "; R^2 = ", Rsq_5
    print "================================================================="

    plt.figure(1, facecolor='White',figsize=[12,6])
    ax1 = plt.subplot2grid((2,4),(0,0))
    ax1.annotate('a', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1.annotate('k='+str(best_k[0])+ 'm$^{-1}$', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)
    ax1.plot(ks,misfit,'-k')
    ax1.set_ylim(0,6)
    ax1.set_xlim(0,5)
    ax1.set_ylabel('Average subplot misfit in LAI')
    ax1.set_xlabel('k / m$^{-1}$')

    ax2 = plt.subplot2grid((2,4),(0,1),rowspan=2)
    ax2.annotate('c', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    for i in range(0,n_subplots):
        ax2.plot(subplot_lidar_profiles[i,:]/1000.,heights,'-',c='k',alpha=0.2,linewidth=1)
    ax2.plot(np.mean(subplot_lidar_profiles,axis=0)/1000.,heights,'-',c='k',linewidth=1)
    ax2.set_ylim(0,80)
    ax2.set_ylabel('Height / m')
    ax2.set_xlabel('Number of returns (x1000)')

    ax3 = plt.subplot2grid((2,4),(1,0))
    ax3.annotate('b', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10) 
    ax3.plot([0,10],[0,10],'-',color='0.5')
    ax3.plot(subplot_LAI,np.sum(subplot_LAD_profiles_native,axis=1),'.',c='k',label = 'Mac-Horn$_{native}$')
    ax3.plot(subplot_LAI,np.sum(best_k_LAD_profiles,axis=1),'.',c='b',label = 'Mac-Horn$_{calibrated}$')
    ax3.plot(subplot_LAI,np.sum(subplot_LAD_profiles_spherical[:,:-1],axis=1),'.',c='r',label = 'Radiative Transfer')
    ax3.set_ylim(0,10)
    ax3.set_xlim(0,10)
    ax3.set_xlabel('LAI$_{Hemisfer}$')
    ax3.set_ylabel('LAI$_{LiDAR}$')
    #plt.legend()


    ax4 = plt.subplot2grid((2,4),(0,2),rowspan=2)
    ax4.annotate('d', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    for i in range(0,n_subplots):
        ax4.plot(best_k_LAD_profiles[i,:],heights,'-',c='b',alpha=0.2,linewidth=1)
    ax4.plot(np.mean(best_k_LAD_profiles,axis=0),heights,'-',c='k',linewidth=1)
    ax4.set_ylim(0,80)
    ax4.set_ylabel('Height / m')
    ax4.set_xlabel('LAD$_{MacArthur-Horn}$ / m$^2$m$^{-1}$')

    ax5 = plt.subplot2grid((2,4),(0,3),rowspan=2,sharex=ax4)
    ax5.annotate('e', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax5.annotate(Plot_name, xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10) 
    for i in range(0,n_subplots):
        ax5.plot(subplot_LAD_profiles_spherical[i,:-1],np.max(heights_rad)-heights_rad[:-1],'-',c='r',alpha=0.2,linewidth=1)
    ax5.plot(np.mean(subplot_LAD_profiles_spherical[:,:-1],axis=0),np.max(heights_rad)-heights_rad[:-1],'-',c='k',linewidth=1)
    ax5.set_ylim(0,80)
    ax5.set_ylabel('Height / m')
    ax5.set_xlabel('LAD$_{rad}$ / m$^2$m$^{-1}$')

    ax5.set_xlim(xmax=1.)

    plt.tight_layout()
    plt.savefig(Plot_name+'_LAD_profile.png')


    plt.figure(2, facecolor='White',figsize=[3,3])
    Rsq = stats.calculate_r_squared(np.sum(best_k_LAD_profiles,axis=1),np.sum(subplot_LAD_profiles_spherical[:,:-1],axis=1))
    ax6 = plt.subplot2grid((1,1),(0,0))
    ax6.annotate('R$^2$='+'%.3f' % Rsq, xy=(0.95,0.05), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='bottom', fontsize=10)
    ax6.plot([0,10],[0,10],'-',color='0.5')
    ax6.plot(np.sum(best_k_LAD_profiles,axis=1),np.sum(subplot_LAD_profiles_spherical[:,:-1],axis=1),'.k')
    ax6.set_ylim(0,10)
    ax6.set_xlim(0,10)
    ax6.set_ylabel('LAI$_{rad}$')
    ax6.set_xlabel('LAI$_{MacArthur-Horn}$')

    plt.tight_layout()
    plt.savefig(Plot_name+'_LAI_lidar_comparison.png')

    #plt.show()


plt.figure(4, facecolor='White',figsize=[9,3])
ax7 = plt.subplot2grid((1,3),(0,0))
ax7.set_xlabel('LAI$_{Hemisfer}$')
ax7.set_ylabel('LAI$_{MacArthur-Horn}$')
ax7.annotate('a - MacArthur-Horn', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    
ax7.plot([0,20],[0,20],'--',color='black',alpha=0.3)

ax8 = plt.subplot2grid((1,3),(0,1), sharex=ax7, sharey=ax7)
ax8.annotate('b - Radiative transfer', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax8.set_xlabel('LAI$_{Hemisfer}$')
ax8.set_ylabel('LAI$_{rad}$')

ax8.plot([0,20],[0,20],'--',color='black',alpha=0.3)

ax9 = plt.subplot2grid((1,3),(0,2), sharex=ax7, sharey=ax7)
ax9.annotate('c - LiDAR comparison', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax9.set_xlabel('LAI$_{MacArthur-Horn}$')
ax9.set_ylabel('LAI$_{rad}$')

ax9.plot([0,20],[0,20],'--',color='black',alpha=0.3)


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




plt.figure(3, facecolor='White',figsize=[4,4])
plot_LAI = np.asarray([4.46, 3.76, 3.65, 3.44, 4.40, 4.05])
plot_LAI_error = np.asarray([0.47, 0.48, 0.52, 0.85, 0.53, 0.39])/np.sqrt(25.)

plot_MacHorn = np.asarray([5.94, 5.24, 3.29, 2.80, 5.81, 6.36])
plot_MacHorn_error = np.asarray([1.18, 1.29, 1.50, 1.42, 1.20, 0.80])/np.sqrt(25.)

plot_rad = np.asarray([3.20, 2.89, 2.74, 2.08, 3.13, 3.61])
plot_rad_error = np.asarray([0.67, 0.52, 0.71, 0.47, 1.55, 1.14])/np.sqrt(25.)

plot_sym = ['v', 'o', 'o', 'o', 'v', 'v']
plot_lab = ['Old-growth', 'Logged', 'Logged', 'Logged','Old-growth', 'Old-growth']



ax7 = plt.subplot2grid((1,1),(0,0))

Rsq1 = stats.calculate_r_squared(plot_LAI,plot_MacHorn)
Rsq2 = stats.calculate_r_squared(plot_LAI,plot_rad)
ax7.annotate('$R^2$$_{MacArthur-Horn}$='+'%.3f' % Rsq1 + '\n$R^2$$_{rad}$='+'%.3f' % Rsq2, xy=(0.95,0.05), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='bottom', fontsize=10)

ax7.plot([0,8],[0,8],'-',color='0.5')
for i in range(0,plot_LAI.size):
    ax7.errorbar(plot_LAI[i],plot_MacHorn[i],xerr=plot_LAI_error[i],yerr=plot_MacHorn_error[i],marker=plot_sym[i],color = 'black',label=plot_lab[i])
    ax7.errorbar(plot_LAI[i],plot_rad[i],xerr=plot_LAI_error[i],yerr=plot_rad_error[i],marker=plot_sym[i],color = 'red',label=plot_lab[i])
ax7.set_ylim(2,7)
ax7.set_xlim(2,7)
ax7.set_ylabel('LAI$_{LiDAR}$')
ax7.set_xlabel('LAI$_{MacArthur-Horn}$')

plt.tight_layout()
plt.savefig('Interplot_LAI_comparison.png')
plt.show()
