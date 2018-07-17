import numpy as np
from matplotlib import pyplot as plt
import LiDAR_io as io
import LiDAR_tools as lidar
import auxilliary_functions as aux
import LiDAR_MacHorn_LAD_profiles as LAD1
import LiDAR_radiative_transfer_LAD_profiles as LAD2
import inventory_based_LAD_profiles as field
import least_squares_fitting as lstsq
from scipy import stats
from osgeo import gdal
import numpy.ma as ma
from matplotlib import rcParams
from datetime import datetime

import fiona
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.ticker as plticker
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import colorbar
from mpl_toolkits.axes_grid1 import make_axes_locatable

from matplotlib.collections import LineCollection

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import sys

# Get perceptually uniform colourmaps
sys.path.append('/home/dmilodow/DataStore_DTM/FOREST2020/EOdata/EO_data_processing/src/')
import geospatial_utility_tools as geo
sys.path.append('/home/dmilodow/DataStore_DTM/FOREST2020/EOdata/EO_data_processing/src/plot_EO_data/colormap/')
import colormaps as cmaps
plt.register_cmap(name='plasma', cmap=cmaps.plasma)
plt.register_cmap(name='viridis', cmap=cmaps.viridis)

# Set up some basiic parameters for the plots
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = 12
rcParams['legend.numpoints'] = 1
axis_size = rcParams['font.size']+2
colour = ['#46E900','#1A2BCE','#E0007F']
n_subplots=25

# generic function to plot an arbitrary line spanning the axes (usually a 1:1 line)
def abline(ax,slope,intercept,colour='0.5',linestyle = '--'):
    x = np.array(ax.get_xlim())
    y = intercept + slope * x
    ax.plot(x, y, color = colour, linestyle = linestyle)
    return 0
def plot_one2one(ax,colour='0.5',linestyle = '--'):
    test = abline(ax,1.,0.,colour=colour,linestyle = linestyle)
    return 0

# function to plot a line with line segments coloured according to third variable.
def plot_colourline(ax,x,y,t,linewidth=10,cmap='plasma'):
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap='plasma',norm=plt.Normalize(t.min(), t.max()))
    lc.set_array(t)
    lc.set_linewidth(linewidth)
    lc.set_zorder(100)
    ax.add_collection(lc)
    return 0

# plot point clouds with canopy profiles
# six rows for six plots (2x old growth, 2x moderately logged, 2x heavily logged)
def plot_point_clouds_and_profiles(figure_name,figure_number, gps_pts_file,plot_point_cloud,heights,heights_rad, lidar_profiles,MacArthurHorn_LAD,MacArthurHorn_LAD_mean,radiative_LAD,radiative_LAD_mean,radiative_DTM_LAD,radiative_DTM_LAD_mean,inventory_LAD):

    max_return=3
    n_subplots = 25
    colour = ['#46E900','#1A2BCE','#E0007F']
    rgb = [[70,233,0],[26,43,206],[224,0,127]]
    labels = ['$1^{st}$', '$2^{nd}$', '$3^{rd}$', '$4^{th}$']
    
    # first up, going to need to find affine transformation to rotate the point cloud for
    # easy plotting
    
    # The GPS coordinates for the plot locations can be used to find this transformation matrix
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
    fig1_plots = ['Belian', 'Seraya','LF','E', 'B North', 'B South']
    Plots =      ['Belian', 'Seraya','LF','E', 'B North', 'B South']
    N_plots = len(fig1_plots)

    plot_point_cloud_display = {}
    
    for pp in range(0, N_plots):
        # first get points for a given plot and build matrices - note that I've reversed xy and xy_prime in this case to reverse the rotation-translation
        mask = plot_coords['plot']==fig1_plots[pp]
        x = plot_coords['x_prime'][mask]
        y = plot_coords['y_prime'][mask]
        x_prime = plot_coords['x'][mask]
        y_prime = plot_coords['y'][mask]
        affine[Plots[pp]]=lstsq.least_squares_affine_matrix(x,y,x_prime,y_prime)
        
        Xi = np.asarray([plot_point_cloud[Plots[pp]][:,0],plot_point_cloud[Plots[pp]][:,1],np.ones(plot_point_cloud[Plots[pp]].shape[0])])
        Xi_prime = np.dot(affine[Plots[pp]],Xi)
        
        plot_point_cloud_display[Plots[pp]] = plot_point_cloud[Plots[pp]].copy()
        plot_point_cloud_display[Plots[pp]][:,0]=Xi_prime[0]
        plot_point_cloud_display[Plots[pp]][:,1]=Xi_prime[1]

    plt.figure(figure_number, facecolor='White',figsize=[9,12])

    # Belian
    ax1a = plt.subplot2grid((6,7),(0,0),colspan=2)
    ax1a.annotate('a - Old growth, MLA01', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1a.set_ylabel('Height / m',fontsize=axis_size)
    plt.gca().set_aspect('equal', adjustable='box-forced')
    
    # Seraya
    ax1b = plt.subplot2grid((6,7),(1,0),sharey=ax1a,sharex=ax1a,colspan=2)
    ax1b.annotate('b - Old growth, MLA02', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1b.set_ylabel('Height / m',fontsize=axis_size)
    plt.gca().set_aspect('equal', adjustable='box-forced')
    
    # LF
    ax1c = plt.subplot2grid((6,7),(2,0),sharey=ax1a,sharex=ax1a,colspan=2)
    ax1c.annotate('c - Moderately logged, SAF04', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1c.set_ylabel('Height / m',fontsize=axis_size)
    plt.gca().set_aspect('equal', adjustable='box-forced')
    
    # E
    ax1d = plt.subplot2grid((6,7),(3,0),sharey=ax1a,sharex=ax1a,colspan=2)
    ax1d.annotate('d - Moderately logged, SAF05', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1d.set_ylabel('Height / m',fontsize=axis_size)
    plt.gca().set_aspect('equal', adjustable='box-forced')
    
    # B North
    ax1e = plt.subplot2grid((6,7),(4,0),sharey=ax1a,sharex=ax1a,colspan=2)
    ax1e.annotate('e - Heavily logged, SAF02', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1e.set_ylabel('Height / m',fontsize=axis_size)
    plt.gca().set_aspect('equal', adjustable='box-forced')
    
    # B South
    ax1f = plt.subplot2grid((6,7),(5,0),sharey=ax1a,sharex=ax1a,colspan=2)
    ax1f.annotate('f - Heavily logged, SAF01', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1f.set_ylabel('Height / m',fontsize=axis_size)
    ax1f.set_xlabel('Horizontal distance / m',fontsize=axis_size)
    plt.gca().set_aspect('equal', adjustable='box-forced')
    
    axes = [ax1a, ax1b, ax1c, ax1d, ax1e, ax1f]
    for pp in range(0,6):
        plot_lidar_pts = plot_point_cloud_display[fig1_plots[pp]]
        for k in range(0,max_return):
            
            mask = np.all((plot_lidar_pts[:,0]>=0,plot_lidar_pts[:,0]<=100,plot_lidar_pts[:,1]>=0,plot_lidar_pts[:,0]<=100,plot_lidar_pts[:,3]==k+1),axis=0)
            points_x = 100-plot_lidar_pts[mask][:,0]
            points_z = plot_lidar_pts[mask][:,2]
            points_y = plot_lidar_pts[mask][:,1]
            
            alpha_max = 0.1
            colours = np.zeros((points_x.size,4))
            colours[:,0]=rgb[k][0]/255.
            colours[:,1]=rgb[k][1]/255.
            colours[:,2]=rgb[k][2]/255.
            
            colours[:,3]=alpha_max*(1-points_x/(points_x.max()+1))
            axes[pp].scatter(points_y,points_z,marker='o',c=colours,edgecolors='none',s=1)
            axes[pp].scatter(0,0,marker='o',c=colours[0,0:3],edgecolors='none',s=1,label=labels[k])

    #x1a.legend(loc=1,fontsize=axis_size)

    #---------------------------------------------------------
    # NOW PLOT PROFILES
    # Belian
    ax2a = plt.subplot2grid((6,7),(0,2),sharey=ax1a)
    ax2a.annotate('LiDAR returns', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=9)
    # - MacHorn
    ax3a = plt.subplot2grid((6,7),(0,3),sharey=ax1a)
    ax3a.annotate('MacArthur-Horn', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=9)
    # - Detto
    ax4a = plt.subplot2grid((6,7),(0,4),sharey=ax1a,sharex=ax3a)
    ax4a.annotate('multi. return\n(Detto)', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=9)
    # - Corrected rad trans
    ax5a = plt.subplot2grid((6,7),(0,5),sharey=ax1a,sharex=ax3a)
    ax5a.annotate('multi. return\n(new)', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=9)
    # - Inventory
    ax6a = plt.subplot2grid((6,7),(0,6),sharey=ax1a,sharex=ax3a)
    ax6a.annotate('crown volume', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=9)

    # Seraya
    ax2b = plt.subplot2grid((6,7),(1,2),sharey=ax1a,sharex=ax2a)
    ax3b = plt.subplot2grid((6,7),(1,3),sharey=ax1a,sharex=ax3a)
    ax4b = plt.subplot2grid((6,7),(1,4),sharey=ax1a,sharex=ax3a)
    ax5b = plt.subplot2grid((6,7),(1,5),sharey=ax1a,sharex=ax3a)
    ax6b = plt.subplot2grid((6,7),(1,6),sharey=ax1a,sharex=ax3a)

    # LF
    ax2c = plt.subplot2grid((6,7),(2,2),sharey=ax1a,sharex=ax2a)
    ax3c = plt.subplot2grid((6,7),(2,3),sharey=ax1a,sharex=ax3a)
    ax4c = plt.subplot2grid((6,7),(2,4),sharey=ax1a,sharex=ax3a)
    ax5c = plt.subplot2grid((6,7),(2,5),sharey=ax1a,sharex=ax3a)
    ax6c = plt.subplot2grid((6,7),(2,6),sharey=ax1a,sharex=ax3a)

    # E
    ax2d = plt.subplot2grid((6,7),(3,2),sharey=ax1a,sharex=ax2a)
    ax3d = plt.subplot2grid((6,7),(3,3),sharey=ax1a,sharex=ax3a)
    ax4d = plt.subplot2grid((6,7),(3,4),sharey=ax1a,sharex=ax3a)
    ax5d = plt.subplot2grid((6,7),(3,5),sharey=ax1a,sharex=ax3a)
    ax6d = plt.subplot2grid((6,7),(3,6),sharey=ax1a,sharex=ax3a)

    # B North
    ax2e = plt.subplot2grid((6,7),(4,2),sharey=ax1a,sharex=ax2a)
    ax3e = plt.subplot2grid((6,7),(4,3),sharey=ax1a,sharex=ax3a)
    ax4e = plt.subplot2grid((6,7),(4,4),sharey=ax1a,sharex=ax3a)
    ax5e = plt.subplot2grid((6,7),(4,5),sharey=ax1a,sharex=ax3a)
    ax6e = plt.subplot2grid((6,7),(4,6),sharey=ax1a,sharex=ax3a)

    # B South
    ax2f = plt.subplot2grid((6,7),(5,2), sharex = ax1a, sharey = ax2a)
    ax2f.set_xlabel('Number of returns\n(x1000)',fontsize=axis_size,horizontalalignment='center')
    # - MacHorn
    ax3f = plt.subplot2grid((6,7),(5,3),sharey=ax1a, sharex=ax3a)
    ax3f.set_xlabel('PAD\n(m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size,horizontalalignment='center')
    # - Detto
    ax4f = plt.subplot2grid((6,7),(5,4),sharey=ax1a,sharex=ax3a)
    ax4f.set_xlabel('PAD\n(m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size,horizontalalignment='center')
    # - Corrected rad trans
    ax5f = plt.subplot2grid((6,7),(5,5),sharey=ax1a,sharex=ax3a)
    ax5f.set_xlabel('PAD\n(m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size,horizontalalignment='center')
    # - Inventory
    ax6f = plt.subplot2grid((6,7),(5,6),sharey=ax1a,sharex=ax3a)
    ax6f.set_xlabel('Crown Volume\n(m$^3$m$^{-2}$m$^{-1}$)',fontsize=axis_size,horizontalalignment='center')

    axes1 = [ax2a,  ax2b, ax2c, ax2d, ax2e, ax2f] 
    axes2 = [ax3a,  ax3b, ax3c, ax3d, ax3e, ax3f] 
    axes3 =  [ax4a,  ax4b, ax4c, ax4d, ax4e, ax4f] 
    axes4 = [ax5a,  ax5b, ax5c, ax5d, ax5e, ax5f] 
    axes5 = [ax6a,  ax6b, ax6c, ax6d, ax6e, ax6f]

    yticklabels=[]
    xticklabels=[]
    xticklabels.append(ax1a.get_xticklabels() + ax1b.get_xticklabels() + ax1c.get_xticklabels() + ax1d.get_xticklabels() + ax1e.get_xticklabels())
    
    for pp in range(0,6):
        Plot_name = fig1_plots[pp]
        
        # plot lidar profile
        return_dist     = np.sum(lidar_profiles[Plot_name],axis=0)
        for k in range(0,max_return):
            axes1[pp].plot(return_dist[:,k]/1000.,np.max(heights_rad)-heights_rad,'-',c=colour[k],linewidth=1)
            
            # plot macarthur horn profile
            for i in range(0,n_subplots):
                axes2[pp].fill_betweenx(heights[2:],0,MacArthurHorn_LAD[Plot_name][i,2:],color=colour[0],alpha=0.01)
            axes2[pp].plot(MacArthurHorn_LAD_mean[Plot_name][2:],heights[2:],'-',c=colour[0],linewidth=2)
            # plot detto profile
            for i in range(0,n_subplots):
                axes3[pp].fill_betweenx(heights_rad[3:],0,radiative_LAD[Plot_name][i,:-3,-1][::-1],color=colour[1],alpha=0.01)
            axes3[pp].plot(radiative_LAD_mean[Plot_name][:-3,1][::-1],heights_rad[3:],'-',c=colour[1],linewidth=2)
                    
            # plot corrective radiative transfer profile
            for i in range(0,n_subplots):
                axes4[pp].fill_betweenx(heights_rad[3:],0,radiative_DTM_LAD[Plot_name][i,:-3,-1][::-1],color=colour[1],alpha=0.01)
            axes4[pp].plot(radiative_DTM_LAD_mean[Plot_name][:-3,1][::-1],heights_rad[3:],'-',c=colour[1],linewidth=2)

            # field inventory
            #for i in range(0,n_subplots):
                #axes5[pp].fill_betweenx(heights[2:],0,inventory_LAD[Plot_name][i,2:],color=colour[2],alpha=0.05)
            #axes5[pp].plot(np.mean(inventory_LAD[Plot_name],axis=0)[2:],heights[2:],'-',c=colour[2],linewidth=2)
            axes5[pp].plot(inventory_LAD[Plot_name][2:],heights[2:],'-',c=colour[2],linewidth=2)

        yticklabels.append(axes1[pp].get_yticklabels())
        yticklabels.append(axes2[pp].get_yticklabels())
        yticklabels.append(axes3[pp].get_yticklabels())
        yticklabels.append(axes4[pp].get_yticklabels())
        yticklabels.append(axes5[pp].get_yticklabels())
            
        if pp < 5:
            xticklabels.append(axes1[pp].get_xticklabels())
            xticklabels.append(axes2[pp].get_xticklabels())
            xticklabels.append(axes3[pp].get_xticklabels())
            xticklabels.append(axes4[pp].get_xticklabels())
            xticklabels.append(axes5[pp].get_xticklabels())

    ax1a.set_xlim(0,100)
    ax1a.set_ylim(0,80)
    ax2a.set_xlim(0,29)
    ax3a.set_xlim(xmin=0,xmax=0.7)

    ax2f.locator_params(axis='x',nbins=5)
    ax3f.locator_params(axis='x',nbins=5)
    ax4f.locator_params(axis='x',nbins=5)
    ax5f.locator_params(axis='x',nbins=5)
    ax6f.locator_params(axis='x',nbins=5)

    plt.setp(yticklabels,visible=False)
    plt.setp(xticklabels,visible=False)
    plt.subplots_adjust(hspace=0.2, wspace = 0.1)

    plt.tight_layout()
    plt.savefig(figure_name)
    plt.show()
    return 0
    
#=======================================================================================
# Compare LiDAR approaches
def compare_LiDAR_PAI(figure_name,figure_number,MacArthurHorn_LAD,MacArthurHorn_LAD_mean,radiative_LAD,radiative_LAD_mean,radiative_DTM_LAD,radiative_DTM_LAD_mean,layer_thickness=1):

    Plots = MacArthurHorn_LAD.keys()
    N_plots = len(Plots)
    N_subplots = np.shape(MacArthurHorn_LAD[Plots[0]])[0]
    mh_1ha = np.zeros(N_plots)
    mh_20m = np.zeros((N_plots,N_subplots))
    rad2_1ha = np.zeros(N_plots)
    rad2_20m = np.zeros((N_plots,N_subplots))
    rad3_1ha = np.zeros(N_plots)
    rad3_20m = np.zeros((N_plots,N_subplots))
    rad2_DTM_1ha = np.zeros(N_plots)
    rad2_DTM_20m = np.zeros((N_plots,N_subplots))
    rad3_DTM_1ha = np.zeros(N_plots)
    rad3_DTM_20m = np.zeros((N_plots,N_subplots))

    for pp in range(0,N_plots):
        print Plots[pp]
        mh_20m[pp,:]=np.nansum(MacArthurHorn_LAD[Plots[pp]],axis=1)*layer_thickness
        rad2_20m[pp,:]=np.nansum(radiative_LAD[Plots[pp]][:,:,1],axis=1)*layer_thickness
        rad3_20m[pp,:]=np.nansum(radiative_LAD[Plots[pp]][:,:,2],axis=1)*layer_thickness
        rad2_DTM_20m[pp,:]=np.nansum(radiative_DTM_LAD[Plots[pp]][:,:,1],axis=1)*layer_thickness
        rad3_DTM_20m[pp,:]=np.nansum(radiative_DTM_LAD[Plots[pp]][:,:,2],axis=1)*layer_thickness
        print rad2_DTM_20m[pp].astype('int')
        print mh_20m[pp]
        
        mh_1ha[pp] = np.nansum(MacArthurHorn_LAD_mean[Plots[pp]])*layer_thickness
        rad2_1ha[pp] = np.nansum(radiative_LAD_mean[Plots[pp]][:,1])*layer_thickness
        rad3_1ha[pp] = np.nansum(radiative_LAD_mean[Plots[pp]][:,2])*layer_thickness
        rad2_DTM_1ha[pp] = np.nansum(radiative_DTM_LAD_mean[Plots[pp]][:,1])*layer_thickness
        rad3_DTM_1ha[pp] = np.nansum(radiative_DTM_LAD_mean[Plots[pp]][:,2])*layer_thickness
        print rad2_DTM_1ha[pp]  
        print mh_1ha[pp]      
        
    # annotate with stats
    r_sq_a = [aux.get_rsquared_annotation(mh_1ha,rad2_1ha), aux.get_rsquared_annotation(mh_1ha,rad3_1ha)]
    r_sq_b = [aux.get_rsquared_annotation(mh_1ha,rad2_DTM_1ha), aux.get_rsquared_annotation(mh_1ha,rad3_DTM_1ha)]

    plt.figure(figure_number, facecolor='White',figsize=[7,4])
    ax1 = plt.subplot2grid((1,2),(0,0))
    ax1.set_xlabel('PAI$_{MacArthur-Horn}$',fontsize=axis_size)
    ax1.set_ylabel('PAI$_{rad}$',fontsize=axis_size)
    ax1.annotate('a - radiative transfer (Detto)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1.annotate('$k_{max}=2$; ' + r_sq_a[0] + '\n$k_{max}=3$; ' + r_sq_a[1], xy=(0.95,0.90), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)
    
    ax1.plot([0,20],[0,20],'--',color='black',alpha=0.3)
    ax1.plot(mh_20m,rad2_20m,'.',color=colour[1],alpha=0.5)
    ax1.plot(mh_20m,rad3_20m,'.',color=colour[2],alpha=0.5)

    x_err=np.std(mh_20m,axis=1)/np.sqrt(n_subplots)
    y_err=np.std(rad2_20m,axis=1)/np.sqrt(n_subplots)
    ax1.errorbar(mh_1ha,rad2_1ha,xerr=x_err,yerr=y_err,marker='o',linestyle='None',color=colour[1])
    y_err=np.std(rad3_20m,axis=1)/np.sqrt(n_subplots)
    ax1.errorbar(mh_1ha,rad3_1ha,xerr=x_err,yerr=y_err,marker='o',linestyle='None',color=colour[2])
    
    ax2 = plt.subplot2grid((1,2),(0,1), sharex = ax1, sharey = ax1)
    ax2.set_xlabel('PAI$_{MacArthur-Horn}$',fontsize=axis_size)
    ax2.set_ylabel('PAI$_{rad}$',fontsize=axis_size)
    ax2.annotate('b - radiative transfer (modified)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax2.annotate('$k_{max}=2$; ' + r_sq_b[0] + '\n$k_{max}=3$; ' + r_sq_b[1], xy=(0.95,0.90), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)

    ax2.plot([0,20],[0,20],'--',color='black',alpha=0.3)
    ax2.plot(mh_20m,rad2_DTM_20m,'.',color=colour[1],alpha=0.5)
    ax2.plot(mh_20m,rad3_DTM_20m,'.',color=colour[2],alpha=0.5)

    y_err=np.std(rad2_DTM_20m,axis=1)/np.sqrt(n_subplots)
    ax2.errorbar(mh_1ha,rad2_DTM_1ha,xerr=x_err,yerr=y_err,marker='o',linestyle='None',color=colour[1],label = '$k_{max}=2$')
    y_err=np.std(rad3_DTM_20m,axis=1)/np.sqrt(n_subplots)
    ax2.errorbar(mh_1ha,rad3_DTM_1ha,xerr=x_err,yerr=y_err,marker='o',linestyle='None',color=colour[2],label = '$k_{max}=3$')

    # configs
    ax2.legend(loc=4)
    ax2.legend(loc=4)
    ax2.set_xlim((0,20))
    ax2.set_ylim((0,20))

    plt.tight_layout()
    plt.savefig(figure_name)
    plt.show()
    return 0


#===============================================================================
# LAI vs. canopy volume
def plot_LAI_vs_inventory(figure_name,figure_number,MacArthurHorn_LAD,MacArthurHorn_LAD_mean,radiative_LAD,radiative_LAD_mean,radiative_DTM_LAD,radiative_DTM_LAD_mean,Inventory_LAD,Inventory_LAI,layer_thickness=1):
    Plots = MacArthurHorn_LAD.keys()
    N_plots = len(Plots)
    N_subplots = np.shape(MacArthurHorn_LAD[Plots[0]])[0]
    vol = np.zeros(N_plots)
    mh = np.zeros(N_plots)
    rad2 = np.zeros(N_plots)
    rad3 = np.zeros(N_plots)
    radDTM2 = np.zeros(N_plots)
    radDTM3 = np.zeros(N_plots)

    for pp in range(0,N_plots):
        print Plots[pp]
        # for 1 ha averages, want to avoid nodata
        # these should be identifiable based on penetration limit
        vol[pp] = np.sum(Inventory_LAD[Plots[pp]])
        mh[pp] = np.nansum(MacArthurHorn_LAD_mean[Plots[pp]])*layer_thickness
        rad2[pp] = np.nansum(radiative_LAD_mean[Plots[pp]][:,1])*layer_thickness
        rad3[pp] = np.nansum(radiative_LAD_mean[Plots[pp]][:,2])*layer_thickness
        radDTM2[pp] = np.nansum(radiative_DTM_LAD_mean[Plots[pp]][:,1])*layer_thickness
        radDTM3[pp] = np.nansum(radiative_DTM_LAD_mean[Plots[pp]][:,2])*layer_thickness
    
    # annotate with stats
    r_sq_a =   aux.get_rsquared_annotation(vol,mh)
    r_sq_b = [aux.get_rsquared_annotation(vol,rad2),aux.get_rsquared_annotation(vol,rad3)]
    r_sq_c = [aux.get_rsquared_annotation(vol,radDTM2),aux.get_rsquared_annotation(vol,radDTM3)]

    plt.figure(figure_number, facecolor='White',figsize=[12,4])
    ax1 = plt.subplot2grid((1,3),(0,0))
    ax1.set_xlabel('crown volume / m$^3$m$^{-2}$',fontsize=axis_size)
    ax1.set_ylabel('PAI',fontsize=axis_size)
    ax1.annotate('a - MacArthur-Horn', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=axis_size)
    ax1.annotate(r_sq_a, xy=(0.95,0.05), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='bottom', fontsize=axis_size)
    #for i in range(0,N_plots):
        #ax1.plot(inventory_LAI[Plots[i]],MacArthurHorn_LAI[Plots[i]],'.',color=colour[0],alpha=0.5)
    # plot linear regression
    m, c, r2, p, x_i, y_i, CI_u, CI_l, PI_u, PI_l = linear_regression(vol,mh,conf=0.95)
    ax1.fill_between(x_i,CI_l,CI_u,color='0.75')
    ax1.plot(x_i,y_i,'-',color='black')

    for i in range(0,N_plots):
        x_err=np.nan#np.std(inventory_LAI[Plots[i]])/np.sqrt(n_subplots)
        y_err=np.nan#np.std(MacArthurHorn_LAI[Plots[i]])/np.sqrt(n_subplots)
        #ax1.errorbar(np.mean(Inventory_LAI[Plots[i]]),np.mean(MacArthurHorn_LAI[Plots[i]]),xerr=x_err,yerr=y_err,marker='o',color='black')
        ax1.errorbar(np.mean(Inventory_LAI[Plots[i]]),np.nansum(MacArthurHorn_LAD_mean[Plots[i]])*layer_thickness,xerr=x_err,yerr=y_err,marker='o',color='black')


    ax2 = plt.subplot2grid((1,3),(0,1), sharex=ax1)#, sharey=ax1)
    ax2.annotate('b - radiative transfer (Detto)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=axis_size)
    ax2.annotate('$k_{max}=2$; ' + r_sq_b[0] + '\n$k_{max}=3$; ' + r_sq_b[1], xy=(0.95,0.05), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='bottom', fontsize=axis_size)
    ax2.set_xlabel('crown volume / m$^3$m$^{-2}$',fontsize=axis_size)
    ax2.set_ylabel('PAI',fontsize=axis_size)
    """
    for k in range(1,3):
        for i in range(0,N_plots):
            ax2.plot(inventory_LAI[Plots[i]],radiative_LAI[Plots[i]][:,k],'.',color=colour[k],alpha=0.5)
    """
    
    m, c, r2, p, x_i, y_i, CI_u, CI_l, PI_u, PI_l = linear_regression(vol,rad2,conf=0.95)
    ax2.fill_between(x_i,CI_l,CI_u,color='0.75')
    ax2.plot(x_i,y_i,'-',color='black')
    
    for k in range(1,3):
        for i in range(0,N_plots):
            x_err=np.nan#np.std(inventory_LAI[Plots[i]])/np.sqrt(n_subplots)
            y_err=np.nan#np.std(radiative_LAI[Plots[i]][:,k])/np.sqrt(n_subplots)
            if i==0:
                leg_label = '$k_{max}=$' + str(k+1) 
                #ax2.errorbar(np.mean(inventory_LAI[Plots[i]]),np.mean(radiative_LAI[Plots[i]][:,k]),xerr=x_err,yerr=y_err,marker='o',color=colour[k],label=leg_label)
                ax2.errorbar(np.mean(Inventory_LAI[Plots[i]]), np.nansum(radiative_LAD_mean[Plots[i]][:,k])*layer_thickness,xerr=x_err,yerr=y_err,marker='o',color=colour[k],label=leg_label)
            else:
                #ax2.errorbar(np.mean(inventory_LAI[Plots[i]]),np.mean(radiative_LAI[Plots[i]][:,k]),xerr=x_err,yerr=y_err,marker='o',color=colour[k])
                ax2.errorbar(np.mean(Inventory_LAI[Plots[i]]), np.nansum(radiative_LAD_mean[Plots[i]][:,k])*layer_thickness,xerr=x_err,yerr=y_err,marker='o',color=colour[k])

    ax3 = plt.subplot2grid((1,3),(0,2), sharex=ax1)#, sharey=ax1)
    ax3.annotate('c - radiative transfer (modified)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=axis_size)
    ax3.annotate('$k_{max}=2$; ' + r_sq_c[0] + '\n$k_{max}=3$; ' + r_sq_c[1], xy=(0.95,0.05), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='bottom', fontsize=axis_size)
    ax3.set_ylabel('PAI',fontsize=axis_size)
    ax3.set_xlabel('crown volume / m$^3$m$^{-2}$',fontsize=axis_size)
    """
    for k in range(1,3):
        for i in range(0,N_plots):
            ax3.plot(inventory_LAI[Plots[i]],radiative_DTM_LAI[Plots[i]][:,k],'.',color=colour[k],alpha=0.5)
    """
    
    m, c, r2, p, x_i, y_i, CI_u, CI_l, PI_u, PI_l = linear_regression(vol,radDTM2,conf=0.95)
    ax3.fill_between(x_i,CI_l,CI_u,color='0.75')
    ax3.plot(x_i,y_i,'-',color='black')
    
    for k in range(1,3):
        for i in range(0,N_plots):
            x_err=np.nan#np.std(inventory_LAI[Plots[i]])/np.sqrt(n_subplots)
            y_err=np.nan#np.std(radiative_DTM_LAI[Plots[i]][:,k])/np.sqrt(n_subplots)
            if i==0:
                leg_label = '$k_{max}=$' + str(k+1) 
                #ax3.errorbar(np.mean(inventory_LAI[Plots[i]]),np.mean(radiative_DTM_LAI[Plots[i]][:,k]),xerr=x_err,yerr=y_err,marker='o',color=colour[k],label=leg_label)
                ax3.errorbar(np.mean(Inventory_LAI[Plots[i]]), np.nansum(radiative_DTM_LAD_mean[Plots[i]][:,k])*layer_thickness,xerr=x_err,yerr=y_err,marker='o',color=colour[k],label=leg_label)
            else:
                #ax3.errorbar(np.mean(inventory_LAI[Plots[i]]),np.mean(radiative_DTM_LAI[Plots[i]][:,k]),xerr=x_err,yerr=y_err,marker='o',color=colour[k])
                ax3.errorbar(np.mean(Inventory_LAI[Plots[i]]), np.nansum(radiative_DTM_LAD_mean[Plots[i]][:,k])*layer_thickness,xerr=x_err,yerr=y_err,marker='o',color=colour[k])
                
    ax3.legend(loc=4)

    #ax1.set_ylim((0,20))
    plt.tight_layout()
    plt.savefig(figure_name)
    plt.show()
    return 0

# LAI vs. basal area
def plot_LAI_vs_basal_area(figure_name,figure_number,MacArthurHorn_LAD,MacArthurHorn_LAD_mean,radiative_LAD,radiative_LAD_mean,radiative_DTM_LAD,radiative_DTM_LAD_mean,BasalArea,layer_thickness=1):


    Plots = MacArthurHorn_LAD.keys()
    N_plots = len(Plots)
    N_subplots = np.shape(MacArthurHorn_LAD[Plots[0]])[0]
    BA = np.zeros(N_plots)
    mh = np.zeros(N_plots)
    rad2 = np.zeros(N_plots)
    rad3 = np.zeros(N_plots)
    radDTM2 = np.zeros(N_plots)
    radDTM3 = np.zeros(N_plots)

    for pp in range(0,N_plots):
        print Plots[pp]
        # for 1 ha averages, want to avoid nodata
        # these should be identifiable based on penetration limit
        BA[pp] = np.mean(BasalArea[Plots[pp]]).copy()
        mh[pp] = np.nansum(MacArthurHorn_LAD_mean[Plots[pp]])*layer_thickness
        rad2[pp] = np.nansum(radiative_LAD_mean[Plots[pp]][:,1])*layer_thickness
        rad3[pp] = np.nansum(radiative_LAD_mean[Plots[pp]][:,2])*layer_thickness
        radDTM2[pp] = np.nansum(radiative_DTM_LAD_mean[Plots[pp]][:,1])*layer_thickness
        radDTM3[pp] = np.nansum(radiative_DTM_LAD_mean[Plots[pp]][:,2])*layer_thickness
    
    # annotate with stats
    r_sq_a =   aux.get_rsquared_annotation(BA,mh)
    r_sq_b = aux.get_rsquared_annotation(BA,radDTM3)

    plt.figure(figure_number, facecolor='White',figsize=[8,4])
    ax1 = plt.subplot2grid((1,2),(0,0))
    ax1.set_xlabel('Basal Area / m$^2$ha$^{-1}$',fontsize=axis_size)
    ax1.set_ylabel('PAI',fontsize=axis_size)
    ax1.annotate('a - MacArthur-Horn', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=axis_size)
    ax1.annotate(r_sq_a, xy=(0.95,0.05), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='bottom', fontsize=axis_size)

    # plot linear regression
    m, c, r2, p, x_i, y_i, CI_u, CI_l, PI_u, PI_l = linear_regression(BA,mh,conf=0.95)
    ax1.fill_between(x_i,CI_l,CI_u,color='0.75')
    ax1.plot(x_i,y_i,'-',color='black')
    for i in range(0,N_plots):
        ax1.plot(BA[i],mh[i],marker=plot_marker[Plots[i]],markerfacecolor=plot_colour[Plots[i]],linestyle='None')

    ax3 = plt.subplot2grid((1,2),(0,1), sharex=ax1)#, sharey=ax1)
    ax3.annotate('b - radiative transfer (modified)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=axis_size)
    ax3.annotate(r_sq_b, xy=(0.95,0.05), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='bottom', fontsize=axis_size)
    ax3.set_ylabel('PAI',fontsize=axis_size)
    ax3.set_xlabel('Basal Area / m$^2$ha$^{-1}$',fontsize=axis_size)
    
    m, c, r2, p, x_i, y_i, CI_u, CI_l, PI_u, PI_l = linear_regression(BA,radDTM3,conf=0.95)
    ax3.fill_between(x_i,CI_l,CI_u,color='0.75')
    ax3.plot(x_i,y_i,'-',color='black')                                                          
    for i in range(0,N_plots):
        ax3.plot(BA[i],radDTM3[i],marker=plot_marker[Plots[i]],markerfacecolor=plot_colour[Plots[i]],linestyle='None')

    ax3.legend(loc=4)

    plt.tight_layout()
    plt.savefig(figure_name)
    plt.show()
    return 0



# Location map
def plot_location_map(figure_name,figure_number):
    
    # create simple colourmap
    cmap_colours = [(1.0,  1.0, 1.0),(170./256.,229./256.,156./256.)]
    cbins= 100
    cmap_name = 'w_gn'
    white_green = LinearSegmentedColormap.from_list(cmap_name, cmap_colours,N=cbins)

    #loc_cb = plticker.MultipleLocator(base=0.2) 
    """
    GFW_file = '/disk/scratch/local.2/dmilodow/GFW/Hansen_GFC2015_treecover2000_10N_110E_resampled.tif'
    mask_file = '/disk/scratch/local.2/dmilodow/GFW/Hansen_GFC2015_datamask_10N_110E_resampled.tif'
    mask = gdal.Open(mask_file).ReadAsArray()
    ds = gdal.Open(GFW_file)
    """
    AGB_file = '/home/dmilodow/DataStore_GCEL/AGB/avitabile/Avitabile_AGB_Map/Avitabile_AGB_Map.tif'
    ds = gdal.Open(AGB_file)
    geoTrans = ds.GetGeoTransform()
    W = 114
    E = 120
    N = 9
    S = 3
    agb,geoT=geo.clip_array_to_bbox(np.asarray(ds.ReadAsArray(),dtype='float'),geoTrans,N,S,W,E)
    agb[agb<0]=np.nan
    array = ma.masked_invalid(agb)
    #array = np.asarray(ds.ReadAsArray(),dtype='float')/100.
    #array[mask!=1]=np.nan
    xres = geoT[1]
    yres = geoT[5]
    x0 = geoT[0]
    y0 = geoT[3]
    rows,cols = array.shape
    lon = np.arange(x0,x0+xres*cols,xres)+ xres * 0.5
    lat = np.arange(y0,y0+yres*rows,yres)+ yres * 0.5
    lons,lats = np.meshgrid(lon, lat)
    plt.figure(figure_number, facecolor='White',figsize=[10,7])
    # Region
    ax1 = plt.subplot2grid((1,1),(0,0))

    m = Basemap(projection='aea', lat_0=5.5, lon_0=116.8, llcrnrlat=3.9, urcrnrlat=7.5,llcrnrlon=115.2, urcrnrlon=119.5, resolution='i')
    m.ax = ax1

    # plot some lines
    parallels = np.arange(0.,12.,2.)
    # labels = [left,right,top,bottom]
    m.drawparallels(parallels,labels=[True,False,False,False],fontsize=8)
    meridians = np.arange(110.,122.,2.)
    m.drawmeridians(meridians,labels=[False,False,False,True],fontsize=8)

    #im = m.pcolormesh(xx,yy, array, cmap=white_green, vmin=0.3, vmax=1)
    im = m.pcolormesh(lons,lats, array, cmap='viridis', latlon=True)

    country_shp = '/home/dmilodow/DataStore_DTM/EOlaboratory/Areas/NaturalEarth/10m_cultural/ne_10m_admin_0_countries'
    shp = m.readshapefile(country_shp, 'countries', drawbounds=False)
    for info, shape in zip(m.countries_info, m.countries):
        if info['name'] == 'Indonesia':
            poly = Polygon(shape, facecolor='0.2',edgecolor='0.2',alpha = 0.8)
            m.ax.add_patch(poly)
        #elif info['name'] == 'Malaysia':
        #    poly = Polygon(shape, facecolor='None',edgecolor='0.5')
        #    m.ax.add_patch(poly)
        elif info['name'] == 'Brunei':
            poly = Polygon(shape, facecolor='0.2',edgecolor='0.2',alpha = 0.8)
            m.ax.add_patch(poly)

    # Plot locations
    plot_shp = 'plot_locations_for_location_map_WGS84'
    pshp = m.readshapefile(plot_shp, 'plot', drawbounds=False)
    for info, point in zip(m.plot_info, m.plot):
        if info['ForestType'] == 'OG':
            col = colour[0]
        elif info['ForestType'] == 'ML':
            col = colour[2]
        else:
            col=colour[1]
        m.ax.plot(point[0], point[1], marker='o', color=col, markersize=15)
    
    # colorbar bits and bobs
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar=plt.colorbar(im, cax=cax)#,extend = 'min')
    #cbar.ax.set_ylabel('fraction tree cover in 2000',fontsize = axis_size)
    cbar.ax.set_ylabel('AGB / Mg ha$^{-1}$',fontsize = 15)
    cbar.solids.set_edgecolor("face")
    #cbar.locator = loc_cb
    #cbar.update_ticks()

    ax1.annotate('Sabah\nMalaysia', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=20, fontweight='bold')

    ax1.annotate('SAFE', xy=(277058,89946), xycoords='data',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=20,color='white', fontweight='bold')
    ax1.annotate('Maliau\nBasin', xy=(186070,100167), xycoords='data',backgroundcolor='none',horizontalalignment='right', verticalalignment='bottom', fontsize=20,color='white', fontweight='bold')
    ax1.annotate('Danum Valley', xy=(290346,122654), xycoords='data',backgroundcolor='none',horizontalalignment='left', verticalalignment='bottom', fontsize=20,color='white', fontweight='bold')
    plt.tight_layout()
    plt.savefig(figure_name)
    #plt.show()
        
    return 0


# Cross-plot canopy layers (new)
def cross_plot_canopy_layers_LiDAR(figure_name,figure_number,heights,heights_rad,MacArthurHorn_LAD,MacArthurHorn_LAD_mean,radiative_LAD,radiative_LAD_mean,radiative_DTM_LAD,radiative_DTM_LAD_mean,inventory_LAD):

    # Process the profiles
    
    # 2 columns (plotting MacHorn on x axis, multi return on y axis)
    # 3 rows (old growth, moderately logged, heavily logged)
    fig = plt.figure(figure_number, facecolor='White',figsize=[6,9])

    # First up - old-growth forest
    Plot_name='Belian'
    x=np.cumsum(MacArthurHorn_LAD_mean[Plot_name][2:][::-1])
    y1=np.cumsum(radiative_LAD_mean[Plot_name][:-3,1])
    y2=np.cumsum(radiative_DTM_LAD_mean[Plot_name][:-3,1])
    t=heights[2:][::-1]
    ax1 = plt.subplot2grid((3,2),(0,0))    
    ax1.annotate('a', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)
    for i in range(0,n_subplots):
        ax1.plot(np.cumsum(MacArthurHorn_LAD[Plot_name][i,2:][::-1]),np.cumsum(radiative_LAD[Plot_name][i,:-3,1]),'-',color='k',linewidth=0.5)
    plot_colourline(ax1,x,y1,t,linewidth=5,cmap='plasma')

    ax2 = plt.subplot2grid((3,2),(0,1),sharex=ax1,sharey=ax1)
    ax2.annotate('b', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)
    for i in range(0,n_subplots):
        ax2.plot(np.cumsum(MacArthurHorn_LAD[Plot_name][i,2:][::-1]),np.cumsum(radiative_DTM_LAD[Plot_name][i,:-3,1]),'-',color='k',linewidth=0.5)
    plot_colourline(ax2,x,y2,t,linewidth=5,cmap='plasma')

    # Next up - Moderately logged forest
    Plot_name='E'
    x=np.cumsum(MacArthurHorn_LAD_mean[Plot_name][2:][::-1])
    y1=np.cumsum(radiative_LAD_mean[Plot_name][:-3,1])
    y2=np.cumsum(radiative_DTM_LAD_mean[Plot_name][:-3,1])
    ax3 = plt.subplot2grid((3,2),(1,0),sharex=ax1,sharey=ax1)
    ax3.annotate('c', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)
    ax3.set_ylabel('cumulative PAD$_{rad trans (Detto)}$ (m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size)
    
    for i in range(0,n_subplots):
        ax3.plot(np.cumsum(MacArthurHorn_LAD[Plot_name][i,2:][::-1]),np.cumsum(radiative_LAD[Plot_name][i,:-3,1]),'-',color='k',linewidth=0.5)
    plot_colourline(ax3,x,y1,t,linewidth=5,cmap='plasma')
    
    ax4 = plt.subplot2grid((3,2),(1,1),sharex=ax1,sharey=ax1)
    ax4.annotate('d', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)
    ax4.set_ylabel('cumulative PAD$_{rad trans (new)}$ (m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size)
    ax4.yaxis.set_label_position("right")
    ax4.yaxis.tick_right()
    ax4.yaxis.set_ticks_position('both')

    for i in range(0,n_subplots):
        ax4.plot(np.cumsum(MacArthurHorn_LAD[Plot_name][i,2:][::-1]),np.cumsum(radiative_DTM_LAD[Plot_name][i,:-3,1]),'-',color='k',linewidth=0.5)
    
    plot_colourline(ax4,x,y2,t,linewidth=5,cmap='plasma')

    # Finally heavily logged forest
    Plot_name = 'B North'
    x=np.cumsum(MacArthurHorn_LAD_mean[Plot_name][2:][::-1])
    y1=np.cumsum(radiative_LAD_mean[Plot_name][:-3,1])
    y2=np.cumsum(radiative_DTM_LAD_mean[Plot_name][:-3,1])
    ax5 = plt.subplot2grid((3,2),(2,0),sharex=ax1,sharey=ax1)    
    ax5.annotate('e', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)
    
    for i in range(0,n_subplots):
        ax5.plot(np.cumsum(MacArthurHorn_LAD[Plot_name][i,2:][::-1]),np.cumsum(radiative_LAD[Plot_name][i,:-3,1]),'-',color='k',linewidth=0.5)
    plot_colourline(ax5,x,y1,t,linewidth=5,cmap='plasma')
    
    ax6 = plt.subplot2grid((3,2),(2,1),sharex=ax1,sharey=ax1)
    ax6.annotate('f', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)

    for i in range(0,n_subplots):
        ax6.plot(np.cumsum(MacArthurHorn_LAD[Plot_name][i,2:][::-1]),np.cumsum(radiative_DTM_LAD[Plot_name][i,:-3,1]),'-',color='k',linewidth=0.5)
        
    plot_colourline(ax6,x,y2,t,linewidth=5,cmap='plasma')

    # tidy up tick labels
    xticklabels=[ax1.get_xticklabels() + ax2.get_xticklabels() + ax3.get_xticklabels() + ax4.get_xticklabels()]
    plt.setp(xticklabels,visible=False)

    plot_one2one(ax1)
    plot_one2one(ax2)
    plot_one2one(ax3)
    plot_one2one(ax4)
    plot_one2one(ax5)
    plot_one2one(ax6)

    # Create dummy subplot for common axis labels
    ax_ = fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    ax_.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    ax_.set_xlabel('cumulative PAD$_{MacArthur-Horn}$ (m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size, labelpad=10) # Use argument `labelpad` to move label downwards.

    # plot colorbar in ax1
    cmap = cmaps.plasma
    norm = colors.Normalize(vmin=0., vmax=80.)
    cbaxes = inset_axes(ax1, width="5%", height="50%", loc=6) 
    cbar= colorbar.ColorbarBase(cbaxes, cmap=cmap,norm=norm, ticks=[0.,40.,80.], orientation='vertical')
    cbar.set_label('Height (m)')
    plt.tight_layout()
    plt.show()
    
    
    return 0


# Plot canopy residuals

def plot_canopy_layer_residuals(figure_name,figure_number,heights,MacArthurHorn_LAD,MacArthurHorn_LAD_mean,radiative_LAD,radiative_LAD_mean,radiative_DTM_LAD,radiative_DTM_LAD_mean,max_return=2):
    
    # 2 columns (plotting MacHorn on x axis, multi return on y axis)
    # 3 rows (old growth, moderately logged, heavily logged)
    fig = plt.figure(figure_number, facecolor='White',figsize=[6,9])

    # First up - old-growth forest
    Plot_name='Belian'
    ids = np.arange(heights.size)
    idmax = ids[MacArthurHorn_LAD_mean[Plot_name]>0][-1]+1
    y=heights[2:idmax]
    x=MacArthurHorn_LAD_mean[Plot_name][2:idmax]
    x1=radiative_LAD_mean[Plot_name][:-1,max_return-1][::-1][2:idmax]
    x2=radiative_DTM_LAD_mean[Plot_name][:-1,max_return-1][::-1][2:idmax]

    ax1 = plt.subplot2grid((3,2),(0,0))    
    ax1.annotate('a', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)
    for i in range(0,n_subplots):
        ax1.plot((radiative_LAD[Plot_name][i,:-1,max_return-1][::-1][2:idmax]-MacArthurHorn_LAD[Plot_name][i,2:idmax]),y,'-',color='0.5',linewidth=0.5,alpha=0.5)
    ax1.plot(x1-x,y,linewidth=2,color='k')

    ax2 = plt.subplot2grid((3,2),(0,1),sharex=ax1,sharey=ax1)
    ax2.annotate('b', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)
    ax2.yaxis.tick_right()
    ax2.yaxis.set_ticks_position('both')
    for i in range(0,n_subplots):
        ax2.plot((radiative_DTM_LAD[Plot_name][i,:-1,max_return-1][::-1][2:idmax]-MacArthurHorn_LAD[Plot_name][i,2:idmax]),y,'-',color='0.5',linewidth=0.5,alpha=0.5)
    ax2.plot(x2-x,y,linewidth=2,color='k')

    # Next up - Moderately logged forest
    Plot_name='E'
    ids = np.arange(heights.size)
    idmax = ids[MacArthurHorn_LAD_mean[Plot_name]>0][-1]+1
    y=heights[2:idmax]
    x=MacArthurHorn_LAD_mean[Plot_name][2:idmax]
    x1=radiative_LAD_mean[Plot_name][:-1,max_return-1][::-1][2:idmax]
    x2=radiative_DTM_LAD_mean[Plot_name][:-1,max_return-1][::-1][2:idmax]
    ax3 = plt.subplot2grid((3,2),(1,0),sharex=ax1,sharey=ax1)
    ax3.annotate('c', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)
    ax3.set_ylabel('height (m)',fontsize=axis_size)
    
    for i in range(0,n_subplots):
        ax3.plot((radiative_LAD[Plot_name][i,:-1,max_return-1][::-1][2:idmax]-MacArthurHorn_LAD[Plot_name][i,2:idmax]),y,'-',color='0.5',linewidth=0.5,alpha=0.5)
    ax3.plot(x1-x,y,linewidth=2,color='k')
    
    ax4 = plt.subplot2grid((3,2),(1,1),sharex=ax1,sharey=ax1)
    ax4.annotate('d', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)
    ax4.set_ylabel('height (m))',fontsize=axis_size)
    ax4.yaxis.set_label_position("right")
    ax4.yaxis.tick_right()
    ax4.yaxis.set_ticks_position('both')

    for i in range(0,n_subplots):
        ax4.plot((radiative_DTM_LAD[Plot_name][i,:-1,max_return-1][::-1][2:idmax]-MacArthurHorn_LAD[Plot_name][i,2:idmax]),y,'-',color='0.5',linewidth=0.5,alpha=0.5)
    ax4.plot(x2-x,y,linewidth=2,color='k')

    # Finally heavily logged forest
    Plot_name = 'B South'
    ids = np.arange(heights.size)
    idmax = ids[MacArthurHorn_LAD_mean[Plot_name]>0][-1]+1
    y=heights[2:idmax]
    x=MacArthurHorn_LAD_mean[Plot_name][2:idmax]
    x1=radiative_LAD_mean[Plot_name][:-1,max_return-1][::-1][2:idmax]
    x2=radiative_DTM_LAD_mean[Plot_name][:-1,max_return-1][::-1][2:idmax]
    ax5 = plt.subplot2grid((3,2),(2,0),sharex=ax1,sharey=ax1)    
    ax5.annotate('e', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)
    ax5.set_xlabel('PAD$_{rad trans-Detto}$-PAD$_{MacArthur-Horn}$\n(m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size)
    
    for i in range(0,n_subplots):
        ax5.plot((radiative_LAD[Plot_name][i,:-1,max_return-1][::-1][2:idmax]-MacArthurHorn_LAD[Plot_name][i,2:idmax]),y,'-',color='0.5',linewidth=0.5,alpha=0.5)
    ax5.plot(x1-x,y,linewidth=2,color='k')
    
    ax6 = plt.subplot2grid((3,2),(2,1),sharex=ax1,sharey=ax1)
    ax6.annotate('f', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)
    ax6.set_xlabel('PAD$_{rad trans-new}$-PAD$_{MacArthur-Horn}$\n(m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size)
    ax6.yaxis.tick_right()
    ax6.yaxis.set_ticks_position('both')
    
    for i in range(0,n_subplots):
        ax6.plot((radiative_DTM_LAD[Plot_name][i,:-1,max_return-1][::-1][2:idmax]-MacArthurHorn_LAD[Plot_name][i,2:idmax]),y,'-',color='0.5',linewidth=0.5,alpha=0.5)
    ax6.plot(x2-x,y,linewidth=2,color='k')

    # tidy up tick labels
    xticklabels=[ax1.get_xticklabels() + ax2.get_xticklabels() + ax3.get_xticklabels() + ax4.get_xticklabels()]
    plt.setp(xticklabels,visible=False)

    plt.tight_layout()
    plt.show()
    
    
# plot PAD-volume ratios
def plot_canopy_layer_PAD_volume_ratio(figure_name,figure_number,heights,MacArthurHorn_LAD,MacArthurHorn_LAD_mean,radiative_LAD,radiative_LAD_mean,radiative_DTM_LAD,radiative_DTM_LAD_mean,inventory_LAD,max_return=2):
    fig = plt.figure(figure_number, facecolor='White',figsize=[8,9])

    # First up - old-growth forest
    Plot_name='Belian'
    ids = np.arange(heights.size)
    idmax = ids[MacArthurHorn_LAD_mean[Plot_name]>0][-1]+1
    y=heights[2:idmax]
    x=np.mean(inventory_LAD[Plot_name],axis=0)[2:idmax]
    x1=MacArthurHorn_LAD_mean[Plot_name][2:idmax]
    x2=radiative_LAD_mean[Plot_name][:-1,max_return-1][::-1][2:idmax]
    x3=radiative_DTM_LAD_mean[Plot_name][:-1,max_return-1][::-1][2:idmax]

    ax1 = plt.subplot2grid((3,3),(0,0))    
    ax1.annotate('a', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)
    ax1.plot(np.log10(x1/x),y,linewidth=2,color='k')

    ax2 = plt.subplot2grid((3,3),(0,1),sharex=ax1,sharey=ax1)
    ax2.annotate('b', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)
    ax2.plot(np.log10(x2/x),y,linewidth=2,color='k')

    ax3 = plt.subplot2grid((3,3),(0,2),sharex=ax1,sharey=ax1)
    ax3.annotate('c', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)
    ax3.yaxis.tick_right()
    ax3.yaxis.set_ticks_position('both')
    ax3.plot(np.log10(x2/x),y,linewidth=2,color='k')
    
    # Next up - Moderately logged forest
    Plot_name='E'
    ids = np.arange(heights.size)
    idmax = ids[MacArthurHorn_LAD_mean[Plot_name]>0][-1]+1
    y=heights[2:idmax]
    x=np.mean(inventory_LAD[Plot_name],axis=0)[2:idmax]
    x1=MacArthurHorn_LAD_mean[Plot_name][2:idmax]
    x2=radiative_LAD_mean[Plot_name][:-1,max_return-1][::-1][2:idmax]
    x3=radiative_DTM_LAD_mean[Plot_name][:-1,max_return-1][::-1][2:idmax]

    ax4 = plt.subplot2grid((3,3),(1,0),sharex=ax1,sharey=ax1)
    ax4.annotate('d', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)
    ax4.set_ylabel('height (m)',fontsize=axis_size)
    ax4.plot(np.log10(x1/x),y,linewidth=2,color='k')
    
    ax5 = plt.subplot2grid((3,3),(1,1),sharex=ax1,sharey=ax1)
    ax5.annotate('e', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)
    ax5.plot(np.log10(x2/x),y,linewidth=2,color='k')

    ax6 = plt.subplot2grid((3,3),(1,2),sharex=ax1,sharey=ax1)
    ax6.annotate('f', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)
    ax6.set_ylabel('height (m))',fontsize=axis_size)
    ax6.yaxis.set_label_position("right")
    ax6.yaxis.tick_right()
    ax6.yaxis.set_ticks_position('both')

    ax6.plot(np.log10(x2/x),y,linewidth=2,color='k')
    # Finally heavily logged forest
    Plot_name = 'B South'
    ids = np.arange(heights.size)
    idmax = ids[MacArthurHorn_LAD_mean[Plot_name]>0][-1]+1
    y=heights[2:idmax]
    x=np.mean(inventory_LAD[Plot_name],axis=0)[2:idmax]
    x1=MacArthurHorn_LAD_mean[Plot_name][2:idmax]
    x2=radiative_LAD_mean[Plot_name][:-1,max_return-1][::-1][2:idmax]
    x3=radiative_DTM_LAD_mean[Plot_name][:-1,max_return-1][::-1][2:idmax]

    ax7 = plt.subplot2grid((3,3),(2,0),sharex=ax1,sharey=ax1)    
    ax7.annotate('g', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)
    ax7.set_xlabel('PAD$_{rad trans-Detto}$-PAD$_{MacArthur-Horn}$\n(m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size)
    ax7.plot(np.log10(x1/x),y,linewidth=2,color='k')
    
    ax8 = plt.subplot2grid((3,3),(2,1),sharex=ax1,sharey=ax1)
    ax8.annotate('h', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)
    ax8.set_xlabel('PAD$_{rad trans-new}$-PAD$_{MacArthur-Horn}$\n(m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size)
    ax8.plot(np.log10(x2/x),y,linewidth=2,color='k')
    
    ax9 = plt.subplot2grid((3,3),(2,2),sharex=ax1,sharey=ax1)
    ax9.annotate('i', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)
    ax9.set_xlabel('PAD$_{rad trans-new}$-PAD$_{MacArthur-Horn}$\n(m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size)
    ax9.yaxis.tick_right()
    ax9.yaxis.set_ticks_position('both')
    ax9.plot(np.log10(x2/x),y,linewidth=2,color='k')

    # tidy up tick labels
    xticklabels=[ax1.get_xticklabels() + ax2.get_xticklabels() + ax3.get_xticklabels() + ax4.get_xticklabels()]
    plt.setp(xticklabels,visible=False)

    plt.tight_layout()
    plt.show()
    
# Plot allometric relationships
# Plot three subplots defining the relationships used to construct the canopy profiles
def plot_allometric_relationships(figure_name,figure_number,field_file,allometry_file):

    DBH_BAAD, H_BAAD, D_BAAD = field.load_BAAD_crown_allometry_data(allometry_file)
    a, b, CF, r_sq, p, H_, PI_u, PI_l = field.log_log_linear_regression(H_BAAD,D_BAAD,conf=0.90)
    
    field_data = field.load_crown_survey_data(field_file)
    a_Ht, b_Ht, CF_Ht, r_sq_Ht, p_Ht, DBH_Ht, PI_u_Ht, PI_l_Ht = field.log_log_linear_regression(field_data['DBH_field'],field_data['Height_field'],conf=0.90)
    
    a_A, b_A, CF_A, r_sq_A, p_A, DBH_A, PI_u_A, PI_l_A = field.log_log_linear_regression(field_data['DBH_field'],field_data['CrownArea'],conf=0.90)

    PI_25_D,PI_75_D = field.calculate_powerlaw_prediction_interval_bootstrap_resampling_residuals(H_,H_BAAD,D_BAAD,conf=0.5,niter=10000)
    PI_05_D,PI_95_D = field.calculate_powerlaw_prediction_interval_bootstrap_resampling_residuals(H_,H_BAAD,D_BAAD,conf=0.9,niter=10000)

    mask = np.all((np.isfinite(field_data['DBH_field']),np.isfinite(field_data['CrownArea'])),axis=0)
    PI_25_A,PI_75_A = field.calculate_powerlaw_prediction_interval_bootstrap_resampling_residuals(DBH_A,field_data['DBH_field'][mask],field_data['CrownArea'][mask],conf=0.5,niter=10000)
    PI_05_A,PI_95_A = field.calculate_powerlaw_prediction_interval_bootstrap_resampling_residuals(DBH_A,field_data['DBH_field'][mask],field_data['CrownArea'][mask],conf=0.9,niter=10000)
    
    mask = np.all((np.isfinite(field_data['DBH_field']),np.isfinite(field_data['Height_field'])),axis=0)
    PI_25_H,PI_75_H = field.calculate_powerlaw_prediction_interval_bootstrap_resampling_residuals(DBH_Ht,field_data['DBH_field'][mask],field_data['Height_field'][mask],conf=0.5,niter=10000)
    PI_05_H,PI_95_H = field.calculate_powerlaw_prediction_interval_bootstrap_resampling_residuals(DBH_Ht,field_data['DBH_field'][mask],field_data['Height_field'][mask],conf=0.9,niter=10000)
    
    D_ = CF*a*H_**b
    CA_ = CF_A*a_A*DBH_A**b_A
    Ht_ = CF_Ht*a_Ht*DBH_Ht**b_Ht
    
    plt.figure(figure_number,facecolor='White',figsize=[9,3])
    ax1 = plt.subplot2grid((1,3),(0,0))
    ax1.set_ylabel('Crown Depth / m')
    ax1.set_xlabel('Height / cm')
    ax1.annotate('a - Height-Crown Depth', xy=(0.08,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    
    ax1.fill_between(H_,PI_05_D,PI_95_D,color='0.95')
    ax1.fill_between(H_,PI_25_D,PI_75_D,color='0.85')
    ax1.plot(H_BAAD,D_BAAD,'.',color='#1A2BCE',alpha=0.2)
    ax1.plot(H_,D_,'-',color='black')
    
    eq = '$D=%.2fH^{%.2f}$\n$r^2=%.2f$' % (CF*a, b,r_sq)
    ax1.annotate(eq, xy=(0.08,0.85), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    
    ax2 = plt.subplot2grid((1,3),(0,1))
    ax2.annotate('b - DBH-Crown Area', xy=(0.08,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax2.set_xlabel('DBH / cm')
    ax2.set_ylabel('Crown Area / m$^2$')

    ax2.fill_between(DBH_A,PI_05_A,PI_95_A,color='0.95')
    ax2.fill_between(DBH_A,PI_25_A,PI_75_A,color='0.85')    
    ax2.plot(field_data['DBH_field'],field_data['CrownArea'],'.',color='#1A2BCE',alpha=0.2)
    ax2.plot(DBH_A,CA_,'-',color='black')
    
    eq = '$A=%.2fDBH^{%.2f}$\n$r^2=%.2f$' % (CF_A*a_A, b_A,r_sq_A)
    ax2.annotate(eq, xy=(0.08,0.85), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax2.axvspan(0,10,color='0.8')
    
    ax3 = plt.subplot2grid((1,3),(0,2))#,sharex=ax3b)
    ax3.annotate('c - DBH-Height', xy=(0.08,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax3.set_xlabel('DBH / cm')
    ax3.set_ylabel('Height / m')

    ax3.fill_between(DBH_Ht,PI_05_H,PI_95_H,color='0.95')
    ax3.fill_between(DBH_Ht,PI_25_H,PI_75_H,color='0.85')    
    ax3.plot(field_data['DBH_field'],field_data['Height_field'],'.',color='#1A2BCE',alpha=0.2)
    ax3.plot(DBH_Ht,Ht_,'-',color='black')
    
    eq = '$H=%.2fDBH^{%.2f}$\n$r^2=%.2f$' % (CF_Ht*a_Ht, b_Ht,r_sq_Ht)
    ax3.annotate(eq, xy=(0.08,0.85), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax3.axvspan(0,10,color='0.8')


    plt.tight_layout()
    plt.show()
    return 0
