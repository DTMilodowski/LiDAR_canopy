import LiDAR_io as io
import LiDAR_tools as lidar
import auxilliary_functions as aux
import LiDAR_MacHorn_LAD_profiles as LAD1
import LiDAR_radiative_transfer_LAD_profiles as LAD2
import inventory_based_LAD_profiles as field
import least_squares_fitting as lstsq
import geospatial_utility_tools as geo
import structural_metrics as metrics

import sys
import numpy as np
from scipy import stats
import pandas as pd
import seaborn as sns
from osgeo import gdal
import numpy.ma as ma

import cartopy
from cartopy import config
import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from matplotlib import pyplot as plt
from matplotlib import rcParams
import matplotlib.ticker as mticker
import matplotlib.ticker as plticker
from matplotlib import colorbar
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Set up some basiic parameters for the plots
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = 8
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

"""
# Plot point cloud example for introduction
"""
# plot point clouds with canopy profiles
def plot_point_cloud(figure_name,figure_number, gps_pts_file,plot_point_cloud):

    max_return=3
    colour = ['#46E900','#1A2BCE','#E0007F']
    rgb = [[70,233,0],[26,43,206],[224,0,127]]
    labels = ['$1^{st}$', '$2^{nd}$', '$3^{rd}$', '$4^{th}$']

    # first up, going to need to find affine transformation to rotate the point cloud for
    # easy plotting

    # The GPS coordinates for the plot locations can be used to find this transformation matrix
    datatype = {'names': ('plot', 'x', 'y', 'x_prime', 'y_prime'), 'formats': ('S32','f16','f16','f16','f16')}
    plot_coords = np.genfromtxt(gps_pts_file, delimiter = ',',dtype=datatype)
    plot_coords['plot'][plot_coords['plot']==b'maliau_belian'] = b'Belian'
    Plot = b'Belian'

    # first get points for a given plot and build matrices - note that I've reversed xy and xy_prime in this case to reverse the rotation-translation
    mask = plot_coords['plot']==Plot
    x = plot_coords['x_prime'][mask]
    y = plot_coords['y_prime'][mask]
    x_prime = plot_coords['x'][mask]
    y_prime = plot_coords['y'][mask]
    affine=lstsq.least_squares_affine_matrix(x,y,x_prime,y_prime)

    Xi = np.asarray([plot_point_cloud[Plot][:,0],plot_point_cloud[Plot][:,1],
                            np.ones(plot_point_cloud[Plot].shape[0])])
    Xi_prime = np.dot(affine,Xi)

    plot_point_cloud_display = plot_point_cloud[Plot].copy()
    plot_point_cloud_display[:,0]=Xi_prime[0]
    plot_point_cloud_display[:,1]=Xi_prime[1]

    plt.figure(figure_number, facecolor='White',figsize=[5,4])

    # Belian
    ax = plt.subplot2grid((1,1),(0,0))
    ax.annotate('Old growth forest, MLA01', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax.set_ylabel('Height / m',fontsize=axis_size)
    plt.gca().set_aspect('equal', adjustable='box-forced')

    plot_lidar_pts = plot_point_cloud_display
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
        ax.scatter(points_y,points_z,marker='o',c=colours,edgecolors='none',s=2)
        ax.scatter(0,0,marker='o',c=colours[0,0:3],edgecolors='none',s=2,
                        label=labels[k])

    ax.set_xlim(0,100)
    ax.set_ylim(0,80)
    plt.tight_layout()
    plt.savefig(figure_name)
    plt.show()
    return 0

"""
# plot point clouds with canopy profiles
# six rows for six plots (2x old growth, 2x moderately logged, 2x heavily logged)
# In old version, plot the Detto model without correction factor as well as the
# other two
"""
def plot_point_clouds_and_profiles_old(figure_name,figure_number, gps_pts_file,
                            plot_point_cloud,heights,heights_rad, lidar_profiles,
                            MacArthurHorn_PAD,MacArthurHorn_PAD_mean,radiative_PAD,radiative_PAD_mean,
                            radiative_DTM_PAD,radiative_DTM_PAD_mean,inventory_PAD):

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
    ax1a = plt.subplot2grid((6,6),(0,0),colspan=2)
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
                axes2[pp].fill_betweenx(heights[2:],0,MacArthurHorn_PAD[Plot_name][i,2:],color=colour[0],alpha=0.01)
            axes2[pp].plot(MacArthurHorn_PAD_mean[Plot_name][2:],heights[2:],'-',c=colour[0],linewidth=2)
            # plot detto profile
            for i in range(0,n_subplots):
                axes3[pp].fill_betweenx(heights_rad[3:],0,radiative_PAD[Plot_name][i,:-3,-1][::-1],color=colour[1],alpha=0.01)
            axes3[pp].plot(radiative_PAD_mean[Plot_name][:-3,1][::-1],heights_rad[3:],'-',c=colour[1],linewidth=2)

            # plot corrective radiative transfer profile
            for i in range(0,n_subplots):
                axes4[pp].fill_betweenx(heights_rad[3:],0,radiative_DTM_PAD[Plot_name][i,:-3,-1][::-1],color=colour[1],alpha=0.01)
            axes4[pp].plot(radiative_DTM_PAD_mean[Plot_name][:-3,1][::-1],heights_rad[3:],'-',c=colour[1],linewidth=2)

            # field inventory
            #for i in range(0,n_subplots):
                #axes5[pp].fill_betweenx(heights[2:],0,inventory_PAD[Plot_name][i,2:],color=colour[2],alpha=0.05)
            #axes5[pp].plot(np.mean(inventory_PAD[Plot_name],axis=0)[2:],heights[2:],'-',c=colour[2],linewidth=2)
            axes5[pp].plot(inventory_PAD[Plot_name][2:],heights[2:],'-',c=colour[2],linewidth=2)

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

def plot_point_clouds_and_profiles(figure_name,figure_number, gps_pts_file,plot_point_cloud,heights,heights_rad,
                                lidar_profiles,MacArthurHorn_PAD,MacArthurHorn_PAD_mean,
                                radiative_DTM_PAD,radiative_DTM_PAD_mean,inventory_PAD):
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
    plot_coords['plot'][plot_coords['plot']==b'danum_1'] = 'DC1'
    plot_coords['plot'][plot_coords['plot']==b'danum_2'] = 'DC2'
    plot_coords['plot'][plot_coords['plot']==b'maliau_belian'] = 'Belian'
    plot_coords['plot'][plot_coords['plot']==b'maliau_seraya'] = 'Seraya'
    plot_coords['plot'][plot_coords['plot']==b'B_north'] = 'B North'
    plot_coords['plot'][plot_coords['plot']==b'B_south'] = 'B South'
    plot_coords['plot'][plot_coords['plot']==b'LFE'] = 'LF'

    affine = {}
    fig1_plots = [b'Belian', b'Seraya',b'LF',b'E', b'B North', b'B South']
    Plots =      [b'Belian', b'Seraya',b'LF',b'E', b'B North', b'B South']
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
    ax1a = plt.subplot2grid((6,6),(0,0),colspan=2)
    ax1a.annotate('a - Old growth, MLA01', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1a.set_ylabel('Height / m',fontsize=axis_size)
    plt.gca().set_aspect('equal', adjustable='box-forced')

    # Seraya
    ax1b = plt.subplot2grid((6,6),(1,0),sharey=ax1a,sharex=ax1a,colspan=2)
    ax1b.annotate('b - Old growth, MLA02', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1b.set_ylabel('Height / m',fontsize=axis_size)
    plt.gca().set_aspect('equal', adjustable='box-forced')

    # LF
    ax1c = plt.subplot2grid((6,6),(2,0),sharey=ax1a,sharex=ax1a,colspan=2)
    ax1c.annotate('c - Moderately logged, SAF04', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1c.set_ylabel('Height / m',fontsize=axis_size)
    plt.gca().set_aspect('equal', adjustable='box-forced')

    # E
    ax1d = plt.subplot2grid((6,6),(3,0),sharey=ax1a,sharex=ax1a,colspan=2)
    ax1d.annotate('d - Moderately logged, SAF05', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1d.set_ylabel('Height / m',fontsize=axis_size)
    plt.gca().set_aspect('equal', adjustable='box-forced')

    # B North
    ax1e = plt.subplot2grid((6,6),(4,0),sharey=ax1a,sharex=ax1a,colspan=2)
    ax1e.annotate('e - Heavily logged, SAF02', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1e.set_ylabel('Height / m',fontsize=axis_size)
    plt.gca().set_aspect('equal', adjustable='box-forced')

    # B South
    ax1f = plt.subplot2grid((6,6),(5,0),sharey=ax1a,sharex=ax1a,colspan=2)
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

    #---------------------------------------------------------
    # NOW PLOT PROFILES
    # Belian
    ax2a = plt.subplot2grid((6,6),(0,2),sharey=ax1a)
    ax2a.annotate('LiDAR returns', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=9)
    # - MacHorn
    ax3a = plt.subplot2grid((6,6),(0,3),sharey=ax1a)
    ax3a.annotate('MacArthur-Horn', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=9)
    # - Detto
    ax4a = plt.subplot2grid((6,6),(0,4),sharey=ax1a,sharex=ax3a)
    ax4a.annotate('multi. return \nrad. trans.', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=9)
    # - Corrected rad trans
    ax5a = plt.subplot2grid((6,6),(0,5),sharey=ax1a,sharex=ax3a)
    ax5a.annotate('crown volume', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=9)

    # Seraya
    ax2b = plt.subplot2grid((6,6),(1,2),sharey=ax1a,sharex=ax2a)
    ax3b = plt.subplot2grid((6,6),(1,3),sharey=ax1a,sharex=ax3a)
    ax4b = plt.subplot2grid((6,6),(1,4),sharey=ax1a,sharex=ax3a)
    ax5b = plt.subplot2grid((6,6),(1,5),sharey=ax1a,sharex=ax3a)

    # LF
    ax2c = plt.subplot2grid((6,6),(2,2),sharey=ax1a,sharex=ax2a)
    ax3c = plt.subplot2grid((6,6),(2,3),sharey=ax1a,sharex=ax3a)
    ax4c = plt.subplot2grid((6,6),(2,4),sharey=ax1a,sharex=ax3a)
    ax5c = plt.subplot2grid((6,6),(2,5),sharey=ax1a,sharex=ax3a)

    # E
    ax2d = plt.subplot2grid((6,6),(3,2),sharey=ax1a,sharex=ax2a)
    ax3d = plt.subplot2grid((6,6),(3,3),sharey=ax1a,sharex=ax3a)
    ax4d = plt.subplot2grid((6,6),(3,4),sharey=ax1a,sharex=ax3a)
    ax5d = plt.subplot2grid((6,6),(3,5),sharey=ax1a,sharex=ax3a)

    # B North
    ax2e = plt.subplot2grid((6,6),(4,2),sharey=ax1a,sharex=ax2a)
    ax3e = plt.subplot2grid((6,6),(4,3),sharey=ax1a,sharex=ax3a)
    ax4e = plt.subplot2grid((6,6),(4,4),sharey=ax1a,sharex=ax3a)
    ax5e = plt.subplot2grid((6,6),(4,5),sharey=ax1a,sharex=ax3a)

    # B South
    ax2f = plt.subplot2grid((6,6),(5,2), sharex = ax1a, sharey = ax2a)
    ax2f.set_xlabel('Number of returns\n(x1000)',fontsize=axis_size,horizontalalignment='center')
    # - MacHorn
    ax3f = plt.subplot2grid((6,6),(5,3),sharey=ax1a, sharex=ax3a)
    ax3f.set_xlabel('PAD\n(m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size,horizontalalignment='center')
    # - Corrected rad trans
    ax4f = plt.subplot2grid((6,6),(5,4),sharey=ax1a,sharex=ax3a)
    ax4f.set_xlabel('PAD\n(m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size,horizontalalignment='center')
    # - Inventory
    ax5f = plt.subplot2grid((6,6),(5,5),sharey=ax1a,sharex=ax3a)
    ax5f.set_xlabel('Crown Volume\n(m$^3$m$^{-2}$m$^{-1}$)',fontsize=axis_size,horizontalalignment='center')

    axes1 = [ax2a,  ax2b, ax2c, ax2d, ax2e, ax2f]
    axes2 = [ax3a,  ax3b, ax3c, ax3d, ax3e, ax3f]
    axes3 =  [ax4a,  ax4b, ax4c, ax4d, ax4e, ax4f]
    axes4 = [ax5a,  ax5b, ax5c, ax5d, ax5e, ax5f]

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
                axes2[pp].fill_betweenx(heights[2:],0,MacArthurHorn_PAD[Plot_name][i,2:],color=colour[0],alpha=0.01)
            axes2[pp].plot(MacArthurHorn_PAD_mean[Plot_name][2:],heights[2:],'-',c=colour[0],linewidth=2)

            # plot corrective radiative transfer profile
            for i in range(0,n_subplots):
                axes3[pp].fill_betweenx(heights_rad[3:],0,radiative_DTM_PAD[Plot_name][i,:-3,-1][::-1],color=colour[1],alpha=0.01)
            axes3[pp].plot(radiative_DTM_PAD_mean[Plot_name][:-3,1][::-1],heights_rad[3:],'-',c=colour[1],linewidth=2)

            # field inventory
            #for i in range(0,n_subplots):
                #axes5[pp].fill_betweenx(heights[2:],0,inventory_PAD[Plot_name][i,2:],color=colour[2],alpha=0.05)
            #axes5[pp].plot(np.mean(inventory_PAD[Plot_name],axis=0)[2:],heights[2:],'-',c=colour[2],linewidth=2)
            axes4[pp].plot(inventory_PAD[Plot_name][2:],heights[2:],'-',c=colour[2],linewidth=2)

        yticklabels.append(axes1[pp].get_yticklabels())
        yticklabels.append(axes2[pp].get_yticklabels())
        yticklabels.append(axes3[pp].get_yticklabels())
        yticklabels.append(axes4[pp].get_yticklabels())

        if pp < 5:
            xticklabels.append(axes1[pp].get_xticklabels())
            xticklabels.append(axes2[pp].get_xticklabels())
            xticklabels.append(axes3[pp].get_xticklabels())
            xticklabels.append(axes4[pp].get_xticklabels())

    ax1a.set_xlim(0,100)
    ax1a.set_ylim(0,80)
    ax2a.set_xlim(0,29)
    ax3a.set_xlim(xmin=0,xmax=0.7)

    ax2f.locator_params(axis='x',nbins=5)
    ax3f.locator_params(axis='x',nbins=5)
    ax4f.locator_params(axis='x',nbins=5)
    ax5f.locator_params(axis='x',nbins=5)

    plt.setp(yticklabels,visible=False)
    plt.setp(xticklabels,visible=False)
    plt.subplots_adjust(hspace=0.2, wspace = 0.1)

    plt.tight_layout()
    plt.savefig(figure_name)
    plt.show()
    return 0

def plot_point_clouds_and_profiles_Danum(figure_name,figure_number, gps_pts_file,plot_point_cloud,heights,heights_rad,
                                lidar_profiles,MacArthurHorn_PAD,MacArthurHorn_PAD_mean,
                                radiative_DTM_PAD,radiative_DTM_PAD_mean,inventory_PAD):
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
    plot_coords['plot'][plot_coords['plot']==b'danum_1'] = 'DC1'
    plot_coords['plot'][plot_coords['plot']==b'danum_2'] = 'DC2'

    affine = {}
    fig1_plots = [b'DC1', b'DC2']
    Plots =      [b'DC1', b'DC2']
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

    plt.figure(figure_number, facecolor='White',figsize=[9,4])

    # DC1
    ax1a = plt.subplot2grid((2,6),(0,0),colspan=2)
    ax1a.annotate('a - Old growth, DAN04', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1a.set_ylabel('Height / m',fontsize=axis_size)
    plt.gca().set_aspect('equal', adjustable='box-forced')

    # DC2
    ax1b = plt.subplot2grid((2,6),(1,0),sharey=ax1a,sharex=ax1a,colspan=2)
    ax1b.annotate('b - Old growth, DAN05', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1b.set_ylabel('Height / m',fontsize=axis_size)
    plt.gca().set_aspect('equal', adjustable='box-forced')

    axes = [ax1a, ax1b]
    for pp in range(0,2):
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

    #---------------------------------------------------------
    # NOW PLOT PROFILES
    # DC1
    ax2a = plt.subplot2grid((2,6),(0,2),sharey=ax1a)
    ax2a.annotate('LiDAR returns', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=9)
    # - MacHorn
    ax3a = plt.subplot2grid((2,6),(0,3),sharey=ax1a)
    ax3a.annotate('MacArthur-Horn', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=9)
    # - Detto
    ax4a = plt.subplot2grid((2,6),(0,4),sharey=ax1a,sharex=ax3a)
    ax4a.annotate('multi. return \nrad. trans.', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=9)
    # - Corrected rad trans
    ax5a = plt.subplot2grid((2,6),(0,5),sharey=ax1a,sharex=ax3a)
    ax5a.annotate('crown volume', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=9)

    # Danum 2
    ax2b = plt.subplot2grid((2,6),(1,2), sharex = ax2a, sharey = ax1a)
    ax2b.set_xlabel('Number of returns',fontsize=axis_size,horizontalalignment='center')
    # - MacHorn
    ax3b = plt.subplot2grid((2,6),(1,3),sharey=ax1a, sharex=ax3a)
    ax3b.set_xlabel('PAD\n(m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size,horizontalalignment='center')
    # - Corrected rad trans
    ax4b = plt.subplot2grid((2,6),(1,4),sharey=ax1a,sharex=ax3a)
    ax4b.set_xlabel('PAD\n(m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size,horizontalalignment='center')
    # - Inventory
    ax5b = plt.subplot2grid((2,6),(1,5),sharey=ax1a,sharex=ax3a)
    ax5b.set_xlabel('Crown Volume\n(m$^3$m$^{-2}$m$^{-1}$)',fontsize=axis_size,horizontalalignment='center')

    axes1 = [ax2a,  ax2b]
    axes2 = [ax3a,  ax3b]
    axes3 =  [ax4a,  ax4b]
    axes4 = [ax5a,  ax5b]

    yticklabels=[]
    xticklabels=[]

    for pp in range(0,2):
        Plot_name = fig1_plots[pp]

        # plot lidar profile
        return_dist     = np.sum(lidar_profiles[Plot_name],axis=0)
        for k in range(0,max_return):
            axes1[pp].plot(return_dist[:,k],np.max(heights_rad)-heights_rad,'-',c=colour[k],linewidth=1)

            # plot macarthur horn profile
            for i in range(0,n_subplots):
                axes2[pp].fill_betweenx(heights[2:],0,MacArthurHorn_PAD[Plot_name][i,2:],color=colour[0],alpha=0.01)
            axes2[pp].plot(MacArthurHorn_PAD_mean[Plot_name][2:],heights[2:],'-',c=colour[0],linewidth=2)

            # plot corrective radiative transfer profile
            for i in range(0,n_subplots):
                axes3[pp].fill_betweenx(heights_rad[3:],0,radiative_DTM_PAD[Plot_name][i,:-3,-1][::-1],color=colour[1],alpha=0.01)
            axes3[pp].plot(radiative_DTM_PAD_mean[Plot_name][:-3,1][::-1],heights_rad[3:],'-',c=colour[1],linewidth=2)

            # field inventory
            axes4[pp].plot(inventory_PAD[Plot_name][2:],heights[2:],'-',c=colour[2],linewidth=2)

        yticklabels.append(axes1[pp].get_yticklabels())
        yticklabels.append(axes2[pp].get_yticklabels())
        yticklabels.append(axes3[pp].get_yticklabels())
        yticklabels.append(axes4[pp].get_yticklabels())

        if pp == 0:
            xticklabels.append(axes1[pp].get_xticklabels())
            xticklabels.append(axes2[pp].get_xticklabels())
            xticklabels.append(axes3[pp].get_xticklabels())
            xticklabels.append(axes4[pp].get_xticklabels())

    ax1a.set_xlim(0,100)
    ax1a.set_ylim(0,80)
    ax2a.set_xlim(0,3500)
    ax3a.set_xlim(xmin=0,xmax=0.7)

    ax2b.locator_params(axis='x',nbins=4)
    ax3b.locator_params(axis='x',nbins=5)
    ax4b.locator_params(axis='x',nbins=5)
    ax5b.locator_params(axis='x',nbins=5)

    plt.setp(yticklabels,visible=False)
    plt.setp(xticklabels,visible=False)
    plt.subplots_adjust(hspace=0.2, wspace = 0.1)

    plt.tight_layout()
    plt.savefig(figure_name)
    plt.show()
    return 0

#=======================================================================================
"""
# Compare LiDAR approaches
"""
def compare_LiDAR_PAI(figure_name,figure_number,MacArthurHorn_PAD,MacArthurHorn_PAD_mean,
                        radiative_PAD,radiative_PAD_mean,
                        radiative_DTM_PAD,radiative_DTM_PAD_mean,layer_thickness=1):
    Plots = MacArthurHorn_PAD.keys()
    N_plots = len(Plots)
    N_subplots = np.shape(MacArthurHorn_PAD[Plots[0]])[0]
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
        print(Plots[pp])
        mh_20m[pp,:]=np.nansum(MacArthurHorn_PAD[Plots[pp]],axis=1)*layer_thickness
        rad2_20m[pp,:]=np.nansum(radiative_PAD[Plots[pp]][:,:,1],axis=1)*layer_thickness
        rad3_20m[pp,:]=np.nansum(radiative_PAD[Plots[pp]][:,:,2],axis=1)*layer_thickness
        rad2_DTM_20m[pp,:]=np.nansum(radiative_DTM_PAD[Plots[pp]][:,:,1],axis=1)*layer_thickness
        rad3_DTM_20m[pp,:]=np.nansum(radiative_DTM_PAD[Plots[pp]][:,:,2],axis=1)*layer_thickness
        print(rad2_DTM_20m[pp].astype('int'))
        print(mh_20m[pp])

        mh_1ha[pp] = np.nansum(MacArthurHorn_PAD_mean[Plots[pp]])*layer_thickness
        rad2_1ha[pp] = np.nansum(radiative_PAD_mean[Plots[pp]][:,1])*layer_thickness
        rad3_1ha[pp] = np.nansum(radiative_PAD_mean[Plots[pp]][:,2])*layer_thickness
        rad2_DTM_1ha[pp] = np.nansum(radiative_DTM_PAD_mean[Plots[pp]][:,1])*layer_thickness
        rad3_DTM_1ha[pp] = np.nansum(radiative_DTM_PAD_mean[Plots[pp]][:,2])*layer_thickness
        print(rad2_DTM_1ha[pp])
        print(mh_1ha[pp])

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
"""
# LAI vs. canopy volume
"""
def plot_LAI_vs_inventory(figure_name,figure_number,MacArthurHorn_PAD,MacArthurHorn_PAD_mean,
                            radiative_PAD,radiative_PAD_mean,radiative_DTM_PAD,
                            radiative_DTM_PAD_mean,Inventory_PAD,Inventory_LAI,layer_thickness=1):
    Plots = MacArthurHorn_PAD.keys()
    N_plots = len(Plots)
    N_subplots = np.shape(MacArthurHorn_PAD[Plots[0]])[0]
    vol = np.zeros(N_plots)
    mh = np.zeros(N_plots)
    rad2 = np.zeros(N_plots)
    rad3 = np.zeros(N_plots)
    radDTM2 = np.zeros(N_plots)
    radDTM3 = np.zeros(N_plots)

    for pp in range(0,N_plots):
        print(Plots[pp])
        # for 1 ha averages, want to avoid nodata
        # these should be identifiable based on penetration limit
        vol[pp] = np.sum(Inventory_PAD[Plots[pp]])
        mh[pp] = np.nansum(MacArthurHorn_PAD_mean[Plots[pp]])*layer_thickness
        rad2[pp] = np.nansum(radiative_PAD_mean[Plots[pp]][:,1])*layer_thickness
        rad3[pp] = np.nansum(radiative_PAD_mean[Plots[pp]][:,2])*layer_thickness
        radDTM2[pp] = np.nansum(radiative_DTM_PAD_mean[Plots[pp]][:,1])*layer_thickness
        radDTM3[pp] = np.nansum(radiative_DTM_PAD_mean[Plots[pp]][:,2])*layer_thickness

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
        ax1.errorbar(np.mean(Inventory_LAI[Plots[i]]),np.nansum(MacArthurHorn_PAD_mean[Plots[i]])*layer_thickness,xerr=x_err,yerr=y_err,marker='o',color='black')


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
                ax2.errorbar(np.mean(Inventory_LAI[Plots[i]]), np.nansum(radiative_PAD_mean[Plots[i]][:,k])*layer_thickness,xerr=x_err,yerr=y_err,marker='o',color=colour[k],label=leg_label)
            else:
                #ax2.errorbar(np.mean(inventory_LAI[Plots[i]]),np.mean(radiative_LAI[Plots[i]][:,k]),xerr=x_err,yerr=y_err,marker='o',color=colour[k])
                ax2.errorbar(np.mean(Inventory_LAI[Plots[i]]), np.nansum(radiative_PAD_mean[Plots[i]][:,k])*layer_thickness,xerr=x_err,yerr=y_err,marker='o',color=colour[k])

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
                ax3.errorbar(np.mean(Inventory_LAI[Plots[i]]), np.nansum(radiative_DTM_PAD_mean[Plots[i]][:,k])*layer_thickness,xerr=x_err,yerr=y_err,marker='o',color=colour[k],label=leg_label)
            else:
                #ax3.errorbar(np.mean(inventory_LAI[Plots[i]]),np.mean(radiative_DTM_LAI[Plots[i]][:,k]),xerr=x_err,yerr=y_err,marker='o',color=colour[k])
                ax3.errorbar(np.mean(Inventory_LAI[Plots[i]]), np.nansum(radiative_DTM_PAD_mean[Plots[i]][:,k])*layer_thickness,xerr=x_err,yerr=y_err,marker='o',color=colour[k])

    ax3.legend(loc=4)

    #ax1.set_ylim((0,20))
    plt.tight_layout()
    plt.savefig(figure_name)
    plt.show()
    return 0

"""
# LAI vs. basal area
"""
def plot_LAI_vs_basal_area(figure_name,figure_number,MacArthurHorn_PAD,MacArthurHorn_PAD_mean,
                            radiative_DTM_PAD,radiative_DTM_PAD_mean,BasalArea,
                            plot_marker,plot_label,plot_colour,layer_thickness=1):


    Plots = list(MacArthurHorn_PAD.keys())
    N_plots = len(Plots)
    N_subplots = np.shape(MacArthurHorn_PAD[Plots[0]])[0]
    BA = np.zeros(N_plots)
    mh = np.zeros(N_plots)
    rad2 = np.zeros(N_plots)
    rad3 = np.zeros(N_plots)
    radDTM2 = np.zeros(N_plots)
    radDTM3 = np.zeros(N_plots)

    for pp in range(0,N_plots):
        print(Plots[pp])
        Plots_str = str(Plots[pp], 'utf-8')
        # for 1 ha averages, want to avoid nodata
        # these should be identifiable based on penetration limit
        BA[pp] = np.mean(BasalArea[Plots_str]).copy()
        mh[pp] = np.nansum(MacArthurHorn_PAD_mean[Plots[pp]])*layer_thickness
        radDTM2[pp] = np.nansum(radiative_DTM_PAD_mean[Plots[pp]][:,1])*layer_thickness
        radDTM3[pp] = np.nansum(radiative_DTM_PAD_mean[Plots[pp]][:,2])*layer_thickness

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
    m, c, r2, p, x_i, y_i, CI_u, CI_l, PI_u, PI_l = field.linear_regression(BA,mh,conf=0.95)
    ax1.fill_between(x_i,CI_l,CI_u,color='0.75')
    ax1.plot(x_i,y_i,'-',color='black')
    for i in range(0,N_plots):
        Plot_str = str(Plots[i], 'utf-8')
        ax1.plot(BA[i],mh[i],marker=plot_marker[Plot_str],
                    color=plot_colour[Plot_str],linestyle='None')

    ax3 = plt.subplot2grid((1,2),(0,1), sharex=ax1)#, sharey=ax1)
    ax3.annotate('b - multi. return\nradiative transfer', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=axis_size)
    ax3.annotate(r_sq_b, xy=(0.5,0.05), xycoords='axes fraction',backgroundcolor='none',
            horizontalalignment='center', verticalalignment='bottom', fontsize=axis_size)
    ax3.set_ylabel('PAI',fontsize=axis_size)
    ax3.set_xlabel('Basal Area / m$^2$ha$^{-1}$',fontsize=axis_size)

    m, c, r2, p, x_i, y_i, CI_u, CI_l, PI_u, PI_l = field.linear_regression(BA,radDTM3,conf=0.95)
    ax3.fill_between(x_i,CI_l,CI_u,color='0.75')
    ax3.plot(x_i,y_i,'-',color='black')
    for i in range(0,N_plots):
        Plot_str = str(Plots[i], 'utf-8')
        ax3.plot(BA[i],radDTM3[i],marker=plot_marker[Plot_str],
                    color=plot_colour[Plot_str],linestyle='None',label=plot_label[Plot_str])

    ax3.legend(loc=4)

    plt.tight_layout()
    plt.savefig(figure_name)
    plt.show()
    return 0

"""
# Location map
"""
def plot_location_map(figure_name,figure_number):

    AGB_file = '/home/dmilodow/DataStore_GCEL/AGB/avitabile/Avitabile_AGB_Map/Avitabile_AGB_Map.tif'
    ds = gdal.Open(AGB_file)
    geoTrans = ds.GetGeoTransform()
    W = 114.; E = 120; N = 9; S = 3
    agb,geoT=geo.clip_array_to_bbox(np.asarray(ds.ReadAsArray(),dtype='float'),geoTrans,N,S,W,E)
    agb[agb<0]=np.nan; array = ma.masked_invalid(agb)

    xres = geoT[1]; yres = geoT[5]; x0 = geoT[0]; y0 = geoT[3]
    rows,cols = array.shape
    lons = np.arange(x0,x0+xres*cols,xres)+ xres * 0.5
    lats = np.arange(y0,y0+yres*rows,yres)+ yres * 0.5
    crs = ccrs.PlateCarree()

    fig = plt.figure(figure_number, facecolor='White',figsize=[10,7])
    # Region
    ax = plt.axes(projection=crs)
    ax.background_patch.set_facecolor('0.8')
    im = plt.pcolormesh(lons,lats, array, cmap='Greens')

    country_shp = '/home/dmilodow/DataStore_DTM/EOlaboratory/Areas/NaturalEarth/10m_cultural/ne_10m_admin_0_boundary_lines_land.shp'
    shp = Reader(country_shp)
    Boundaries_left = [i for i in shp.records()
                    if i.attributes['adm0_left']=='Malaysia']
    Boundaries_right = [i for i in shp.records()
                    if i.attributes['adm0_right']=='Malaysia']
    for b1 in Boundaries_right:
        if b1.attributes['adm0_left'] != 'Singapore':
            sp = ShapelyFeature(b1.geometry, crs,  facecolor = 'none', edgecolor='black')
            ax.add_feature(sp)
    for b2 in Boundaries_left:
        sp = ShapelyFeature(b2.geometry, crs,  facecolor = 'none', edgecolor='black')
        ax.add_feature(sp)

    states_shp = '/home/dmilodow/DataStore_DTM/EOlaboratory/Areas/NaturalEarth/10m_cultural/ne_10m_admin_1_states_provinces_lines_shp.shp'
    shp = Reader(states_shp)
    Boundaries = [i for i in shp.records()
                    if i.attributes['adm0_name']=='Malaysia']
    for b in Boundaries:
        ax.plot(b.geometry.coords.xy[0],b.geometry.coords.xy[1], ":", color='0.1',lw=1)

    # Plot locations
    plot_shp = 'plot_locations_for_location_map_WGS84.shp'
    pshp = Reader(plot_shp)
    for pp in pshp.records():
        col=colour[1]
        if pp.attributes['ForestType'] == 'OG':
            col = colour[0]
        elif pp.attributes['ForestType'] == 'ML':
            col=colour[2]
        ax.plot(pp.geometry.coords.xy[0],pp.geometry.coords.xy[1], marker='o', mfc=col,mec='white', mew=0.5, markersize=6)

    ax.set_xlim(left=W,right=E);ax.set_ylim(top=N,bottom=S)

    # Grid lines
    parallels = np.arange(0.,12.,2.); meridians = np.arange(110.,122.,2.)

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', linestyle='--')
    gl.xlabels_bottom = False; gl.ylabels_left = False
    gl.xlocator = mticker.FixedLocator(meridians); gl.ylocator = mticker.FixedLocator(parallels)
    gl.xformatter = LONGITUDE_FORMATTER; gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'color': '0.5'}; gl.ylabel_style = {'color': '0.5'}

    # colorbar bits and bobs
    divider = make_axes_locatable(ax)
    cax = divider.new_horizontal(size="5%", pad=0.8, axes_class=plt.Axes)
    fig.add_axes(cax)
    cbar=plt.colorbar(im, cax=cax)#,extend = 'min')
    cbar.ax.set_ylabel('AGB / Mg ha$^{-1}$',fontsize = 15)
    cbar.solids.set_edgecolor("face")

    # Labels
    ax.annotate('SABAH', xy=(117,5.7), xycoords='data',
            backgroundcolor='none',horizontalalignment='center', verticalalignment='center',
            fontsize=12, style='italic')
    ax.annotate('SARAWAK', xy=(114.7,3.56), xycoords='data',
            backgroundcolor='none',horizontalalignment='center', verticalalignment='center',
            fontsize=12, style='italic')
    ax.annotate('BRUNEI', xy=(114.4,4.74), xycoords='data',
            backgroundcolor='none',horizontalalignment='center', verticalalignment='center',
            fontsize=12, style='italic')
    ax.annotate('N. KALIMANTAN', xy=(116.5,3.67), xycoords='data',
            backgroundcolor='none',horizontalalignment='center', verticalalignment='center',
            fontsize=12, style='italic')
    ax.annotate('SAFE', xy=(117.7,4.6), xycoords='data',backgroundcolor='none',
            horizontalalignment='left', verticalalignment='top', fontsize=12,
            color='black', fontweight='bold')
    ax.annotate('Maliau\nBasin', xy=(116.9,4.8), xycoords='data',
        backgroundcolor='none',horizontalalignment='right', verticalalignment='bottom',
        fontsize=12,color='black', fontweight='bold')
    ax.annotate('Danum Valley', xy=(117.8,5.), xycoords='data',
        backgroundcolor='none',horizontalalignment='left', verticalalignment='bottom',
        fontsize=12,color='black', fontweight='bold')

    plt.savefig(figure_name)
    plt.show()
    return 0


# Cross-plot canopy layers (new)
def cross_plot_canopy_layers_LiDAR(figure_name,figure_number,heights,heights_rad,MacArthurHorn_PAD,MacArthurHorn_PAD_mean,radiative_PAD,radiative_PAD_mean,radiative_DTM_PAD,radiative_DTM_PAD_mean,inventory_PAD):

    # Process the profiles

    # 2 columns (plotting MacHorn on x axis, multi return on y axis)
    # 3 rows (old growth, moderately logged, heavily logged)
    fig = plt.figure(figure_number, facecolor='White',figsize=[6,9])

    # First up - old-growth forest
    Plot_name='Belian'
    x=np.cumsum(MacArthurHorn_PAD_mean[Plot_name][2:][::-1])
    y1=np.cumsum(radiative_PAD_mean[Plot_name][:-3,1])
    y2=np.cumsum(radiative_DTM_PAD_mean[Plot_name][:-3,1])
    t=heights[2:][::-1]
    ax1 = plt.subplot2grid((3,2),(0,0))
    ax1.annotate('a', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)
    for i in range(0,n_subplots):
        ax1.plot(np.cumsum(MacArthurHorn_PAD[Plot_name][i,2:][::-1]),np.cumsum(radiative_PAD[Plot_name][i,:-3,1]),'-',color='k',linewidth=0.5)
    plot_colourline(ax1,x,y1,t,linewidth=5,cmap='plasma')

    ax2 = plt.subplot2grid((3,2),(0,1),sharex=ax1,sharey=ax1)
    ax2.annotate('b', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)
    for i in range(0,n_subplots):
        ax2.plot(np.cumsum(MacArthurHorn_PAD[Plot_name][i,2:][::-1]),np.cumsum(radiative_DTM_PAD[Plot_name][i,:-3,1]),'-',color='k',linewidth=0.5)
    plot_colourline(ax2,x,y2,t,linewidth=5,cmap='plasma')

    # Next up - Moderately logged forest
    Plot_name='E'
    x=np.cumsum(MacArthurHorn_PAD_mean[Plot_name][2:][::-1])
    y1=np.cumsum(radiative_PAD_mean[Plot_name][:-3,1])
    y2=np.cumsum(radiative_DTM_PAD_mean[Plot_name][:-3,1])
    ax3 = plt.subplot2grid((3,2),(1,0),sharex=ax1,sharey=ax1)
    ax3.annotate('c', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)
    ax3.set_ylabel('cumulative PAD$_{rad trans (Detto)}$ (m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size)

    for i in range(0,n_subplots):
        ax3.plot(np.cumsum(MacArthurHorn_PAD[Plot_name][i,2:][::-1]),np.cumsum(radiative_PAD[Plot_name][i,:-3,1]),'-',color='k',linewidth=0.5)
    plot_colourline(ax3,x,y1,t,linewidth=5,cmap='plasma')

    ax4 = plt.subplot2grid((3,2),(1,1),sharex=ax1,sharey=ax1)
    ax4.annotate('d', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)
    ax4.set_ylabel('cumulative PAD$_{rad trans (new)}$ (m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size)
    ax4.yaxis.set_label_position("right")
    ax4.yaxis.tick_right()
    ax4.yaxis.set_ticks_position('both')

    for i in range(0,n_subplots):
        ax4.plot(np.cumsum(MacArthurHorn_PAD[Plot_name][i,2:][::-1]),np.cumsum(radiative_DTM_PAD[Plot_name][i,:-3,1]),'-',color='k',linewidth=0.5)

    plot_colourline(ax4,x,y2,t,linewidth=5,cmap='plasma')

    # Finally heavily logged forest
    Plot_name = 'B North'
    x=np.cumsum(MacArthurHorn_PAD_mean[Plot_name][2:][::-1])
    y1=np.cumsum(radiative_PAD_mean[Plot_name][:-3,1])
    y2=np.cumsum(radiative_DTM_PAD_mean[Plot_name][:-3,1])
    ax5 = plt.subplot2grid((3,2),(2,0),sharex=ax1,sharey=ax1)
    ax5.annotate('e', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)

    for i in range(0,n_subplots):
        ax5.plot(np.cumsum(MacArthurHorn_PAD[Plot_name][i,2:][::-1]),np.cumsum(radiative_PAD[Plot_name][i,:-3,1]),'-',color='k',linewidth=0.5)
    plot_colourline(ax5,x,y1,t,linewidth=5,cmap='plasma')

    ax6 = plt.subplot2grid((3,2),(2,1),sharex=ax1,sharey=ax1)
    ax6.annotate('f', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)

    for i in range(0,n_subplots):
        ax6.plot(np.cumsum(MacArthurHorn_PAD[Plot_name][i,2:][::-1]),np.cumsum(radiative_DTM_PAD[Plot_name][i,:-3,1]),'-',color='k',linewidth=0.5)

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

"""
# Plot canopy residuals for comparison between the Mac-Horn and multi-return
# radiative transfer approachself.
# Old version also includes the uncorrected Detto model
"""
def plot_canopy_layer_residuals_old(figure_name,figure_number,heights,MacArthurHorn_PAD,MacArthurHorn_PAD_mean,radiative_PAD,radiative_PAD_mean,radiative_DTM_PAD,radiative_DTM_PAD_mean,max_return=2):

    # 2 columns (plotting MacHorn on x axis, multi return on y axis)
    # 3 rows (old growth, moderately logged, heavily logged)
    fig = plt.figure(figure_number, facecolor='White',figsize=[6,9])

    # First up - old-growth forest
    Plot_name='Belian'
    ids = np.arange(heights.size)
    idmax = ids[MacArthurHorn_PAD_mean[Plot_name]>0][-1]+1
    y=heights[2:idmax]
    x=MacArthurHorn_PAD_mean[Plot_name][2:idmax]
    x1=radiative_PAD_mean[Plot_name][:-1,max_return-1][::-1][2:idmax]
    x2=radiative_DTM_PAD_mean[Plot_name][:-1,max_return-1][::-1][2:idmax]

    ax1 = plt.subplot2grid((3,2),(0,0))
    ax1.annotate('a', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)
    for i in range(0,n_subplots):
        ax1.plot((radiative_PAD[Plot_name][i,:-1,max_return-1][::-1][2:idmax]-MacArthurHorn_PAD[Plot_name][i,2:idmax]),y,'-',color='0.5',linewidth=0.5,alpha=0.5)
    ax1.plot(x1-x,y,linewidth=2,color='k')

    ax2 = plt.subplot2grid((3,2),(0,1),sharex=ax1,sharey=ax1)
    ax2.annotate('b', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)
    ax2.yaxis.tick_right()
    ax2.yaxis.set_ticks_position('both')
    for i in range(0,n_subplots):
        ax2.plot((radiative_DTM_PAD[Plot_name][i,:-1,max_return-1][::-1][2:idmax]-MacArthurHorn_PAD[Plot_name][i,2:idmax]),y,'-',color='0.5',linewidth=0.5,alpha=0.5)
    ax2.plot(x2-x,y,linewidth=2,color='k')

    # Next up - Moderately logged forest
    Plot_name='E'
    ids = np.arange(heights.size)
    idmax = ids[MacArthurHorn_PAD_mean[Plot_name]>0][-1]+1
    y=heights[2:idmax]
    x=MacArthurHorn_PAD_mean[Plot_name][2:idmax]
    x1=radiative_PAD_mean[Plot_name][:-1,max_return-1][::-1][2:idmax]
    x2=radiative_DTM_PAD_mean[Plot_name][:-1,max_return-1][::-1][2:idmax]
    ax3 = plt.subplot2grid((3,2),(1,0),sharex=ax1,sharey=ax1)
    ax3.annotate('c', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)
    ax3.set_ylabel('height (m)',fontsize=axis_size)

    for i in range(0,n_subplots):
        ax3.plot((radiative_PAD[Plot_name][i,:-1,max_return-1][::-1][2:idmax]-MacArthurHorn_PAD[Plot_name][i,2:idmax]),y,'-',color='0.5',linewidth=0.5,alpha=0.5)
    ax3.plot(x1-x,y,linewidth=2,color='k')

    ax4 = plt.subplot2grid((3,2),(1,1),sharex=ax1,sharey=ax1)
    ax4.annotate('d', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)
    ax4.set_ylabel('height (m))',fontsize=axis_size)
    ax4.yaxis.set_label_position("right")
    ax4.yaxis.tick_right()
    ax4.yaxis.set_ticks_position('both')

    for i in range(0,n_subplots):
        ax4.plot((radiative_DTM_PAD[Plot_name][i,:-1,max_return-1][::-1][2:idmax]-MacArthurHorn_PAD[Plot_name][i,2:idmax]),y,'-',color='0.5',linewidth=0.5,alpha=0.5)
    ax4.plot(x2-x,y,linewidth=2,color='k')

    # Finally heavily logged forest
    Plot_name = 'B South'
    ids = np.arange(heights.size)
    idmax = ids[MacArthurHorn_PAD_mean[Plot_name]>0][-1]+1
    y=heights[2:idmax]
    x=MacArthurHorn_PAD_mean[Plot_name][2:idmax]
    x1=radiative_PAD_mean[Plot_name][:-1,max_return-1][::-1][2:idmax]
    x2=radiative_DTM_PAD_mean[Plot_name][:-1,max_return-1][::-1][2:idmax]
    ax5 = plt.subplot2grid((3,2),(2,0),sharex=ax1,sharey=ax1)
    ax5.annotate('e', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)
    ax5.set_xlabel('PAD$_{rad trans-Detto}$-PAD$_{MacArthur-Horn}$\n(m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size)

    for i in range(0,n_subplots):
        ax5.plot((radiative_PAD[Plot_name][i,:-1,max_return-1][::-1][2:idmax]-MacArthurHorn_PAD[Plot_name][i,2:idmax]),y,'-',color='0.5',linewidth=0.5,alpha=0.5)
    ax5.plot(x1-x,y,linewidth=2,color='k')

    ax6 = plt.subplot2grid((3,2),(2,1),sharex=ax1,sharey=ax1)
    ax6.annotate('f', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=12)
    ax6.set_xlabel('PAD$_{rad trans-new}$-PAD$_{MacArthur-Horn}$\n(m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size)
    ax6.yaxis.tick_right()
    ax6.yaxis.set_ticks_position('both')

    for i in range(0,n_subplots):
        ax6.plot((radiative_DTM_PAD[Plot_name][i,:-1,max_return-1][::-1][2:idmax]-MacArthurHorn_PAD[Plot_name][i,2:idmax]),y,'-',color='0.5',linewidth=0.5,alpha=0.5)
    ax6.plot(x2-x,y,linewidth=2,color='k')

    # tidy up tick labels
    xticklabels=[ax1.get_xticklabels() + ax2.get_xticklabels() + ax3.get_xticklabels() + ax4.get_xticklabels()]
    plt.setp(xticklabels,visible=False)

    plt.tight_layout()
    plt.show()

def plot_canopy_layer_residuals(figure_name,figure_number,heights,MacArthurHorn_PAD,MacArthurHorn_PAD_mean,
                                radiative_DTM_PAD,radiative_DTM_PAD_mean,max_return=2):

    # 3 columns (old growth, moderately logged, heavily logged)
    fig = plt.figure(figure_number, facecolor='White',figsize=[8,4])

    # First up - old-growth forest
    Plot_name=b'Belian'
    ids = np.arange(heights.size)
    idmax = ids[MacArthurHorn_PAD_mean[Plot_name]>0][-1]+1
    y=heights[2:idmax]
    x=MacArthurHorn_PAD_mean[Plot_name][2:idmax]
    x1=radiative_DTM_PAD_mean[Plot_name][:-1,max_return-1][::-1][2:idmax]

    ax1 = plt.subplot2grid((1,3),(0,0))
    ax1.annotate('a - MLA-01', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',
                    horizontalalignment='left', verticalalignment='top', fontsize=12)
    for i in range(0,n_subplots):
        ax1.plot((radiative_DTM_PAD[Plot_name][i,:-1,max_return-1][::-1][2:idmax]-MacArthurHorn_PAD[Plot_name][i,2:idmax]),
                    y,'-',color='0.5',linewidth=0.5,alpha=0.5)
    ax1.plot(x1-x,y,linewidth=2,color='k')

    # Next up - Moderately logged forest
    Plot_name=b'E'
    ids = np.arange(heights.size)
    idmax = ids[MacArthurHorn_PAD_mean[Plot_name]>0][-1]+1
    y=heights[2:idmax]
    x=MacArthurHorn_PAD_mean[Plot_name][2:idmax]
    x1=radiative_DTM_PAD_mean[Plot_name][:-1,max_return-1][::-1][2:idmax]
    ax2 = plt.subplot2grid((1,3),(0,1),sharex=ax1,sharey=ax1)
    ax2.annotate('b - SAF-03', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',
                    horizontalalignment='left', verticalalignment='top', fontsize=12)
    for i in range(0,n_subplots):
        ax2.plot((radiative_DTM_PAD[Plot_name][i,:-1,max_return-1][::-1][2:idmax]-MacArthurHorn_PAD[Plot_name][i,2:idmax]),
                    y,'-',color='0.5',linewidth=0.5,alpha=0.5)
    ax2.plot(x1-x,y,linewidth=2,color='k')
    ax2.set_xlabel('PAD$_{rad trans}$ - PAD$_{MacArthur-Horn}$ / (m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size)


    # Finally heavily logged forest
    Plot_name = b'B South'
    ids = np.arange(heights.size)
    idmax = ids[MacArthurHorn_PAD_mean[Plot_name]>0][-1]+1
    y=heights[2:idmax]
    x=MacArthurHorn_PAD_mean[Plot_name][2:idmax]
    x1=radiative_DTM_PAD_mean[Plot_name][:-1,max_return-1][::-1][2:idmax]
    ax3 = plt.subplot2grid((1,3),(0,2),sharex=ax1,sharey=ax1)
    ax3.annotate('c - SAF-01', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',
                horizontalalignment='left', verticalalignment='top', fontsize=12)
    for i in range(0,n_subplots):
        ax3.plot((radiative_DTM_PAD[Plot_name][i,:-1,max_return-1][::-1][2:idmax]-MacArthurHorn_PAD[Plot_name][i,2:idmax]),
                y,'-',color='0.5',linewidth=0.5,alpha=0.5)
    ax3.plot(x1-x,y,linewidth=2,color='k')
    # tidy up tick labels
    ax3.set_ylabel('height (m))',fontsize=axis_size)
    ax3.yaxis.set_label_position("right")
    ax3.yaxis.tick_right()
    ax3.yaxis.set_ticks_position('both')
    yticklabels=[ax2.get_yticklabels()]
    plt.setp(yticklabels,visible=False)
    ax1.yaxis.set_ticks_position('both')
    ax2.yaxis.set_ticks_position('both')
    ax1.set_xlim(-2,2)

    plt.tight_layout()
    plt.savefig(figure_name)
    plt.show()

    return 0

"""
# plot PAD-volume ratios
"""
def plot_canopy_layer_PAD_volume_ratio(figure_name,figure_number,heights,MacArthurHorn_PAD,MacArthurHorn_PAD_mean,radiative_PAD,radiative_PAD_mean,radiative_DTM_PAD,radiative_DTM_PAD_mean,inventory_PAD,max_return=2):
    fig = plt.figure(figure_number, facecolor='White',figsize=[8,9])

    # First up - old-growth forest
    Plot_name='Belian'
    ids = np.arange(heights.size)
    idmax = ids[MacArthurHorn_PAD_mean[Plot_name]>0][-1]+1
    y=heights[2:idmax]
    x=np.mean(inventory_PAD[Plot_name],axis=0)[2:idmax]
    x1=MacArthurHorn_PAD_mean[Plot_name][2:idmax]
    x2=radiative_PAD_mean[Plot_name][:-1,max_return-1][::-1][2:idmax]
    x3=radiative_DTM_PAD_mean[Plot_name][:-1,max_return-1][::-1][2:idmax]

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
    idmax = ids[MacArthurHorn_PAD_mean[Plot_name]>0][-1]+1
    y=heights[2:idmax]
    x=np.mean(inventory_PAD[Plot_name],axis=0)[2:idmax]
    x1=MacArthurHorn_PAD_mean[Plot_name][2:idmax]
    x2=radiative_PAD_mean[Plot_name][:-1,max_return-1][::-1][2:idmax]
    x3=radiative_DTM_PAD_mean[Plot_name][:-1,max_return-1][::-1][2:idmax]

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
    idmax = ids[MacArthurHorn_PAD_mean[Plot_name]>0][-1]+1
    y=heights[2:idmax]
    x=np.mean(inventory_PAD[Plot_name],axis=0)[2:idmax]
    x1=MacArthurHorn_PAD_mean[Plot_name][2:idmax]
    x2=radiative_PAD_mean[Plot_name][:-1,max_return-1][::-1][2:idmax]
    x3=radiative_DTM_PAD_mean[Plot_name][:-1,max_return-1][::-1][2:idmax]

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

"""
# Plot allometric relationships
# Plot three subplots defining the relationships used to construct the canopy profiles
"""
def plot_allometric_relationships(figure_name,figure_number,field_file,
                                    allometry_file,niter=10000):

    DBH_BAAD, H_BAAD, D_BAAD = field.load_BAAD_crown_allometry_data(allometry_file)
    #BAAD = field.load_BAAD_allometry_data(allometry_file,filter_dbh=True,filter_status='None')
    #H_BAAD = BAAD['ht']
    #D_BAAD = BAAD['cd']
    a, b, CF, r_sq, p, H_, PI_u, PI_l = field.log_log_linear_regression(H_BAAD,D_BAAD,conf=0.90)

    field_data = field.load_crown_survey_data(field_file)
    a_Ht, b_Ht, CF_Ht, r_sq_Ht, p_Ht, DBH_Ht, PI_u_Ht, PI_l_Ht = field.log_log_linear_regression(field_data['DBH_field'],field_data['Height_field'],conf=0.90)

    a_A, b_A, CF_A, r_sq_A, p_A, DBH_A, PI_u_A, PI_l_A = field.log_log_linear_regression(field_data['DBH_field'],field_data['CrownArea'],conf=0.90)

    PI_25_D,PI_75_D = field.calculate_powerlaw_prediction_interval_bootstrap_resampling_residuals(H_,H_BAAD,D_BAAD,conf=0.5,niter=niter)
    PI_05_D,PI_95_D = field.calculate_powerlaw_prediction_interval_bootstrap_resampling_residuals(H_,H_BAAD,D_BAAD,conf=0.9,niter=niter)

    mask = np.all((np.isfinite(field_data['DBH_field']),np.isfinite(field_data['CrownArea'])),axis=0)
    PI_25_A,PI_75_A = field.calculate_powerlaw_prediction_interval_bootstrap_resampling_residuals(DBH_A,field_data['DBH_field'][mask],field_data['CrownArea'][mask],conf=0.5,niter=niter)
    PI_05_A,PI_95_A = field.calculate_powerlaw_prediction_interval_bootstrap_resampling_residuals(DBH_A,field_data['DBH_field'][mask],field_data['CrownArea'][mask],conf=0.9,niter=niter)

    mask = np.all((np.isfinite(field_data['DBH_field']),np.isfinite(field_data['Height_field'])),axis=0)
    PI_25_H,PI_75_H = field.calculate_powerlaw_prediction_interval_bootstrap_resampling_residuals(DBH_Ht,field_data['DBH_field'][mask],field_data['Height_field'][mask],conf=0.5,niter=niter)
    PI_05_H,PI_95_H = field.calculate_powerlaw_prediction_interval_bootstrap_resampling_residuals(DBH_Ht,field_data['DBH_field'][mask],field_data['Height_field'][mask],conf=0.9,niter=niter)

    D_ = CF*a*H_**b
    CA_ = CF_A*a_A*DBH_A**b_A
    Ht_ = CF_Ht*a_Ht*DBH_Ht**b_Ht

    plt.figure(figure_number,facecolor='White',figsize=[9,3])
    ax1 = plt.subplot2grid((1,3),(0,0))
    ax1.set_ylabel('Crown Depth / m')
    ax1.set_xlabel('Height / m')
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

    ax1.set_xlim(xmin=0)
    ax2.set_xlim(xmin=0)
    ax3.set_xlim(xmin=0)
    ax1.set_ylim(ymin=0)
    ax2.set_ylim(ymin=0)
    ax3.set_ylim(ymin=0)
    plt.tight_layout()
    plt.savefig(figure_name)
    plt.show()
    return 0


"""
# Plot transmittance ratio - the ratio of observed vs expected number of subsequent
# LiDAR returns for a single pulse
"""
def plot_transmittance_ratio(figure_number,figure_name,pts):
    # A figure illustrating transmittance ratio between successive returns
    plt.figure(figure_number, facecolor='White',figsize=[3,3])
    ax = plt.subplot2grid((1,1),(0,0))
    ax.set_xlabel('return number, k')
    ax.set_ylabel(r'$\gamma$')
    for i in range(0,4):
        if i==0:
            ax.plot(1,1,'o',color='#1A2BCE')
        else:
            N_veg = float(np.all((pts[:,3]==i,pts[:,4]==1),axis=0).sum())
            N_i = float((pts[:,3]==i+1).sum())
            ax.plot(i+1,N_i/N_veg,'o',color='#1A2BCE')

    ax.set_ylim(0,1.1)
    ax.set_xlim(0.2,4.8)
    plt.xticks(np.arange(4)+1)
    plt.tight_layout()
    plt.savefig(figure_name)
    plt.show()
    return 0

"""
# Plot comparison of LiDAR derived profiles: Mac-Horn, Detto, Detto-corrected
"""
def plot_LiDAR_profiles_comparison(figure_name,figure_number,heights,heights_rad,
                        lidar_profiles,MacArthurHorn_PAD,MacArthurHorn_PAD_mean,
                        radiative_PAD,radiative_PAD_mean,
                        radiative_DTM_PAD,radiative_DTM_PAD_mean):
    max_return=3
    n_subplots = 25
    fig_plots = [b'Belian',b'E',b'B South']

    plt.figure(figure_number, facecolor='White',figsize=[8,8])

    # Belian
    ax1a = plt.subplot2grid((3,4),(0,0))
    ax1a.annotate('a -  MLA01\nOld growth', xy=(0.05,0.95), xycoords='axes fraction',
                backgroundcolor='none',horizontalalignment='left', verticalalignment='top',
                fontsize=10)
    ax1a.annotate('Return profile', xy=(0.95,0.05), xycoords='axes fraction',
                backgroundcolor='none',horizontalalignment='right',
                verticalalignment='bottom', fontsize=9)
    ax2a = plt.subplot2grid((3,4),(0,1),sharey=ax1a)
    ax2a.annotate('MacArthur-Horn', xy=(0.95,0.95), xycoords='axes fraction',
                backgroundcolor='none',horizontalalignment='right',
                verticalalignment='top', fontsize=9)
    ax3a = plt.subplot2grid((3,4),(0,2),sharey=ax1a,sharex=ax2a)
    ax3a.annotate('multi. return \nrad. trans.\npre-correction', xy=(0.95,0.95),
                xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right',
                verticalalignment='top', fontsize=9)
    ax4a = plt.subplot2grid((3,4),(0,3),sharey=ax1a,sharex=ax2a)
    ax4a.annotate('multi. return \nrad. trans.\npost-correction', xy=(0.95,0.95),
                xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right',
                verticalalignment='top', fontsize=9)

    # E
    ax1b = plt.subplot2grid((3,4),(1,0),sharey=ax1a,sharex=ax1a)
    ax1b.annotate('b - SAF05\nModerately logged', xy=(0.05,0.95), xycoords='axes fraction',
                backgroundcolor='none',horizontalalignment='left', verticalalignment='top',
                fontsize=10)
    ax2b = plt.subplot2grid((3,4),(1,1),sharey=ax1a,sharex=ax2a)
    ax3b = plt.subplot2grid((3,4),(1,2),sharey=ax1a,sharex=ax3a)
    ax4b = plt.subplot2grid((3,4),(1,3),sharey=ax1a,sharex=ax3a)

    # B South
    ax1c = plt.subplot2grid((3,4),(2,0),sharey=ax1a,sharex=ax1a)
    ax1c.annotate('c - SAF01\nHeavily logged', xy=(0.05,0.95), xycoords='axes fraction',
                backgroundcolor='none',horizontalalignment='left', verticalalignment='top',
                fontsize=10)
    ax2c = plt.subplot2grid((3,4),(2,1),sharey=ax1a,sharex=ax2a)
    ax3c = plt.subplot2grid((3,4),(2,2),sharey=ax1a,sharex=ax3a)
    ax4c = plt.subplot2grid((3,4),(2,3),sharey=ax1a,sharex=ax3a)

    ax1a.set_ylabel('Height / m',fontsize=axis_size)
    ax1b.set_ylabel('Height / m',fontsize=axis_size)
    ax1c.set_ylabel('Height / m',fontsize=axis_size)
    ax1c.set_xlabel('Number of returns\n(x1000)',fontsize=axis_size,horizontalalignment='center')
    ax2c.set_xlabel('PAD\n(m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size,horizontalalignment='center')
    ax3c.set_xlabel('PAD\n(m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size,horizontalalignment='center')
    ax4c.set_xlabel('PAD\n(m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size,horizontalalignment='center')
    #---------------------------------------------------------

    axes1 = [ax1a,  ax1b, ax1c]
    axes2 = [ax2a,  ax2b, ax2c]
    axes3 =  [ax3a,  ax3b, ax3c]
    axes4 = [ax4a,  ax4b, ax4c]

    yticklabels=[]
    xticklabels=[]

    for pp in range(0,3):
        Plot_name = fig_plots[pp]

        # plot lidar profile
        return_dist     = np.sum(lidar_profiles[Plot_name],axis=0)
        for k in range(0,max_return):
            axes1[pp].plot(return_dist[:,k]/1000.,np.max(heights_rad)-heights_rad,
            '-',c=colour[k], linewidth=1)

        # plot macarthur horn profile
        for i in range(0,n_subplots):
            axes2[pp].fill_betweenx(heights[2:],0,MacArthurHorn_PAD[Plot_name][i,2:],
                                    color=colour[0],alpha=0.02)
        axes2[pp].plot(MacArthurHorn_PAD_mean[Plot_name][2:],heights[2:],'-',c=colour[0],linewidth=2)

        # plot original Detto radiative transfer profile
        for i in range(0,n_subplots):
            axes3[pp].fill_betweenx(heights_rad[3:],0,radiative_PAD[Plot_name][i,:-3,-1][::-1],
                                    color=colour[1],alpha=0.02)
        axes3[pp].plot(radiative_PAD_mean[Plot_name][:-3,1][::-1],heights_rad[3:],
                                    '-',c=colour[1],linewidth=2)

        # plot corrective radiative transfer profile
        for i in range(0,n_subplots):
            axes4[pp].fill_betweenx(heights_rad[3:],0,radiative_DTM_PAD[Plot_name][i,:-3,-1][::-1],
                                    color=colour[1],alpha=0.02)
        axes4[pp].plot(radiative_DTM_PAD_mean[Plot_name][:-3,1][::-1],heights_rad[3:],
                                    '-',c=colour[1],linewidth=2)


        yticklabels.append(axes2[pp].get_yticklabels())
        yticklabels.append(axes3[pp].get_yticklabels())
        yticklabels.append(axes4[pp].get_yticklabels())

        if pp < 2:
            xticklabels.append(axes1[pp].get_xticklabels())
            xticklabels.append(axes2[pp].get_xticklabels())
            xticklabels.append(axes3[pp].get_xticklabels())
            xticklabels.append(axes4[pp].get_xticklabels())

    ax1a.set_ylim(0,80)
    ax1a.set_xlim(0,29)
    ax2a.set_xlim(xmin=0,xmax=0.7)

    ax2c.locator_params(axis='x',nbins=5)
    ax3c.locator_params(axis='x',nbins=5)
    ax4c.locator_params(axis='x',nbins=5)

    plt.setp(yticklabels,visible=False)
    plt.setp(xticklabels,visible=False)
    plt.subplots_adjust(hspace=0.2, wspace = 0.1)

    plt.tight_layout()
    plt.savefig(figure_name)
    plt.show()

    return 0

"""
# Plot comparison of LiDAR derived profiles using stochastic radiative transfer
# assuming different leaf angle distributions
"""
def plot_leaf_angle_distribution_profile_comparison(figure_name,figure_number,heights,
                        spherical_PAD_mean,erectophile_PAD_mean,planophile_PAD_mean):
    max_return=2
    fig_plots = [b'Belian',b'E',b'B South']

    plt.figure(figure_number, facecolor='White',figsize=[7,4])

    # Belian
    ax1 = plt.subplot2grid((1,3),(0,0))
    ax1.annotate('a - MLA01\nOld growth', xy=(0.05,0.95), xycoords='axes fraction',
                backgroundcolor='none',horizontalalignment='left', verticalalignment='top',
                fontsize=10)
    ax2 = plt.subplot2grid((1,3),(0,1),sharey=ax1,sharex=ax1)
    ax2.annotate('b - SAF05\nmoderately logged', xy=(0.05,0.95), xycoords='axes fraction',
                backgroundcolor='none',horizontalalignment='left',
                verticalalignment='top', fontsize=9)
    ax3 = plt.subplot2grid((1,3),(0,2),sharey=ax1,sharex=ax1)
    ax3.annotate('c - SAF01\nheavily logged', xy=(0.05,0.95), xycoords='axes fraction',
                backgroundcolor='none',horizontalalignment='left',
                verticalalignment='top', fontsize=9)

    ax1.set_ylabel('Height / m',fontsize=axis_size)
    ax1.set_xlabel('PAD\n(m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size,horizontalalignment='center')
    ax2.set_xlabel('PAD\n(m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size,horizontalalignment='center')
    ax3.set_xlabel('PAD\n(m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size,horizontalalignment='center')
    #---------------------------------------------------------

    axes = [ax1,  ax2, ax3]
    yticklabels=[];yticklabels.append(ax2.get_yticklabels()); yticklabels.append(ax3.get_yticklabels())

    for pp in range(0,3):
        Plot_name = fig_plots[pp]
        if pp<2:
            axes[pp].plot(spherical_PAD_mean[Plot_name][:-2,max_return-1][::-1],heights[2:],'-',c=colour[0],linewidth=1)
            axes[pp].plot(planophile_PAD_mean[Plot_name][:-2,max_return-1][::-1],heights[2:],'-',c=colour[1],linewidth=1)
            axes[pp].plot(erectophile_PAD_mean[Plot_name][:-2,max_return-1][::-1],heights[2:],'-',c=colour[2],linewidth=1)

        else:
            axes[pp].plot(spherical_PAD_mean[Plot_name][:-2,max_return-1][::-1],heights[2:],'-',c=colour[0],linewidth=1,label='spherical')
            axes[pp].plot(planophile_PAD_mean[Plot_name][:-2,max_return-1][::-1],heights[2:],'-',c=colour[1],linewidth=1,label='planophile')
            axes[pp].plot(erectophile_PAD_mean[Plot_name][:-2,max_return-1][::-1],heights[2:],'-',c=colour[2],linewidth=1,label='erectophile')

    ax1.set_ylim(0,80)
    ax2.set_xlim(xmin=0,xmax=1)

    ax2.locator_params(axis='x',nbins=5)
    ax3.locator_params(axis='x',nbins=5)
    ax1.locator_params(axis='x',nbins=5)

    ax3.legend(loc = 'center right')

    plt.setp(yticklabels,visible=False)
    plt.subplots_adjust(hspace=0.2, wspace = 0.1)

    plt.tight_layout()
    plt.savefig(figure_name)
    plt.show()

    return 0



"""
# Plot example crown model
"""
def plot_canopy_model(figure_number,figure_name,Plot_name,field_data,angle,
                    a_ht, b_ht, CF_ht, a_A, b_A, CF_A, a_D, b_D, CF_D,
                    clip=True,buff=10):
    # set up crown model etc
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

    Ht,Area,Depth = field.calculate_crown_dimensions(field_data['DBH_field'][mask],
                    field_data['Height_field'][mask],field_data['CrownArea'][mask],
                    a_ht, b_ht, CF_ht, a_A, b_A, CF_A, a_D, b_D, CF_D)

    Rmax = np.sqrt(Area/np.pi)

    canopy = field.generate_3D_ellipsoid_canopy(x,y,z,x0,y0,Ht,Depth,Rmax)

    from scipy.ndimage import rotate as rot
    canopy_ = rot(canopy,angle)
    if clip:
        row_test = np.sum(np.sum(canopy_,axis=2),axis=1)
        col_test = np.sum(np.sum(canopy_,axis=2),axis=0)
        print(col_test.shape,row_test.shape)
        row_min = np.argwhere(row_test>0.001)[0][0]
        row_max = np.argwhere(row_test>0.001)[-1][0]
        col_min = np.argwhere(col_test>0.001)[0][0]
        col_max = np.argwhere(col_test>0.001)[-1][0]
        print(row_min,row_max,col_min,col_max)
        canopy_ = canopy_[row_min:row_max,col_min:col_max]

    dim_y = canopy_.shape[0]
    dim_z = canopy_.shape[2]

    # Now plot up the crown model
    fig, (ax1, ax2) = plt.subplots(2,1, gridspec_kw = {'height_ratios':[dim_y,dim_z]},figsize=(8,8),num=figure_number)

    ax1.annotate('a', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',
                horizontalalignment='left', verticalalignment='top', fontsize=10,color='white')
    ax2.annotate('b', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',
                horizontalalignment='left', verticalalignment='top',fontsize=10,color='white')
    ax1.set_ylabel('Horizontal distance / m',fontsize=axis_size)
    ax1.set_xlabel('Horizontal distance / m',fontsize=axis_size)
    ax2.set_ylabel('Height / m',fontsize=axis_size)
    ax2.set_xlabel('Horizontal distance / m',fontsize=axis_size)

    im1 = ax1.imshow(np.sum(canopy_,axis=2), cmap = 'viridis',vmin=0,origin = 'lower')
    im2 = ax2.imshow(np.sum(canopy_,axis=0).transpose()/100., cmap = 'viridis',vmin=0,origin = 'lower')

    # colorbar bits and bobs
    divider = make_axes_locatable(ax1)
    cax1 = divider.new_horizontal(size="5%", pad=0.8, axes_class=plt.Axes)
    fig.add_axes(cax1)
    cbar1=plt.colorbar(im1, cax=cax1)
    cbar1.ax.set_ylabel('Crown volume / m$^3$m$^{-2}$',fontsize = 10)
    cbar1.solids.set_edgecolor("face")

    divider = make_axes_locatable(ax2)
    cax2 = divider.new_horizontal(size="5%", pad=0.8, axes_class=plt.Axes)
    fig.add_axes(cax2)
    cbar2=plt.colorbar(im2, cax=cax2)
    cbar2.ax.set_ylabel('Crown volume density / m$^3$m$^{-2}$m$^{-1}$',fontsize = 10)
    cbar2.solids.set_edgecolor("face")
    plt.savefig(figure_name)
    plt.show()
    return 0

"""
Plot distributions of PAD for different canopy subdivisions
- understory: <15 m
- lower canopy: 15-30 m
- mid canopy: 30-45 m
- upper canopy: 45 m
"""
def plot_PAD_distributions_for_canopy_subdivisions(figure_name,figure_number,
                            heights,heights_rad,MacArthurHorn_PAD,
                            radiative_PAD,radiative_DTM_PAD):

    max_k = 2
    # Create pandas dataframes for upper canopy, mid canopy lower canopy and
    # understory PAD
    uc = np.zeros(8*25*2)
    lc = np.zeros(8*25*2)
    mc = np.zeros(8*25*2)
    und = np.zeros(8*25*2)
    plot = []
    method = np.zeros(8*25*2)
    plot_name_1 = [b'Belian',b'Seraya',b'DC1',b'DC2',b'LF',b'E',b'B North',b'B South']
    plot_name = ['MLA01','MLA02','DAN04','DAN05','SAF04','SAF03','SAF02','SAF01']
    plot_colour = ['#46E900','#46E900','#46E900','#46E900','#1A2BCE','#1A2BCE','#E0007F','#E0007F']
    cc=0
    for pp in range(0,8):
        for ss in range(0,25):
            for ii in range(0,2):
                if ii==0:
                    PAD = MacArthurHorn_PAD[plot_name_1[pp]][ss]
                    uc_mask = heights>=45
                    mc_mask = np.all((heights>=30,heights<45),axis=0)
                    lc_mask = np.all((heights>=15,heights<30),axis=0)
                    und_mask = heights<15

                    uc[cc]=np.nansum(PAD[uc_mask])
                    mc[cc]=np.nansum(PAD[mc_mask])
                    lc[cc]=np.nansum(PAD[lc_mask])
                    und[cc]=np.nansum(PAD[und_mask])
                    """
                    uc[cc]=np.nanmean(PAD[uc_mask])
                    mc[cc]=np.nanmean(PAD[mc_mask])
                    lc[cc]=np.nanmean(PAD[lc_mask])
                    und[cc]=np.nanmean(PAD[und_mask])
                    """

                else:
                    PAD = radiative_DTM_PAD[plot_name_1[pp]][ss,:,max_k-1]
                    uc_mask = heights_rad>=45
                    mc_mask = np.all((heights_rad>=30,heights_rad<45),axis=0)
                    lc_mask = np.all((heights_rad>=15,heights_rad<30),axis=0)
                    und_mask = heights_rad<15

                    uc[cc]=np.nansum(PAD[uc_mask[::-1]])
                    mc[cc]=np.nansum(PAD[mc_mask[::-1]])
                    lc[cc]=np.nansum(PAD[lc_mask[::-1]])
                    und[cc]=np.nansum(PAD[und_mask[::-1]])

                plot.append(plot_name[pp])
                method[cc]=ii
                cc+=1

    df = pd.DataFrame({'plot' : plot,'upper' : uc,'mid' : mc,
                            'lower':lc,'understory':und})
    df_mh = df[method==0]
    df_rad = df[method==1]

    fig = plt.figure(figure_number, facecolor='White',figsize=(7,8))

    axa = plt.subplot2grid((4,2),(0,0))
    axa.set_title('MacArthur-Horn', fontsize=10)
    axa.annotate('a - upper canopy (>45 m)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    sns.violinplot(x='plot',y='upper',data=df_mh,inner="box",linewidth=0.5,
                    scale='width',hue="plot",palette = plot_colour,dodge=False,saturation=0.7,cut=0)
    axa.set_xlabel('')
    axa.set_ylabel('')

    axb = plt.subplot2grid((4,2),(0,1),sharex=axa,sharey=axa)
    axb.set_title('multi return\nrad. trans.', fontsize=10)
    axb.annotate('b', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    sns.violinplot(x='plot',y='upper',data=df_rad,inner="box",linewidth=0.5,
                    scale='width',hue="plot",palette = plot_colour,dodge=False,saturation=0.7,cut=0)
    axb.set_xlabel('')
    axb.set_ylabel('')

    axc = plt.subplot2grid((4,2),(1,0),sharex=axa,sharey=axa)
    axc.annotate('c - mid canopy (30-45 m)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    sns.violinplot(x='plot',y='mid',data=df_mh,inner="box",linewidth=0.5,
                    scale='width',hue="plot",palette = plot_colour,dodge=False,saturation=0.7,cut=0)
    axc.set_xlabel('')
    axc.set_ylabel('')

    axd = plt.subplot2grid((4,2),(1,1),sharex=axa,sharey=axa)
    axd.annotate('d', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    sns.violinplot(x='plot',y='mid',data=df_rad,inner="box",linewidth=0.5,
                    scale='width',hue="plot",palette = plot_colour,dodge=False,saturation=0.7,cut=0)
    axd.set_xlabel('')
    axd.set_ylabel('')

    axe = plt.subplot2grid((4,2),(2,0),sharex=axa,sharey=axa)
    axe.annotate('e - lower canopy (15-30 m)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    sns.violinplot(x='plot',y='lower',data=df_mh,inner="box",linewidth=0.5,
                    scale='width',hue="plot",palette = plot_colour,dodge=False,saturation=0.7,cut=0)
    axe.set_xlabel('')
    axe.set_ylabel('')

    axf = plt.subplot2grid((4,2),(2,1),sharex=axa,sharey=axa)
    axf.annotate('f', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    sns.violinplot(x='plot',y='lower',data=df_rad,inner="box",linewidth=0.5,
                    scale='width',hue="plot",palette = plot_colour,dodge=False,saturation=0.7,cut=0)
    axf.set_xlabel('')
    axf.set_ylabel('')

    axg = plt.subplot2grid((4,2),(3,0),sharex=axa,sharey=axa)
    axg.annotate('g - understory (<15 m)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    sns.violinplot(x='plot',y='understory',data=df_mh,inner="box",linewidth=0.5,
                    scale='width',hue="plot",palette = plot_colour,dodge=False,saturation=0.7,cut=0)
    axg.set_xlabel('')
    axg.set_ylabel('')

    axh = plt.subplot2grid((4,2),(3,1),sharex=axa,sharey=axa)
    axh.annotate('h', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    sns.violinplot(x='plot',y='understory',data=df_rad,inner="box",linewidth=0.5,
                    scale='width',hue="plot",palette = plot_colour,dodge=False,saturation=0.7,cut=0)
    axh.set_xlabel('')
    axh.set_ylabel('')

    axes =[axa,axb,axc,axd,axe,axf,axg,axh]
    yticklabels = axb.get_yticklabels() + axd.get_yticklabels() + axf.get_yticklabels() + axh.get_yticklabels()
    plt.setp(yticklabels,visible=False)

    for ii,ax in enumerate(axes):
        ax.yaxis.set_ticks_position('both')
        ax.legend_.remove()
        for item in ax.get_xticklabels():
            item.set_rotation(90)
        if ii<6:
            plt.setp(ax.get_xticklabels(),visible=False)

    fig.text(0.04, 0.5, 'contribution to PAI m$^2$m$^{-2}$', va='center', rotation='vertical', fontsize=10)
    plt.subplots_adjust(wspace = 0.2,hspace=0.1)
    plt.savefig(figure_name)

    plt.show()
    return 0

"""
# Plot subcanopy environments
# Plot subplot level profiles showing the cumulative distribution of effective
# plant area as you move down through the canopy with depth. Depth is defined
# relative to the maximum canopy height in a given subplot. Subplot profiles are
# plotted across four subplots, corresponding to old growth - Maliau Basin, old
# growth - Danum VAlley, Moderately logged - SAFE and heavily logged - SAFE
# Optional arguments: method (0 = MacArthur-Horn; 1 = Radiative transfer)
"""
def plot_cumulative_PAD_vs_depth(figure_name,figure_number,MacArthurHorn_PAD,
                        radiative_DTM_PAD, method=0):
    max_return=2
    fig_plots = [[b'Belian',b'Seraya'],[b'DC1',b'DC2'],[b'LF',b'E'],[b'B North',b'B South']]

    plt.figure(figure_number, facecolor='White',figsize=[6,8])

    # Belian
    ax1 = plt.subplot2grid((2,2),(0,0))
    ax1.annotate('a - Maliau Basi; Old growth', xy=(0.95,0.95), xycoords='axes fraction',
                backgroundcolor='none',horizontalalignment='right', verticalalignment='top',
                fontsize=10)
    ax2 = plt.subplot2grid((2,2),(0,1),sharey=ax1,sharex=ax1)
    ax2.annotate('b - Danum Valley; Old growth', xy=(0.95,0.95), xycoords='axes fraction',
                backgroundcolor='none',horizontalalignment='right',
                verticalalignment='top', fontsize=9)
    ax3 = plt.subplot2grid((2,2),(1,0),sharey=ax1,sharex=ax1)
    ax3.annotate('c - SAFE; moderately logged', xy=(0.95,0.95), xycoords='axes fraction',
                backgroundcolor='none',horizontalalignment='right',
                verticalalignment='top', fontsize=9)
    ax4 = plt.subplot2grid((2,2),(1,1),sharey=ax1,sharex=ax1)
    ax4.annotate('d - SAFE; heavily logged', xy=(0.95,0.95), xycoords='axes fraction',
                backgroundcolor='none',horizontalalignment='right',
                verticalalignment='top', fontsize=9)

    ax1.set_ylabel('Depth in canopy / m',fontsize=axis_size)
    ax3.set_ylabel('Depth in canopy / m',fontsize=axis_size)
    ax3.set_xlabel('cumulative PAD\n(m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size,horizontalalignment='center')
    ax4.set_xlabel('cumulative PAD\n(m$^2$m$^{-2}$m$^{-1}$)',fontsize=axis_size,horizontalalignment='center')
    #---------------------------------------------------------

    axes = [ax1,  ax2, ax3, ax4]
    yticklabels=[];yticklabels.append(ax2.get_yticklabels()); yticklabels.append(ax4.get_yticklabels())
    xticklabels=[];xticklabels.append(ax1.get_xticklabels()); xticklabels.append(ax2.get_xticklabels())

    for pp in range(0,4):
        Plot_name = fig_plots[pp]
        depth = np.arange(0,80)
        if pp<2:
            col = colour[0]
        elif pp<3:
            col = colour[1]
        else:
            col=colour[2]
        for plot in fig_plots[pp]:
            if method == 0:
                PAD = MacArthurHorn_PAD[plot][:,2:][:,::-1]
                for ss in range(0,25):
                    cumulative_PAD = np.cumsum(PAD[ss])
                    mask = cumulative_PAD>0
                    axes[pp].plot(cumulative_PAD[mask],depth[:mask.sum()],'-',color=col,
                                    alpha=0.25,linewidth=1)
            elif method ==1:
                PAD = radiative_DTM_PAD[plot][:,:-2,max_return-1]
                for ss in range(0,25):
                    cumulative_PAD = np.cumsum(PAD[ss])
                    mask = cumulative_PAD>0
                    axes[pp].plot(cumulative_PAD[mask],depth[:mask.sum()],'-',color=col,
                                    alpha=0.25,linewidth=1)

    ax1.set_ylim(0,80)
    ax1.set_xlim(0,13)
    ax1.invert_yaxis()

    for ax in axes:
        ax.locator_params(axis='x',nbins=5)

    plt.setp(yticklabels,visible=False)
    plt.setp(xticklabels,visible=False)
    plt.subplots_adjust(hspace=0.2, wspace = 0.1)

    plt.tight_layout()
    plt.savefig(figure_name)
    plt.show()


"""
Plot canopy Shannon Index and PAI
Violin plot for canopy shannon index and PAI for different methods
"""
def plot_PAI_Shannon_Index_distributions(figure_name,figure_number,PAD):
    n_plots = len(PAD.keys())
    Shannon = np.zeros(n_subplots*n_plots)
    PAI = np.zeros(n_subplots*n_plots)
    subplots = np.zeros(n_subplots*n_plots)
    plots = []
    region = []
    forest_type = []
    for pp,plot in enumerate(PAD.keys()):
        subplots[pp*n_subplots:(pp+1)*n_subplots]= np.arange(n_subplots)+1
        plots+=[plot]*n_subplots
        if plot in [b'Belian',b'Seraya']:
            region += ['Maliau\nBasin']*n_subplots
            forest_type += ['OG']*n_subplots
        elif plot in [b'DC1',b'DC2']:
            region += ['Danum\nValley']*n_subplots
            forest_type += ['OG']*n_subplots
        elif plot in [b'LF',b'E']:
            region += ['SAFE\nmoderately\nlogged']*n_subplots
            forest_type += ['ML']*n_subplots
        else:
            region += ['SAFE\nheavily\nlogged']*n_subplots
            forest_type += ['HL']*n_subplots
        P = PAD[plot][:,2:]
        for ii in range(P.shape[1]):
            P[np.isnan(P[:,ii])]=np.mean(P[np.isfinite(P[:,ii])])
        for ii in range(n_subplots):
            Shannon[pp*n_subplots+ii]=metrics.calculate_Shannon_index(P[ii])
            PAI[pp*n_subplots+ii]=np.nansum(P[ii])
    df = pd.DataFrame({'Plot':plots,'Subplot':subplots,'ShannonIndex':Shannon,
                        'PAI':PAI,'Region':region,'forest_type':forest_type})
    palette_colour = ['#46E900','#1A2BCE','#E0007F']

    # Now plot the distributions for the two approaches
    # distributions = plot
    # hue = region
    # facet = method
    fig = plt.figure(figure_number, facecolor='White',figsize=(7,3))

    axa = plt.subplot2grid((1,2),(0,0))
    axa.annotate('a - PAI', xy=(0.95,0.95), xycoords='axes fraction',
                    backgroundcolor='none',horizontalalignment='right',
                    verticalalignment='top', fontsize=10)
    sns.violinplot(x='Region',y='PAI',data=df,inner="box",linewidth=0.5,
                    scale='width',hue="forest_type",hue_order=['OG','ML','HL'],
                    order=['Maliau\nBasin', 'Danum\nValley',
                    'SAFE\nmoderately\nlogged','SAFE\nheavily\nlogged'],
                    palette = palette_colour,dodge=False,saturation=0.7,cut=0)
    axa.set_xlabel('')
    axa.set_ylabel('PAI / m$^2$m$^{-2}$')

    axb = plt.subplot2grid((1,2),(0,1))
    axb.annotate('b - Shannon Index', xy=(0.95,0.95), xycoords='axes fraction',
                    backgroundcolor='none',horizontalalignment='right',
                    verticalalignment='top', fontsize=10)
    sns.violinplot(x='Region',y='ShannonIndex',data=df,inner="box",linewidth=0.5,
                    scale='width',hue="forest_type",hue_order=['OG','ML','HL'],
                    order=['Maliau\nBasin', 'Danum\nValley',
                    'SAFE\nmoderately\nlogged','SAFE\nheavily\nlogged'],
                    palette = palette_colour,dodge=False,saturation=0.7,cut=0)
    axb.set_xlabel('')
    axb.set_ylabel('Shannon Index')

    axes =[axa,axb]
    for ii,ax in enumerate(axes):
        ax.yaxis.set_ticks_position('both')
        for item in ax.get_xticklabels():
            item.set_rotation(90)


    axa.legend_.remove()

    plt.subplots_adjust(wspace = 0.4,hspace=0.1,bottom=0.3)
    plt.savefig(figure_name)
    plt.show()
    return 0


"""
Plot histograms of cumulative PAD
Summarise subcanopy environments by calculating the volume occupied by
sub-canopy environments with different levels of overstory vegetation,
characterised by differing levels of PAD. This plot should highlight which
environmental niches are particularly sensitive to disturbance across the
disturbance gradient
"""
def plot_cumulative_PAD_histograms(figure_name,figure_number,PAD,heights,
                        horizontal_resolution = 20):

    vertical_resolution = 1
    fig_plots = [[b'Belian',b'Seraya'],[b'DC1',b'DC2'],[b'LF',b'E'],[b'B North',b'B South']]

    n_plots = len(PAD.keys())
    n_layers = heights.size-2
    cumulative_PAD = np.zeros(n_subplots*n_plots*n_layers)
    plots = []
    region = []
    forest_type = []

    # Loop through the plots calculating cumulative PAD for each subplot
    for pp,plot in enumerate(PAD.keys()):

        plots+=[plot]*n_subplots*n_layers
        if plot in [b'Belian',b'Seraya']:
            region += ['Maliau\nBasin']*n_subplots*n_layers
            forest_type += ['OG']*n_subplots*n_layers
        elif plot in [b'DC1',b'DC2']:
            region += ['Danum\nValley']*n_subplots*n_layers
            forest_type += ['OG']*n_subplots*n_layers
        elif plot in [b'LF',b'E']:
            region += ['SAFE\nmoderately\nlogged']*n_subplots*n_layers
            forest_type += ['ML']*n_subplots*n_layers
        else:
            region += ['SAFE\nheavily\nlogged']*n_subplots*n_layers
            forest_type += ['HL']*n_subplots*n_layers

        P = PAD[plot][:,2:]
        for ii in range(P.shape[1]):
            P[np.isnan(P[:,ii])]=np.mean(P[np.isfinite(P[:,ii])])
        for ii in range(n_subplots):
            cumulative_PAD[pp*n_subplots*n_layers+ii*n_layers:pp*n_subplots*n_layers+(ii+1)*n_layers] = np.cumsum(P[ii,::-1])
    # mask out areas above canopy (cuulative PAD=0)
    mask = cumulative_PAD>0
    df = pd.DataFrame({'Plot':np.asarray(plots)[mask],'cumulative_PAD':cumulative_PAD[mask],
                        'Region':np.asarray(region)[mask],'forest_type':np.asarray(forest_type)[mask]})

    palette_colour = ['#46E900','#1A2BCE','#E0007F']

    # Now plot cumulative PAD distributions by region, coloured by forest type
    fig = plt.figure('cal/val random',figsize=(4,3))
    fig.clf()
    ax = fig.add_subplot(1,1,1)
    ax = sns.distplot(df['cumulative_PAD'][df['Region']=='Maliau\nBasin'],
                        color = '#30A000',label='Maliau Basin',
                        kde=False,norm_hist=False,bins=np.arange(0,14,0.5),
                        ax=ax,hist_kws={'alpha':1})
    ax = sns.distplot(df['cumulative_PAD'][df['Region']=='Danum\nValley'],
                        color = palette_colour[0],label='Danum Valley',
                        kde=False,norm_hist=False,bins=np.arange(0,14,0.5),
                        ax=ax,hist_kws={'alpha':1})
    ax = sns.distplot(df['cumulative_PAD'][df['Region']=='SAFE\nmoderately\nlogged'],
                        color = palette_colour[1],label='SAFE (moderately logged)',
                        kde=False,norm_hist=False,bins=np.arange(0,14,0.5),
                        ax=ax,hist_kws={'alpha':1})
    ax = sns.distplot(df['cumulative_PAD'][df['Region']=='SAFE\nheavily\nlogged'],
                        color = palette_colour[2],label='SAFE (heavily logged)',
                        kde=False,norm_hist=False,bins=np.arange(0,14,0.5),
                        ax=ax,hist_kws={'alpha':1})
    ax = sns.distplot(df['cumulative_PAD'][df['Region']=='SAFE\nheavily\nlogged'],
                        color = palette_colour[2],
                        hist_kws={"histtype": "step", "linewidth": 1.2,'alpha':1},
                        kde=False,norm_hist=False,bins=np.arange(0,14,0.5),ax=ax)
    ax = sns.distplot(df['cumulative_PAD'][df['Region']=='SAFE\nmoderately\nlogged'],
                        color = palette_colour[1],
                        hist_kws={"histtype": "step", "linewidth": 1.2,'alpha':1},
                        ax=ax,kde=False,norm_hist=False,bins=np.arange(0,14,0.5))
    ax = sns.distplot(df['cumulative_PAD'][df['Region']=='Danum\nValley'],
                        color = palette_colour[0],
                        hist_kws={"histtype": "step", "linewidth": 1.2,'alpha':1}, ax=ax,
                        kde=False,norm_hist=False,bins=np.arange(0,14,0.5),)
    ax = sns.distplot(df['cumulative_PAD'][df['Region']=='Maliau\nBasin'],
                        color = '#30A000',
                        hist_kws={"histtype": "step", "linewidth": 1.2,'alpha':1},
                        kde=False,norm_hist=False,bins=np.arange(0,14,0.5),ax=ax)
    ax.set_xlabel('Overlying plant area / m$^2$m$^{-2}$')
    ax.set_ylabel('Subcanopy volume / m$^3$ m$^{-2}$')
    y_ticks = ax.get_yticks()
    ax.set_yticklabels(['{:3.0f}'.format(i*vertical_resolution/(2.*n_subplots)) for i in y_ticks])
    ax.set_xlim((0,14))
    ax.legend()
    fig.tight_layout()
    fig.show()
