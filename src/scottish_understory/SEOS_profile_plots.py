"""
SEOS_lidar_plots.py
--------------------------------------------------------------------------------
Plots for SEOS LiDAR surveys illustrating canopy profiles
"""
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rcParams
import matplotlib.ticker as mticker
import matplotlib.ticker as plticker
from matplotlib import colorbar
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import seaborn as sns
import pandas as pd
sns.set()

import sys
sys.path.append('../')
import LiDAR_MacHorn_LAD_profiles as LAD

colour = ['#46E900','#1A2BCE','#E0007F']
cmap = cm.get_cmap('viridis')
scale = np.arange(0.,4.)
scale /=3.3
colour_loop = cmap(scale)

"""
plot_point_cloud_and_profiles
Subplots showing:
- point cloud transect
- point cloud distribution with height
- estimated canopy profile, 1 ha average
- selection of individual profiles (5m x 5m for SEOS UAV surveys) with confidence
  intervals from sensitivity analysis
"""
def plot_point_cloud_and_profiles(figure_name, pts, PAD_profiles,heights,return_profile):

    # filter points to 5m transect * 100 m
    pts[:,0]=pts[:,0]-np.min(pts[:,0])
    pts[:,1]=pts[:,1]-np.min(pts[:,1])
    pts=pts[pts[:,1]<=50]
    pts=pts[pts[:,1]>=45]

    # plot things up
    # - point cloud
    fig=plt.figure(figsize=[8,8])
    axes=[]
    axes.append(fig.add_subplot(2,1,1))
    axes.append(fig.add_subplot(2,3,4))
    axes.append(fig.add_subplot(2,3,5,sharey=axes[0]))
    axes.append(fig.add_subplot(2,3,6,sharey=axes[0]))
    axes[0].scatter(pts[:,0],pts[:,2],marker='.',s=0.05,c=pts[:,2],cmap='viridis',alpha=0.1)
    # - point distribution
    axes[1].plot(return_profile/10.**4,heights,'-',color=colour[0])
    # - canopy profile (1 ha)
    profile_sp = np.nanmean(PAD_profiles,axis=1)
    profile = np.nanmedian(profile_sp,axis=0)
    perc = np.nanpercentile(profile_sp,[2.5,97.5],axis=0)
    axes[2].fill_betweenx(heights[2:],perc[0][2:],perc[1][2:],color=colour[0],
                            alpha=0.3)
    axes[2].plot(profile[2:],heights[2:],'-',c=colour[0],linewidth=1.5)
    # canopy profile (4 different individual profiles)
    for ii in range(0,4):
        profile_idx = np.random.randint(PAD_profiles.shape[1])
        profile_sp = PAD_profiles[:,profile_idx,:]
        profile = np.nanmedian(profile_sp,axis=0)
        perc = np.nanpercentile(profile_sp,[2.5,97.5],axis=0)
        axes[3].fill_betweenx(heights[2:],perc[0][2:],perc[1][2:],
                                color=colour_loop[ii],alpha=0.3)
        axes[3].plot(profile[2:],heights[2:],'-',c=colour_loop[ii],linewidth=1.5)

    # label up axes
    axes[0].set_aspect('equal')
    axes[0].set_ylabel('height / m')
    axes[1].set_ylabel('height / m')
    x_labels = ['distance / m','return density / m$^{-2}$', 'PAD / m$^2$m$^{-3}$', 'PAD / m$^2$m$^{-3}$']
    annotations = ['point cloud','return profile','1 ha average\ncanopy profile','individual 5m x 5m\ncanopy profile']

    for ii,ax in enumerate(axes):
        ax.annotate(annotations[ii],xy=(0.5,0.95), xycoords='axes fraction',
                    backgroundcolor='none',ha='center', va='top')
        ax.set_xlabel(x_labels[ii])
        ax.locator_params(axis='y',nbins=5)
        if ii > 0:
            ax.locator_params(axis='x',nbins=5)
        #if ii > 1:
            #ax.set_yticklabels(ax.get_yticklabels(),visible=False)
    axes[2].set_yticklabels(axes[2].get_yticklabels(),visible=False)
    axes[3].set_yticklabels(axes[3].get_yticklabels(),visible=False)
    fig.tight_layout()
    fig.show()
    fig.savefig(figure_name)
