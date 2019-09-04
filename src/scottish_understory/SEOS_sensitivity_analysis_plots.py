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

"""
Plot sensitivity of 1ha profiles to spatial resolution
"""
def plot_profile_sensitivity_resolution(figure_name,heights,PAD_profiles):

    fig,axes = plt.subplots(1,5,sharey="all",sharex="all",figsize=[8,3])

    annotations = ['2 m','5 m','10 m','20 m','50 m']
    # loop through the subplots
    keys = ['2m','5m','10m','20m','50m']
    for ii,ax in enumerate(axes):
        ax.set_title(annotations[ii])

        profile_sp = np.nanmean(PAD_profiles[keys[ii]],axis=1)
        profile = np.nanmedian(profile_sp,axis=0)
        perc = np.nanpercentile(profile_sp,[2.5,97.5],axis=0)
        ax.fill_betweenx(heights[2:],perc[0][2:],perc[1][2:],color=colour[0],alpha=0.3)
        ax.plot(profile[2:],heights[2:],'-',c=colour[0],linewidth=1.5)
        print('%.3f' % np.sum(profile[2:]))

    axes[0].set_ylim(0,40)
    axes[0].set_xlim(0,0.35)
    for ax in axes:
        ax.locator_params(axis='x',nbins=4)

    plt.subplots_adjust(hspace=0.4, wspace = 0.4, bottom = 0.2)

    axes[0].set_ylabel('height / m',fontsize=axis_size)
    fig.text(0.5, 0.05, 'PAD / m$^2$m$^{-3}$',fontsize=axis_size, ha='center')

    fig.savefig(figure_name)
    fig.show()
    return 0


def plot_profile_sensitivity_density(figure_name,heights,PAD_profiles):

    fig,axes = plt.subplots(1,5,sharey="all",sharex="all",figsize=[8,3])
    annotations = ['20 shots m$^{-2}$','50 shots m$^{-2}$',
                    '100 shots m$^{-2}$','500 shots m$^{-2}$',
                    '1000 shots m$^{-2}$']
    # loop through the subplots
    keys = ['20','50','100','500','1000']
    for ii,ax in enumerate(axes):
        ax.set_title(annotations[ii])
        profile_sp = np.nanmean(PAD_profiles[keys[ii]],axis=1)
        profile = np.nanmedian(profile_sp,axis=0)
        perc = np.nanpercentile(profile_sp,[2.5,97.5],axis=0)
        ax.fill_betweenx(heights[2:],perc[0][2:],perc[1][2:],color=colour[0],alpha=0.3)
        ax.plot(profile[2:],heights[2:],'-',c=colour[0],linewidth=1.5)
        print('%.3f' % np.sum(profile[2:]))

    axes[0].set_ylim(0,40)
    axes[0].set_xlim(0,0.35)
    for ax in axes:
        ax.locator_params(axis='x',nbins=5)

    plt.subplots_adjust(hspace=0.4, wspace = 0.4, bottom = 0.2)

    axes[0].set_ylabel('height / m')
    fig.text(0.5, 0.05, 'PAD / m$^2$m$^{-3}$', ha='center')
    fig.savefig(figure_name)
    fig.show()
    return 0


def plot_profile_sensitivity_density_individual_profile(figure_name,heights,PAD_profiles, profile_idx=np.nan):
    if np.isnan(profile_idx):
        profile_idx = np.random.randint(PAD_profiles['20'].shape[1])


    fig,axes = plt.subplots(1,5,sharey="all",sharex="all",figsize=[8,3])
    annotations = ['20 shots m$^{-2}$','50 shots m$^{-2}$',
                    '100 shots m$^{-2}$','500 shots m$^{-2}$',
                    '1000 shots m$^{-2}$']
    # loop through the subplots
    keys = ['20','50','100','500','1000']
    for ii,ax in enumerate(axes):
        ax.set_title(annotations[ii])
        profile_sp = PAD_profiles[keys[ii]][:,profile_idx,:]
        profile = np.nanmedian(profile_sp,axis=0)
        perc = np.nanpercentile(profile_sp,[2.5,97.5],axis=0)
        ax.fill_betweenx(heights[2:],perc[0][2:],perc[1][2:],color=colour[0],alpha=0.3)
        ax.plot(profile[2:],heights[2:],'-',c=colour[0],linewidth=1.5)
        print('%.3f' % np.sum(profile[2:]))

    axes[0].set_ylim(0,40)
    axes[0].set_xlim(0,1.05)
    for ax in axes:
        ax.locator_params(axis='x',nbins=3)

    plt.subplots_adjust(hspace=0.4, wspace = 0.4, bottom = 0.2)

    axes[0].set_ylabel('height / m')
    fig.text(0.5, 0.05, 'PAD / m$^2$m$^{-3}$', ha='center')
    fig.savefig(figure_name)
    fig.show()
    return 0
