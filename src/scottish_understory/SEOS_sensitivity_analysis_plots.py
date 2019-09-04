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
def plot_profile_sensitivity_resolution(figure_number,figure_name,heights,PAD_profiles_profile):

    ## Belian
    fig = plt.figure(figure_number, facecolor='White',figsize=[8,4])
    axa = plt.subplot2grid((1,5),(0,0))
    axb = plt.subplot2grid((1,5),(0,1),sharex=axa,sharey=axa)
    axc = plt.subplot2grid((1,5),(0,2),sharex=axa,sharey=axa)
    axd = plt.subplot2grid((1,5),(0,3),sharex=axa,sharey=axa)
    axe = plt.subplot2grid((1,5),(0,4),sharex=axa,sharey=axa)

    annotations = ['a - 2 m','b - 5 m','c - 10 m','d - 20 m','e - 50 m']
    axes = [axa,axb,axc,axd,axe]
    # loop through the subplots
    keys = ['2m','5m','10m','20m','50m']
    for ii,ax in enumerate(axes):
        ax.annotate(annotations[ii], xy=(0.05,0.95), xycoords='axes fraction',
                    backgroundcolor='none',ha='left', va='top', fontsize=10)
        #rad_sp = np.nanmean(PAD_profiles_rad[keys[ii]],axis=1)
        profile_sp = np.nanmean(PAD_profiles_profile[keys[ii]],axis=1)
        profile = np.nanmedian(profile_sp,axis=0)
        perc = np.nanpercentile(profile_sp,[2.5,97.5],axis=0)
        ax.fill_betweenx(heights[2:],perc[0][2:],perc[1][2:],color=colour[0],alpha=0.3)

        if ii==4:
            ax.plot(profile[2:],heights[2:],'-',c=colour[0],linewidth=1.5, label = 'MacArthur-Horn')
            #ax.legend(loc=9, bbox_to_anchor=(-0.1, -0.15), ncol=3)
        else:
            ax.plot(profile[2:],heights[2:],'-',c=colour[0],linewidth=1.5)
        print('%.3f' % np.sum(profile[2:]))
    yticklabels = axb.get_yticklabels() + axc.get_yticklabels() + axd.get_yticklabels() + axe.get_yticklabels()
    plt.setp(yticklabels,visible=False)

    axa.set_ylim(0,40)
    axa.set_xlim(0,0.38)
    for ax in axes:
        ax.locator_params(axis='x',nbins=4)
        ax.yaxis.set_ticks_position('both')

    plt.subplots_adjust(hspace=0.4, wspace = 0.2, bottom = 0.2)

    axa.set_ylabel('height / m',fontsize=axis_size)
    fig.text(0.5, 0.05, 'PAD / m$^2$m$^{-3}$',fontsize=axis_size, ha='center')

    fig.savefig(figure_name)
    #fig.show()
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

        if ii==4:
            ax.plot(profile[2:],heights[2:],'-',c=colour[0],linewidth=1.5, label = 'MacArthur-Horn')
        else:
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

        if ii==4:
            ax.plot(profile[2:],heights[2:],'-',c=colour[0],linewidth=1.5, label = 'MacArthur-Horn')
            #ax.legend(loc=9, bbox_to_anchor=(-0.1, -0.15), ncol=3)
        else:
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
