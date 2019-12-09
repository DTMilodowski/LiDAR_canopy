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

rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = 8
rcParams['legend.numpoints'] = 1
axis_size = rcParams['font.size']+2
colour = ['#46E900','#1A2BCE','#E0007F']
n_subplots=25

import sys
sys.path.append('/home/dmilodow/DataStore_DTM/FOREST2020/EOdata/EO_data_processing/src/plot_EO_data/colormap/')
import colormaps as cmaps
plt.register_cmap(name='viridis', cmap=cmaps.viridis)
plt.register_cmap(name='plasma', cmap=cmaps.plasma)
plt.set_cmap(cmaps.viridis)

# code to get trio of nice colourblind friendly colours
cmap = cm.get_cmap('plasma')
scale = np.arange(0.,5.)
scale /=4.
cmap_colour = cmap(scale)

"""
Plot sensitivity of 1ha profiles to spatial resolution
"""
def plot_profile_sensitivity_resolution_full(figure_number,figure_name,heights,PAD_profiles_MH_Belian,
                            PAD_profiles_MH_BNorth,PAD_profiles_MH_E,PAD_profiles_rad_Belian,
                            PAD_profiles_rad_BNorth,PAD_profiles_rad_E,PAD_profiles_wt_Belian,
                            PAD_profiles_wt_BNorth,PAD_profiles_wt_E):
    ## Belian
    fig = plt.figure(figure_number, facecolor='White',figsize=[8,9])
    ax1a = plt.subplot2grid((3,5),(0,0))
    ax1a.set_title('MLA01 - old growth', fontsize=10)
    ax1a.set_ylabel('height / m',fontsize=axis_size)
    ax1a.annotate('a - 5 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1b = plt.subplot2grid((3,5),(0,1),sharex=ax1a,sharey=ax1a)
    ax1b.annotate('b - 10 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1c = plt.subplot2grid((3,5),(0,2),sharex=ax1a,sharey=ax1a)
    ax1c.annotate('c - 20 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1d = plt.subplot2grid((3,5),(0,3),sharex=ax1a,sharey=ax1a)
    ax1d.annotate('d - 50 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1e = plt.subplot2grid((3,5),(0,4),sharex=ax1a,sharey=ax1a)
    ax1e.annotate('e - 100 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ## E
    ax1f = plt.subplot2grid((3,5),(1,0),sharex=ax1a,sharey=ax1a)
    ax1f.set_title('SAF03 - moderately logged', fontsize=10)
    ax1f.set_ylabel('height / m',fontsize=axis_size)
    ax1f.annotate('f - 5 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1g = plt.subplot2grid((3,5),(1,1),sharex=ax1a,sharey=ax1a)
    ax1g.annotate('g - 10 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1h = plt.subplot2grid((3,5),(1,2),sharex=ax1a,sharey=ax1a)
    ax1h.annotate('h - 20 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1i = plt.subplot2grid((3,5),(1,3),sharex=ax1a,sharey=ax1a)
    ax1i.annotate('i - 50 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1j = plt.subplot2grid((3,5),(1,4),sharex=ax1a,sharey=ax1a)
    ax1j.annotate('j - 100 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ## B North
    ax1k = plt.subplot2grid((3,5),(2,0),sharex=ax1a,sharey=ax1a)
    ax1k.set_title('SAF02 - heavily logged', fontsize=10)
    ax1k.set_ylabel('height / m',fontsize=axis_size)
    ax1k.annotate('k - 5 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1l = plt.subplot2grid((3,5),(2,1),sharex=ax1a,sharey=ax1a)
    ax1l.annotate('l - 10 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1m = plt.subplot2grid((3,5),(2,2),sharex=ax1a,sharey=ax1a)
    ax1m.annotate('m - 20 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1n = plt.subplot2grid((3,5),(2,3),sharex=ax1a,sharey=ax1a)
    ax1n.annotate('n - 50 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1o = plt.subplot2grid((3,5),(2,4),sharex=ax1a,sharey=ax1a)
    ax1o.annotate('o - 100 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    # loop through the subplots
    sp = [ax1a,ax1b,ax1c,ax1d,ax1e,ax1f,ax1g,ax1h,ax1i,ax1j,ax1k,ax1l,ax1m,ax1n,ax1o]
    pkeys = ['5m','10m','20m','50m','100m']
    for pp in range(0,len(sp)):
        skip_rad = False
        if pp<5:
            if pkeys[pp]=='5m':
                skip_rad=True
            else:
                rad_sp = np.nanmean(PAD_profiles_rad_Belian[pkeys[pp]]['25'],axis=1)
            MH_sp = np.nanmean(PAD_profiles_MH_Belian[pkeys[pp]]['25'],axis=1)
            wt_sp = np.nanmean(PAD_profiles_wt_Belian[pkeys[pp]]['25'],axis=1)
        elif pp<10:
            if pkeys[pp-5]=='5m':
                skip_rad=True
            else:
                rad_sp = np.nanmean(PAD_profiles_rad_E[pkeys[pp-5]]['25'],axis=1)
            MH_sp = np.nanmean(PAD_profiles_MH_E[pkeys[pp-5]]['25'],axis=1)
            wt_sp = np.nanmean(PAD_profiles_wt_E[pkeys[pp-5]]['25'],axis=1)
        else:
            if pkeys[pp-10]=='5m':
                skip_rad=True
            else:
                rad_sp = np.nanmean(PAD_profiles_rad_BNorth[pkeys[pp-10]]['25'],axis=1)
            MH_sp = np.nanmean(PAD_profiles_MH_BNorth[pkeys[pp-10]]['25'],axis=1)
            wt_sp = np.nanmean(PAD_profiles_wt_BNorth[pkeys[pp-10]]['25'],axis=1)

        wt_sp[~np.isfinite(wt_sp)]=np.nan
        wt = np.nanmedian(wt_sp,axis=0)
        wtperc = np.nanpercentile(wt_sp,[2.5,97.5],axis=0)
        sp[pp].fill_betweenx(heights[2:],wtperc[0][2:],wtperc[1][2:],color=colour[2],alpha=0.3)
        if skip_rad == False:
            rad = np.nanmedian(rad_sp,axis=0)
            radperc = np.nanpercentile(rad_sp,[2.5,97.5],axis=0)
            sp[pp].fill_betweenx(heights[2:],radperc[0][2:],radperc[1][2:],color=colour[1],alpha=0.3)

        MH = np.nanmedian(MH_sp,axis=0)
        MHperc = np.nanpercentile(MH_sp,[2.5,97.5],axis=0)
        sp[pp].fill_betweenx(heights[2:],MHperc[0][2:],MHperc[1][2:],color=colour[0],alpha=0.3)

        if pp==14:
            sp[pp].plot(MH[2:],heights[2:],'-',c=colour[0],linewidth=1.5, label = 'M1a')
            sp[pp].plot(wt[2:],heights[2:],'-',c=colour[2],linewidth=1.5, label = 'M1b')
            sp[pp].plot(rad[2:],heights[2:],'-',c=colour[1],linewidth=1.5, label = 'M2')
            sp[pp].legend(loc=9, bbox_to_anchor=(-0.1, -0.15), ncol=3)
        else:
            sp[pp].plot(MH[2:],heights[2:],'-',c=colour[0],linewidth=1.5)
            sp[pp].plot(wt[2:],heights[2:],'-',c=colour[2],linewidth=1.5)
            if skip_rad == False:
                sp[pp].plot(rad[2:],heights[2:],'-',c=colour[1],linewidth=1.5)

    yticklabels = ax1b.get_yticklabels() + ax1c.get_yticklabels() + ax1d.get_yticklabels() + ax1e.get_yticklabels() + \
                  ax1g.get_yticklabels() + ax1h.get_yticklabels() + ax1i.get_yticklabels() + ax1j.get_yticklabels() + \
                  ax1l.get_yticklabels() + ax1m.get_yticklabels() + ax1n.get_yticklabels() + ax1o.get_yticklabels()
    xticklabels = ax1a.get_xticklabels() + ax1b.get_xticklabels() + ax1c.get_xticklabels() + ax1d.get_xticklabels() + ax1e.get_xticklabels() +  \
                  ax1f.get_xticklabels() + ax1g.get_xticklabels() + ax1h.get_xticklabels() + ax1i.get_xticklabels() + ax1j.get_xticklabels()

    plt.setp(yticklabels,visible=False)

    ax1a.set_ylim(0,89)
    ax1a.set_xlim(0,0.55)
    for ax in sp:
        ax.locator_params(axis='x',nbins=5)
        ax.yaxis.set_ticks_position('both')

    plt.subplots_adjust(hspace=0.4, wspace = 0.2, bottom = 0.2)

    fig.text(0.5, 0.15, 'PAD / m$^2$m$^{-3}$',fontsize=axis_size, ha='center')

    plt.savefig(figure_name)
    #plt.show()
    return 0


def plot_profile_sensitivity_resolution(figure_number,figure_name,heights,PAD_profiles_MH_Belian,
                            PAD_profiles_MH_BNorth,PAD_profiles_MH_E,PAD_profiles_rad_Belian,
                            PAD_profiles_rad_BNorth,PAD_profiles_rad_E):
    ## Belian
    fig = plt.figure(figure_number, facecolor='White',figsize=[8,9])
    ax1a = plt.subplot2grid((3,5),(0,0))
    ax1a.set_title('MLA01 - old growth', fontsize=10)
    ax1a.set_ylabel('height / m',fontsize=axis_size)
    ax1a.annotate('a - 5 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1b = plt.subplot2grid((3,5),(0,1),sharex=ax1a,sharey=ax1a)
    ax1b.annotate('b - 10 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1c = plt.subplot2grid((3,5),(0,2),sharex=ax1a,sharey=ax1a)
    ax1c.annotate('c - 20 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1d = plt.subplot2grid((3,5),(0,3),sharex=ax1a,sharey=ax1a)
    ax1d.annotate('d - 50 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1e = plt.subplot2grid((3,5),(0,4),sharex=ax1a,sharey=ax1a)
    ax1e.annotate('e - 100 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ## E
    ax1f = plt.subplot2grid((3,5),(1,0),sharex=ax1a,sharey=ax1a)
    ax1f.set_title('SAF03 - moderately logged', fontsize=10)
    ax1f.set_ylabel('height / m',fontsize=axis_size)
    ax1f.annotate('f - 5 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1g = plt.subplot2grid((3,5),(1,1),sharex=ax1a,sharey=ax1a)
    ax1g.annotate('g - 10 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1h = plt.subplot2grid((3,5),(1,2),sharex=ax1a,sharey=ax1a)
    ax1h.annotate('h - 20 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1i = plt.subplot2grid((3,5),(1,3),sharex=ax1a,sharey=ax1a)
    ax1i.annotate('i - 50 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1j = plt.subplot2grid((3,5),(1,4),sharex=ax1a,sharey=ax1a)
    ax1j.annotate('j - 100 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ## B North
    ax1k = plt.subplot2grid((3,5),(2,0),sharex=ax1a,sharey=ax1a)
    ax1k.set_title('SAF02 - heavily logged', fontsize=10)
    ax1k.set_ylabel('height / m',fontsize=axis_size)
    ax1k.annotate('k - 5 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1l = plt.subplot2grid((3,5),(2,1),sharex=ax1a,sharey=ax1a)
    ax1l.annotate('l - 10 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1m = plt.subplot2grid((3,5),(2,2),sharex=ax1a,sharey=ax1a)
    ax1m.annotate('m - 20 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1n = plt.subplot2grid((3,5),(2,3),sharex=ax1a,sharey=ax1a)
    ax1n.annotate('n - 50 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1o = plt.subplot2grid((3,5),(2,4),sharex=ax1a,sharey=ax1a)
    ax1o.annotate('o - 100 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    # loop through the subplots
    sp = [ax1a,ax1b,ax1c,ax1d,ax1e,ax1f,ax1g,ax1h,ax1i,ax1j,ax1k,ax1l,ax1m,ax1n,ax1o]
    pkeys = ['5m','10m','20m','50m','100m']
    for pp in range(0,len(sp)):

        if pp<5:
            rad_sp = np.nanmean(PAD_profiles_rad_Belian[pkeys[pp]]['25'],axis=1)
            MH_sp = np.nanmean(PAD_profiles_MH_Belian[pkeys[pp]]['25'],axis=1)
        elif pp<10:
            rad_sp = np.nanmean(PAD_profiles_rad_E[pkeys[pp-5]]['25'],axis=1)
            MH_sp = np.nanmean(PAD_profiles_MH_E[pkeys[pp-5]]['25'],axis=1)
        else:
            rad_sp = np.nanmean(PAD_profiles_rad_BNorth[pkeys[pp-10]]['25'],axis=1)
            MH_sp = np.nanmean(PAD_profiles_MH_BNorth[pkeys[pp-10]]['25'],axis=1)

        rad = np.nanmedian(rad_sp,axis=0)
        radperc = np.nanpercentile(rad_sp,[2.5,97.5],axis=0)
        sp[pp].fill_betweenx(heights[2:],radperc[0][2:],radperc[1][2:],color=colour[1],alpha=0.3)

        MH = np.nanmedian(MH_sp,axis=0)
        MHperc = np.nanpercentile(MH_sp,[2.5,97.5],axis=0)
        sp[pp].fill_betweenx(heights[2:],MHperc[0][2:],MHperc[1][2:],color=colour[0],alpha=0.3)

        if pp==17:
            sp[pp].plot(MH[2:],heights[2:],'-',c=colour[0],linewidth=1.5, label = 'MacArthur-Horn')
            sp[pp].plot(rad[2:],heights[2:],'-',c=colour[1],linewidth=1.5, label = 'multi. return\nrad trans')
            sp[pp].legend(loc=9, bbox_to_anchor=(-0.1, -0.15), ncol=3)
        else:
            sp[pp].plot(MH[2:],heights[2:],'-',c=colour[0],linewidth=1.5)
            sp[pp].plot(rad[2:],heights[2:],'-',c=colour[1],linewidth=1.5)

    yticklabels = ax1b.get_yticklabels() + ax1c.get_yticklabels() + ax1d.get_yticklabels() + ax1e.get_yticklabels() + \
                  ax1g.get_yticklabels() + ax1h.get_yticklabels() + ax1i.get_yticklabels() + ax1j.get_yticklabels() + \
                  ax1l.get_yticklabels() + ax1m.get_yticklabels() + ax1n.get_yticklabels() + ax1o.get_yticklabels()
    xticklabels = ax1a.get_xticklabels() + ax1b.get_xticklabels() + ax1c.get_xticklabels() + ax1d.get_xticklabels() + ax1e.get_xticklabels() +  \
                  ax1f.get_xticklabels() + ax1g.get_xticklabels() + ax1h.get_xticklabels() + ax1i.get_xticklabels() + ax1j.get_xticklabels()

    plt.setp(yticklabels,visible=False)

    ax1a.set_ylim(0,89)
    ax1a.set_xlim(0,0.55)
    for ax in sp:
        ax.locator_params(axis='x',nbins=5)
        ax.yaxis.set_ticks_position('both')

    plt.subplots_adjust(hspace=0.4, wspace = 0.2, bottom = 0.2)

    fig.text(0.5, 0.15, 'PAD / m$^2$m$^{-3}$',fontsize=axis_size, ha='center')

    plt.savefig(figure_name)
    #plt.show()
    return 0

def plot_profile_sensitivity_resolution_old(figure_number,figure_name,heights,PAD_profiles_MH_Belian,
                            PAD_profiles_MH_BNorth,PAD_profiles_MH_E,PAD_profiles_rad_Belian,
                            PAD_profiles_rad_BNorth,PAD_profiles_rad_E):
    ## Belian
    fig = plt.figure(figure_number, facecolor='White',figsize=[8,9])
    ax1a = plt.subplot2grid((3,6),(0,0))
    ax1a.set_title('MLA01 - old growth', fontsize=10)
    ax1a.set_ylabel('height / m',fontsize=axis_size)
    ax1a.annotate('a - 2 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1b = plt.subplot2grid((3,6),(0,1),sharex=ax1a,sharey=ax1a)
    ax1b.annotate('b - 5 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1c = plt.subplot2grid((3,6),(0,2),sharex=ax1a,sharey=ax1a)
    ax1c.annotate('c - 10 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1d = plt.subplot2grid((3,6),(0,3),sharex=ax1a,sharey=ax1a)
    ax1d.annotate('d - 20 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1e = plt.subplot2grid((3,6),(0,4),sharex=ax1a,sharey=ax1a)
    ax1e.annotate('e - 50 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1f = plt.subplot2grid((3,6),(0,5),sharex=ax1a,sharey=ax1a)
    ax1f.annotate('f - 100 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ## E
    ax1g = plt.subplot2grid((3,6),(1,0),sharex=ax1a,sharey=ax1a)
    ax1g.set_title('SAF03 - moderately logged', fontsize=10)
    ax1g.set_ylabel('height / m',fontsize=axis_size)
    ax1g.annotate('g - 2 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1h = plt.subplot2grid((3,6),(1,1),sharex=ax1a,sharey=ax1a)
    ax1h.annotate('h - 5 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1i = plt.subplot2grid((3,6),(1,2),sharex=ax1a,sharey=ax1a)
    ax1i.annotate('i - 10 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1j = plt.subplot2grid((3,6),(1,3),sharex=ax1a,sharey=ax1a)
    ax1j.annotate('j - 20 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1k = plt.subplot2grid((3,6),(1,4),sharex=ax1a,sharey=ax1a)
    ax1k.annotate('k - 50 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1l = plt.subplot2grid((3,6),(1,5),sharex=ax1a,sharey=ax1a)
    ax1l.annotate('l - 100 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ## B North
    ax1m = plt.subplot2grid((3,6),(2,0),sharex=ax1a,sharey=ax1a)
    ax1m.set_title('SAF02 - heavily logged', fontsize=10)
    ax1m.set_ylabel('height / m',fontsize=axis_size)
    ax1m.annotate('m - 2 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1n = plt.subplot2grid((3,6),(2,1),sharex=ax1a,sharey=ax1a)
    ax1n.annotate('n - 5 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1o = plt.subplot2grid((3,6),(2,2),sharex=ax1a,sharey=ax1a)
    ax1o.annotate('o - 10 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1p = plt.subplot2grid((3,6),(2,3),sharex=ax1a,sharey=ax1a)
    ax1p.annotate('p - 20 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1q = plt.subplot2grid((3,6),(2,4),sharex=ax1a,sharey=ax1a)
    ax1q.annotate('q - 50 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    ax1r = plt.subplot2grid((3,6),(2,5),sharex=ax1a,sharey=ax1a)
    ax1r.annotate('r - 100 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',ha='left', va='top', fontsize=10)

    # loop through the subplots
    sp = [ax1a,ax1b,ax1c,ax1d,ax1e,ax1f,ax1g,ax1h,ax1i,ax1j,ax1k,ax1l,ax1m,ax1n,ax1o,ax1p,ax1q,ax1r]
    pkeys = ['2m','5m','10m','20m','50m','100m']
    for pp in range(0,len(sp)):

        if pp<6:
            rad_sp = np.nanmean(PAD_profiles_rad_Belian[pkeys[pp]]['40'],axis=1)
            MH_sp = np.nanmean(PAD_profiles_MH_Belian[pkeys[pp]]['40'],axis=1)
        elif pp<12:
            rad_sp = np.nanmean(PAD_profiles_rad_E[pkeys[pp-6]]['40'],axis=1)
            MH_sp = np.nanmean(PAD_profiles_MH_E[pkeys[pp-6]]['40'],axis=1)
        else:
            rad_sp = np.nanmean(PAD_profiles_rad_BNorth[pkeys[pp-12]]['40'],axis=1)
            MH_sp = np.nanmean(PAD_profiles_MH_BNorth[pkeys[pp-12]]['40'],axis=1)

        rad = np.nanmedian(rad_sp,axis=0)
        radperc = np.nanpercentile(rad_sp,[2.5,97.5],axis=0)
        sp[pp].fill_betweenx(heights[2:],radperc[0][2:],radperc[1][2:],color=colour[1],alpha=0.3)

        MH = np.nanmedian(MH_sp,axis=0)
        MHperc = np.nanpercentile(MH_sp,[2.5,97.5],axis=0)
        sp[pp].fill_betweenx(heights[2:],MHperc[0][2:],MHperc[1][2:],color=colour[0],alpha=0.3)

        if pp==17:
            sp[pp].plot(MH[2:],heights[2:],'-',c=colour[0],linewidth=1.5, label = 'MacArthur-Horn')
            sp[pp].plot(rad[2:],heights[2:],'-',c=colour[1],linewidth=1.5, label = 'multi. return\nrad trans')
            sp[pp].legend(loc=9, bbox_to_anchor=(-0.1, -0.15), ncol=3)
        else:
            sp[pp].plot(MH[2:],heights[2:],'-',c=colour[0],linewidth=1.5)
            sp[pp].plot(rad[2:],heights[2:],'-',c=colour[1],linewidth=1.5)

    yticklabels = ax1b.get_yticklabels() + ax1c.get_yticklabels() + ax1d.get_yticklabels() + ax1e.get_yticklabels() + ax1f.get_yticklabels() + \
                  ax1h.get_yticklabels() + ax1i.get_yticklabels() + ax1j.get_yticklabels() + ax1k.get_yticklabels() + ax1l.get_yticklabels() + \
                  ax1n.get_yticklabels() + ax1o.get_yticklabels() + ax1p.get_yticklabels() + ax1q.get_yticklabels() + ax1r.get_yticklabels()
    xticklabels = ax1a.get_xticklabels() + ax1b.get_xticklabels() + ax1c.get_xticklabels() + ax1d.get_xticklabels() + ax1e.get_xticklabels() + ax1f.get_xticklabels() + \
                  ax1g.get_xticklabels() + ax1h.get_xticklabels() + ax1i.get_xticklabels() + ax1j.get_xticklabels() + ax1k.get_xticklabels() + ax1l.get_xticklabels()

    plt.setp(yticklabels,visible=False)

    ax1a.set_ylim(0,89)
    ax1a.set_xlim(0,0.55)
    for ax in sp:
        ax.locator_params(axis='x',nbins=5)
        ax.yaxis.set_ticks_position('both')

    plt.subplots_adjust(hspace=0.4, wspace = 0.2, bottom = 0.2)

    fig.text(0.5, 0.15, 'PAD / m$^2$m$^{-3}$',fontsize=axis_size, ha='center')

    plt.savefig(figure_name)
    #plt.show()
    return 0

"""
Plot sensitivity of 1ha profiles to point density
"""
def plot_profile_sensitivity_density(figure_number,figure_name,heights,PAD_profiles_MH_Belian,
                            PAD_profiles_MH_BNorth,PAD_profiles_MH_E,PAD_profiles_rad_Belian,
                            PAD_profiles_rad_BNorth,PAD_profiles_rad_E,PAD_profiles_wt_Belian,
                            PAD_profiles_wt_BNorth,PAD_profiles_wt_E):

    fig = plt.figure(figure_number, facecolor='White',figsize=[7,9])
    ax2a = plt.subplot2grid((3,5),(0,0))
    ax2a.set_title('MLA01 - old growth', fontsize=10)
    ax2a.set_ylabel('height / m',fontsize=axis_size)
    ax2a.annotate('a - 5 shots m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=8)

    ax2b = plt.subplot2grid((3,5),(0,1),sharex=ax2a,sharey=ax2a)
    ax2b.annotate('b - 10 shots m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=8)

    ax2c = plt.subplot2grid((3,5),(0,2),sharex=ax2a,sharey=ax2a)
    ax2c.annotate('c - 15 shots m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=8)

    ax2d = plt.subplot2grid((3,5),(0,3),sharex=ax2a,sharey=ax2a)
    ax2d.annotate('d - 20 shots m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=8)

    ax2e = plt.subplot2grid((3,5),(0,4),sharex=ax2a,sharey=ax2a)
    ax2e.annotate('e - 25 shots m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=8)

    ax2f = plt.subplot2grid((3,5),(1,0),sharex=ax2a,sharey=ax2a)
    ax2f.set_title('SAF03 - moderately logged', fontsize=10)
    ax2f.set_ylabel('height / m',fontsize=axis_size)
    ax2f.annotate('f - 5 shots m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=8)

    ax2g = plt.subplot2grid((3,5),(1,1),sharex=ax2a,sharey=ax2a)
    ax2g.annotate('g - 10 shots m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=8)

    ax2h = plt.subplot2grid((3,5),(1,2),sharex=ax2a,sharey=ax2a)
    ax2h.annotate('h - 15 shots m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=8)

    ax2i = plt.subplot2grid((3,5),(1,3),sharex=ax2a,sharey=ax2a)
    ax2i.annotate('i - 20 shots m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=8)

    ax2j = plt.subplot2grid((3,5),(1,4),sharex=ax2a,sharey=ax2a)
    ax2j.annotate('j - 25 shots m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=8)

    ax2k = plt.subplot2grid((3,5),(2,0),sharex=ax2a,sharey=ax2a)
    ax2k.set_title('SAF02 - heavily logged', fontsize=10)
    ax2k.set_ylabel('height / m',fontsize=axis_size)
    ax2k.annotate('k - 5 shots m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=8)

    ax2l = plt.subplot2grid((3,5),(2,1),sharex=ax2a,sharey=ax2a)
    ax2l.annotate('l - 10 shots m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=8)

    ax2m = plt.subplot2grid((3,5),(2,2),sharex=ax2a,sharey=ax2a)
    ax2m.annotate('m - 15 shots m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=8)

    ax2n = plt.subplot2grid((3,5),(2,3),sharex=ax2a,sharey=ax2a)
    ax2n.annotate('n - 20 shots m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=8)

    ax2o = plt.subplot2grid((3,5),(2,4),sharex=ax2a,sharey=ax2a)
    ax2o.annotate('o - 25 shots m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=8)

    # loop through the subplots
    sp = [ax2a,ax2b,ax2c,ax2d,ax2e,ax2f,ax2g,ax2h,ax2i,ax2j,ax2k,ax2l,ax2m,ax2n,ax2o]
    pkeys = ['5','10','15','20','25']
    for pp in range(0,len(sp)):
        if pp<5:
            rad_sp = np.nanmean(PAD_profiles_rad_Belian['20m'][pkeys[pp]],axis=1)
            MH_sp = np.nanmean(PAD_profiles_MH_Belian['20m'][pkeys[pp]],axis=1)
            wt_sp = np.nanmean(PAD_profiles_wt_Belian['20m'][pkeys[pp]],axis=1)
        elif pp<10:
            rad_sp = np.nanmean(PAD_profiles_rad_E['20m'][pkeys[pp-5]],axis=1)
            MH_sp = np.nanmean(PAD_profiles_MH_E['20m'][pkeys[pp-5]],axis=1)
            wt_sp = np.nanmean(PAD_profiles_wt_E['20m'][pkeys[pp-5]],axis=1)
        else:
            rad_sp = np.nanmean(PAD_profiles_rad_BNorth['20m'][pkeys[pp-10]],axis=1)
            MH_sp = np.nanmean(PAD_profiles_MH_BNorth['20m'][pkeys[pp-10]],axis=1)
            wt_sp = np.nanmean(PAD_profiles_wt_BNorth['20m'][pkeys[pp-10]],axis=1)

        rad = np.nanmedian(rad_sp,axis=0)
        radperc = np.nanpercentile(rad_sp,[2.5,97.5],axis=0)
        sp[pp].fill_betweenx(heights[2:],radperc[0][2:],radperc[1][2:],color=colour[1],alpha=0.3)

        MH = np.median(MH_sp,axis=0)
        MHperc = np.percentile(MH_sp,[2.5,97.5],axis=0)
        sp[pp].fill_betweenx(heights[2:],MHperc[0][2:],MHperc[1][2:],color=colour[0],alpha=0.3)

        wt_sp[~np.isfinite(wt_sp)]=np.nan
        wt = np.median(wt_sp,axis=0)
        wtperc = np.percentile(wt_sp,[2.5,97.5],axis=0)
        sp[pp].fill_betweenx(heights[2:],wtperc[0][2:],wtperc[1][2:],color=colour[2],alpha=0.3)

        if pp==14:
            sp[pp].plot(MH[2:],heights[2:],'-',c=colour[0],linewidth=1.5, label = 'M1a')
            sp[pp].plot(wt[2:],heights[2:],'-',c=colour[2],linewidth=1.5, label = 'M1b')
            sp[pp].plot(rad[2:],heights[2:],'-',c=colour[1],linewidth=1.5, label = 'M2')
            sp[pp].legend(loc=9, bbox_to_anchor=(-0.1, -0.15), ncol=3)
        else:
            sp[pp].plot(MH[2:],heights[2:],'-',c=colour[0],linewidth=1.5)
            sp[pp].plot(wt[2:],heights[2:],'-',c=colour[2],linewidth=1.5)
            sp[pp].plot(rad[2:],heights[2:],'-',c=colour[1],linewidth=1.5)

    ax2a.set_ylim(0,89);ax2a.set_xlim(0,0.55)

    yticklabels = ax2b.get_yticklabels() + ax2c.get_yticklabels() + ax2d.get_yticklabels() + ax2e.get_yticklabels() + \
                  ax2g.get_yticklabels() + ax2h.get_yticklabels() + ax2i.get_yticklabels() + ax2j.get_yticklabels() + \
                  ax2l.get_yticklabels() + ax2m.get_yticklabels() + ax2n.get_yticklabels() + ax2o.get_yticklabels()
    xticklabels = ax2a.get_xticklabels() + ax2b.get_xticklabels() + ax2c.get_xticklabels() + ax2d.get_xticklabels() + ax2e.get_xticklabels() + \
                  ax2f.get_xticklabels() + ax2g.get_xticklabels() + ax2h.get_xticklabels() + ax2i.get_xticklabels() + ax2j.get_xticklabels()

    plt.setp(yticklabels,visible=False);plt.setp(xticklabels,visible=False)

    for ax in sp:
        ax.locator_params(axis='x',nbins=5)
        ax.yaxis.set_ticks_position('both')

    plt.subplots_adjust(hspace=0.4, wspace = 0.2, bottom = 0.2)

    fig.text(0.5, 0.15, 'PAD / m$^2$m$^{-3}$',fontsize=axis_size, ha='center')

    plt.savefig(figure_name)
    #plt.show()
    return 0

"""
Plot sensitivity of 1ha profiles to point density
"""
def plot_profile_sensitivity_point_density_old(figure_number,figure_name,heights,PAD_profiles_MH_Belian,
                            PAD_profiles_MH_BNorth,PAD_profiles_MH_E,PAD_profiles_rad_Belian,
                            PAD_profiles_rad_BNorth,PAD_profiles_rad_E):

    fig = plt.figure(figure_number, facecolor='White',figsize=[7,9])
    ax2a = plt.subplot2grid((3,5),(0,0))
    ax2a.set_title('MLA01 - old growth', fontsize=10)
    ax2a.set_ylabel('height / m',fontsize=axis_size)
    ax2a.annotate('a - 5 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=8)

    ax2b = plt.subplot2grid((3,5),(0,1),sharex=ax2a,sharey=ax2a)
    ax2b.annotate('b - 10 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=8)

    ax2c = plt.subplot2grid((3,5),(0,2),sharex=ax2a,sharey=ax2a)
    ax2c.annotate('c - 20 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=8)

    ax2d = plt.subplot2grid((3,5),(0,3),sharex=ax2a,sharey=ax2a)
    ax2d.annotate('d - 30 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=8)

    ax2e = plt.subplot2grid((3,5),(0,4),sharex=ax2a,sharey=ax2a)
    ax2e.annotate('e - 40 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=8)

    ax2f = plt.subplot2grid((3,5),(1,0),sharex=ax2a,sharey=ax2a)
    ax2f.set_title('SAF03 - moderately logged', fontsize=10)
    ax2f.set_ylabel('height / m',fontsize=axis_size)
    ax2f.annotate('f - 5 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=8)

    ax2g = plt.subplot2grid((3,5),(1,1),sharex=ax2a,sharey=ax2a)
    ax2g.annotate('g - 10 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=8)

    ax2h = plt.subplot2grid((3,5),(1,2),sharex=ax2a,sharey=ax2a)
    ax2h.annotate('h - 20 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=8)

    ax2i = plt.subplot2grid((3,5),(1,3),sharex=ax2a,sharey=ax2a)
    ax2i.annotate('i - 30 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=8)

    ax2j = plt.subplot2grid((3,5),(1,4),sharex=ax2a,sharey=ax2a)
    ax2j.annotate('j - 40 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=8)

    ax2k = plt.subplot2grid((3,5),(2,0),sharex=ax2a,sharey=ax2a)
    ax2k.set_title('SAF02 - heavily logged', fontsize=10)
    ax2k.set_ylabel('height / m',fontsize=axis_size)
    ax2k.annotate('k - 5 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=8)

    ax2l = plt.subplot2grid((3,5),(2,1),sharex=ax2a,sharey=ax2a)
    ax2l.annotate('l - 10 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=8)

    ax2m = plt.subplot2grid((3,5),(2,2),sharex=ax2a,sharey=ax2a)
    ax2m.annotate('m - 20 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=8)

    ax2n = plt.subplot2grid((3,5),(2,3),sharex=ax2a,sharey=ax2a)
    ax2n.annotate('n - 30 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=8)

    ax2o = plt.subplot2grid((3,5),(2,4),sharex=ax2a,sharey=ax2a)
    ax2o.annotate('o - 40 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=8)

    # loop through the subplots
    sp = [ax2a,ax2b,ax2c,ax2d,ax2e,ax2f,ax2g,ax2h,ax2i,ax2j,ax2k,ax2l,ax2m,ax2n,ax2o]
    pkeys = ['5','10','20','30','40']
    for pp in range(0,len(sp)):
        if pp<5:
            rad_sp = np.nanmean(PAD_profiles_rad_Belian['20m'][pkeys[pp]],axis=1)
            MH_sp = np.nanmean(PAD_profiles_MH_Belian['20m'][pkeys[pp]],axis=1)
        elif pp<10:
            rad_sp = np.nanmean(PAD_profiles_rad_E['20m'][pkeys[pp-5]],axis=1)
            MH_sp = np.nanmean(PAD_profiles_MH_E['20m'][pkeys[pp-5]],axis=1)
        else:
            rad_sp = np.nanmean(PAD_profiles_rad_BNorth['20m'][pkeys[pp-10]],axis=1)
            MH_sp = np.nanmean(PAD_profiles_MH_BNorth['20m'][pkeys[pp-10]],axis=1)

        rad = np.nanmedian(rad_sp,axis=0)
        radperc = np.nanpercentile(rad_sp,[2.5,97.5],axis=0)
        sp[pp].fill_betweenx(heights[2:],radperc[0][2:],radperc[1][2:],color=colour[1],alpha=0.3)

        MH = np.median(MH_sp,axis=0)
        MHperc = np.percentile(MH_sp,[2.5,97.5],axis=0)
        sp[pp].fill_betweenx(heights[2:],MHperc[0][2:],MHperc[1][2:],color=colour[0],alpha=0.3)

        if pp==14:
            sp[pp].plot(MH[2:],heights[2:],'-',c=colour[0],linewidth=1.5, label = 'MacArthur-Horn')
            sp[pp].plot(rad[2:],heights[2:],'-',c=colour[1],linewidth=1.5, label = 'multi return\nrad trans')
            sp[pp].legend(loc=9, bbox_to_anchor=(-0.1, -0.15), ncol=3)#(loc='upper right')
        else:
            sp[pp].plot(MH[2:],heights[2:],'-',c=colour[0],linewidth=1.5)
            sp[pp].plot(rad[2:],heights[2:],'-',c=colour[1],linewidth=1.5)

    ax2a.set_ylim(0,89);ax2a.set_xlim(0,0.55)

    yticklabels = ax2b.get_yticklabels() + ax2c.get_yticklabels() + ax2d.get_yticklabels() + ax2e.get_yticklabels() + \
                  ax2g.get_yticklabels() + ax2h.get_yticklabels() + ax2i.get_yticklabels() + ax2j.get_yticklabels() + \
                  ax2l.get_yticklabels() + ax2m.get_yticklabels() + ax2n.get_yticklabels() + ax2o.get_yticklabels()
    xticklabels = ax2a.get_xticklabels() + ax2b.get_xticklabels() + ax2c.get_xticklabels() + ax2d.get_xticklabels() + ax2e.get_xticklabels() + \
                  ax2f.get_xticklabels() + ax2g.get_xticklabels() + ax2h.get_xticklabels() + ax2i.get_xticklabels() + ax2j.get_xticklabels()

    plt.setp(yticklabels,visible=False);plt.setp(xticklabels,visible=False)

    for ax in sp:
        ax.locator_params(axis='x',nbins=5)
        ax.yaxis.set_ticks_position('both')

    plt.subplots_adjust(hspace=0.4, wspace = 0.2, bottom = 0.2)

    fig.text(0.5, 0.15, 'PAD / m$^2$m$^{-3}$',fontsize=axis_size, ha='center')

    plt.savefig(figure_name)
    #plt.show()
    return 0

"""
Plot fraction of unsampled layers in profile
"""
def plot_penetration_limits(figure_number,figure_name,heights,penetration_lim_Belian,
                            penetration_lim_E,penetration_lim_BNorth):
    ## Belian
    fig =plt.figure(figure_number, facecolor='White',figsize=[7,3])
    ax1a = plt.subplot2grid((1,3),(0,0))
    ax1a.set_ylabel('height / m',fontsize=axis_size)
    ax1a.annotate('a - MLA01\n(old growth)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    #ax1a.set_xlabel('occluded voxels / %',fontsize=axis_size)

    ax1b = plt.subplot2grid((1,3),(0,1),sharex=ax1a,sharey=ax1a)
    ax1b.annotate('b - SAF03\n(moderately logged)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1b.set_xlabel('percentage of voxels occluded voxels below\npenetration limit of first returns / %',fontsize=axis_size)

    ax1c = plt.subplot2grid((1,3),(0,2),sharex=ax1a,sharey=ax1a)
    ax1c.annotate('c - SAF02\n(heavily logged)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    #ax1c.set_xlabel('occluded voxels / %',fontsize=axis_size)

    # loop through the subplots
    sp = [ax1a,ax1b,ax1c]
    pkeys = ['5m','10m','20m']
    palette_colour = ['#E0007F','#1A2BCE','#46E900']
    for pp in range(0,len(pkeys)):

        lim_OG = np.nanmean(penetration_lim_Belian[pkeys[pp]]['25'],axis=1)*100.
        lim_ML = np.nanmean(penetration_lim_E[pkeys[pp]]['25'],axis=1)*100.
        lim_HL = np.nanmean(penetration_lim_BNorth[pkeys[pp]]['25'],axis=1)*100.

        limperc_OG = np.nanpercentile(lim_OG,[2.5,50,97.5],axis=0)
        limperc_ML = np.nanpercentile(lim_ML,[2.5,50,97.5],axis=0)
        limperc_HL = np.nanpercentile(lim_HL,[2.5,50,97.5],axis=0)

        ax1a.fill_betweenx(heights[2:],limperc_OG[0][2:],limperc_OG[2][2:],color=palette_colour[pp],alpha=0.3)
        ax1a.plot(limperc_OG[1][2:],heights[2:],'-',c=palette_colour[pp],linewidth=1.5)
        ax1b.fill_betweenx(heights[2:],limperc_ML[0][2:],limperc_ML[2][2:],color=palette_colour[pp],alpha=0.3)
        ax1b.plot(limperc_ML[1][2:],heights[2:],'-',c=palette_colour[pp],linewidth=1.5)
        ax1c.fill_betweenx(heights[2:],limperc_HL[0][2:],limperc_HL[2][2:],color=palette_colour[pp],alpha=0.3)
        ax1c.plot(limperc_HL[1][2:],heights[2:],'-',c=palette_colour[pp],linewidth=1.5,label=pkeys[pp])

    ax1a.set_ylim(0,89)
    ax1a.set_xlim(0,100)

    yticklabels = ax1b.get_yticklabels() + ax1c.get_yticklabels()
    plt.setp(yticklabels,visible=False)

    for ax in sp:
        ax.locator_params(axis='x',nbins=5)
        ax.yaxis.set_ticks_position('both')

    plt.subplots_adjust(hspace=0.4, wspace = 0.3,bottom=0.2)

    ax1c.legend(loc=5,title='grid resolution')
    plt.savefig(figure_name)
    #plt.show()
    return 0

"""
Plot individual profile sensitivity to resolution
"""
def plot_profile_sensitivity_to_resolution_individual_CI(figure_number,figure_name,heights,
                            PAD_profiles_MH,PAD_profiles_rad,PAD_profiles_wt):
    fig = plt.figure(figure_number, facecolor='White',figsize=[8,4])
    ax3a = plt.subplot2grid((1,5),(0,0))
    ax3a.set_ylabel('height / m',fontsize=axis_size)
    ax3a.annotate('a - 5 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

    ax3b = plt.subplot2grid((1,5),(0,1),sharex=ax3a,sharey=ax3a)
    ax3b.annotate('b - 10 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

    ax3c = plt.subplot2grid((1,5),(0,2),sharex=ax3a,sharey=ax3a)
    ax3c.annotate('c - 20 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

    ax3d = plt.subplot2grid((1,5),(0,3),sharex=ax3a,sharey=ax3a)
    ax3d.annotate('d - 50 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

    ax3e = plt.subplot2grid((1,5),(0,4),sharex=ax3a,sharey=ax3a)
    ax3e.annotate('e - 100 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

    # loop through the subplots
    sp = [ax3a,ax3b,ax3c,ax3d,ax3e]
    pkeys = ['5m','10m','20m','50m','100m']
    for pp in range(0,len(sp)):
        MHl = np.nanpercentile(PAD_profiles_MH[pkeys[pp]]['25'],2.5,axis=0)
        MH25 = np.nanpercentile(PAD_profiles_MH[pkeys[pp]]['25'],25,axis=0)
        MH75 = np.nanpercentile(PAD_profiles_MH[pkeys[pp]]['25'],75,axis=0)
        MHu = np.nanpercentile(PAD_profiles_MH[pkeys[pp]]['25'],97.5,axis=0)
        MH_95CI_i = ((MHu-MHl)/np.nanmean(PAD_profiles_MH[pkeys[pp]]['25'],axis=0))*100.
        MH_50CI_i = ((MH75-MH25)/np.nanmean(PAD_profiles_MH[pkeys[pp]]['25'],axis=0))*100.
        MH_95CI = np.nansum(MH_95CI_i,axis=0)/np.sum(np.isfinite(MH_95CI_i),axis=0)
        MH_50CI = np.nansum(MH_50CI_i,axis=0)/np.sum(np.isfinite(MH_50CI_i),axis=0)
        sp[pp].plot(MH_95CI[2:],heights[2:],':',c=colour[0],linewidth=1)
        sp[pp].plot(MH_50CI[2:],heights[2:],'-',c=colour[0],linewidth=1)
        if pp==2:
            sp[pp].plot(MH_50CI[2:],heights[2:],'-',c=colour[0],linewidth=1,label='M1a 50% CI')
            sp[pp].plot(MH_95CI[2:],heights[2:],':',c=colour[0],linewidth=1,label='M1a 95% CI')
        else:
            sp[pp].plot(MH_95CI[2:],heights[2:],':',c=colour[0],linewidth=1)
            sp[pp].plot(MH_50CI[2:],heights[2:],'-',c=colour[0],linewidth=1)

        PAD_profiles_wt[pkeys[pp]]['25'][~np.isfinite(PAD_profiles_wt[pkeys[pp]]['25'])]=np.nan
        wtl = np.nanpercentile(PAD_profiles_wt[pkeys[pp]]['25'],2.5,axis=0)
        wt25 = np.nanpercentile(PAD_profiles_wt[pkeys[pp]]['25'],25,axis=0)
        wt75 = np.nanpercentile(PAD_profiles_wt[pkeys[pp]]['25'],75,axis=0)
        wtu = np.nanpercentile(PAD_profiles_wt[pkeys[pp]]['25'],97.5,axis=0)
        wt_95CI_i = ((wtu-wtl)/np.nanmean(PAD_profiles_wt[pkeys[pp]]['25'],axis=0))*100.
        wt_50CI_i = ((wt75-wt25)/np.nanmean(PAD_profiles_wt[pkeys[pp]]['25'],axis=0))*100.
        wt_95CI = np.nansum(wt_95CI_i,axis=0)/np.sum(np.isfinite(wt_95CI_i),axis=0)
        wt_50CI = np.nansum(wt_50CI_i,axis=0)/np.sum(np.isfinite(wt_50CI_i),axis=0)
        sp[pp].plot(wt_95CI[2:],heights[2:],':',c=colour[2],linewidth=1)
        sp[pp].plot(wt_50CI[2:],heights[2:],'-',c=colour[2],linewidth=1)
        if pp==2:
            sp[pp].plot(wt_50CI[2:],heights[2:],'-',c=colour[2],linewidth=1,label='M1b 50% CI')
            sp[pp].plot(wt_95CI[2:],heights[2:],':',c=colour[2],linewidth=1,label='M1b 95% CI')
        else:
            sp[pp].plot(wt_95CI[2:],heights[2:],':',c=colour[2],linewidth=1)
            sp[pp].plot(wt_50CI[2:],heights[2:],'-',c=colour[2],linewidth=1)

        if pkeys[pp]!= '5m':
            radl = np.nanpercentile(PAD_profiles_rad[pkeys[pp]]['25'],2.5,axis=0)
            rad25 = np.nanpercentile(PAD_profiles_rad[pkeys[pp]]['25'],25,axis=0)
            rad75 = np.nanpercentile(PAD_profiles_rad[pkeys[pp]]['25'],75,axis=0)
            radu = np.nanpercentile(PAD_profiles_rad[pkeys[pp]]['25'],97.5,axis=0)
            rad_95CI_i = ((radu-radl)/np.nanmean(PAD_profiles_rad[pkeys[pp]]['25'],axis=0))*100.
            rad_50CI_i = ((rad75-rad25)/np.nanmean(PAD_profiles_rad[pkeys[pp]]['25'],axis=0))*100.
            rad_95CI = np.nansum(rad_95CI_i,axis=0)/np.sum(np.isfinite(rad_95CI_i),axis=0)
            rad_50CI = np.nansum(rad_50CI_i,axis=0)/np.sum(np.isfinite(rad_50CI_i),axis=0)
            if pp==2:
                sp[pp].plot(rad_50CI[2:],heights[2:],'-',c=colour[1],linewidth=1,label='M2 50% CI')
                sp[pp].plot(rad_95CI[2:],heights[2:],':',c=colour[1],linewidth=1,label='M2 95% CI')
                lgd = sp[pp].legend(loc=9, bbox_to_anchor=(0.5, -0.3), ncol=3)
            else:
                sp[pp].plot(rad_95CI[2:],heights[2:],':',c=colour[1],linewidth=1)
                sp[pp].plot(rad_50CI[2:],heights[2:],'-',c=colour[1],linewidth=1)

    yticklabels = ax3b.get_yticklabels() + ax3c.get_yticklabels() + ax3d.get_yticklabels() + ax3e.get_yticklabels()
    plt.setp(yticklabels,visible=False)

    ax3a.set_ylim(0,89)
    ax3a.set_xlim(0,520)

    for ax in sp:
        ax.locator_params(axis='x',nbins=5)
        ax.yaxis.set_ticks_position('both')

    plt.subplots_adjust(hspace=0.4, wspace = 0.2, bottom = 0.4)

    fig.text(0.5, 0.3, "relative CIs for individual profiles / %",
                    fontsize=axis_size, ha='center')

    plt.savefig(figure_name)
    #plt.show()
    return 0

"""
Plot individual profile sensitivity to resolution
"""
def plot_profile_sensitivity_to_point_density_individual_CI(figure_number,figure_name,heights,
                            PAD_profiles_MH,PAD_profiles_rad, PAD_profiles_wt):
    plt.figure(figure_number, facecolor='White',figsize=[8,4])
    axa = plt.subplot2grid((1,5),(0,0))
    axa.set_ylabel('height / m',fontsize=axis_size)
    axa.annotate('a - 5 shots m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',
                    backgroundcolor='none',horizontalalignment='left',
                    verticalalignment='top', fontsize=8)
    axb = plt.subplot2grid((1,5),(0,1),sharex=axa,sharey=axa)
    axb.annotate('b - 10 shots m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',
                    backgroundcolor='none',horizontalalignment='left',
                    verticalalignment='top', fontsize=8)
    axc = plt.subplot2grid((1,5),(0,2),sharex=axa,sharey=axa)
    axc.annotate('c - 15 shots m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',
                    backgroundcolor='none',horizontalalignment='left',
                    verticalalignment='top', fontsize=8)
    axc.set_xlabel("relative CIs for individual profiles / %",fontsize = axis_size)
    axd = plt.subplot2grid((1,5),(0,3),sharex=axa,sharey=axa)
    axd.annotate('d - 20 shots m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',
                    backgroundcolor='none',horizontalalignment='left',
                    verticalalignment='top', fontsize=8)
    axe = plt.subplot2grid((1,5),(0,4),sharex=axa,sharey=axa)
    axe.annotate('e - 25 shots m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',
                    backgroundcolor='none',horizontalalignment='left',
                    verticalalignment='top', fontsize=8)
    # loop through the subplots
    sp = [axa,axb,axc,axd,axe]
    pkeys = ['5','10','15','20','25']
    for pp in range(0,len(sp)):
        MHl = np.nanpercentile(PAD_profiles_MH['20m'][pkeys[pp]],2.5,axis=0)
        MH25 = np.nanpercentile(PAD_profiles_MH['20m'][pkeys[pp]],25,axis=0)
        MH75 = np.nanpercentile(PAD_profiles_MH['20m'][pkeys[pp]],75,axis=0)
        MHu = np.nanpercentile(PAD_profiles_MH['20m'][pkeys[pp]],97.5,axis=0)
        MH_95CI_i = ((MHu-MHl)/np.nanmean(PAD_profiles_MH['20m'][pkeys[pp]],axis=0))*100.
        MH_50CI_i = ((MH75-MH25)/np.nanmean(PAD_profiles_MH['20m'][pkeys[pp]],axis=0))*100.
        MH_95CI = np.nansum(MH_95CI_i,axis=0)/np.sum(np.isfinite(MH_95CI_i),axis=0)
        MH_50CI = np.nansum(MH_50CI_i,axis=0)/np.sum(np.isfinite(MH_50CI_i),axis=0)
        sp[pp].plot(MH_95CI[2:],heights[2:],':',c=colour[0],linewidth=1)
        sp[pp].plot(MH_50CI[2:],heights[2:],'-',c=colour[0],linewidth=1)
        if pp==2:
            sp[pp].plot(MH_50CI[2:],heights[2:],'-',c=colour[0],linewidth=1,label='M1a 50% CI')
            sp[pp].plot(MH_95CI[2:],heights[2:],':',c=colour[0],linewidth=1,label='M1a 95% CI')
        else:
            sp[pp].plot(MH_95CI[2:],heights[2:],':',c=colour[0],linewidth=1)
            sp[pp].plot(MH_50CI[2:],heights[2:],'-',c=colour[0],linewidth=1)

        PAD_profiles_wt['20m'][pkeys[pp]][~np.isfinite(PAD_profiles_wt['20m'][pkeys[pp]])]=np.nan
        wtl = np.nanpercentile(PAD_profiles_wt['20m'][pkeys[pp]],2.5,axis=0)
        wt25 = np.nanpercentile(PAD_profiles_wt['20m'][pkeys[pp]],25,axis=0)
        wt75 = np.nanpercentile(PAD_profiles_wt['20m'][pkeys[pp]],75,axis=0)
        wtu = np.nanpercentile(PAD_profiles_wt['20m'][pkeys[pp]],97.5,axis=0)
        wt_95CI_i = ((wtu-wtl)/np.nanmean(PAD_profiles_wt['20m'][pkeys[pp]],axis=0))*100.
        wt_50CI_i = ((wt75-wt25)/np.nanmean(PAD_profiles_wt['20m'][pkeys[pp]],axis=0))*100.
        wt_95CI = np.nansum(wt_95CI_i,axis=0)/np.sum(np.isfinite(wt_95CI_i),axis=0)
        wt_50CI = np.nansum(wt_50CI_i,axis=0)/np.sum(np.isfinite(wt_50CI_i),axis=0)
        sp[pp].plot(wt_95CI[2:],heights[2:],':',c=colour[2],linewidth=1)
        sp[pp].plot(wt_50CI[2:],heights[2:],'-',c=colour[2],linewidth=1)
        if pp==2:
            sp[pp].plot(wt_50CI[2:],heights[2:],'-',c=colour[2],linewidth=1,label='M1b 50% CI')
            sp[pp].plot(wt_95CI[2:],heights[2:],':',c=colour[2],linewidth=1,label='M1b 95% CI')
        else:
            sp[pp].plot(wt_95CI[2:],heights[2:],':',c=colour[2],linewidth=1)
            sp[pp].plot(wt_50CI[2:],heights[2:],'-',c=colour[2],linewidth=1)


        radl = np.nanpercentile(PAD_profiles_rad['20m'][pkeys[pp]],2.5,axis=0)
        rad25 = np.nanpercentile(PAD_profiles_rad['20m'][pkeys[pp]],25,axis=0)
        rad75 = np.nanpercentile(PAD_profiles_rad['20m'][pkeys[pp]],75,axis=0)
        radu = np.nanpercentile(PAD_profiles_rad['20m'][pkeys[pp]],97.5,axis=0)
        rad_95CI_i = ((radu-radl)/np.nanmean(PAD_profiles_rad['20m'][pkeys[pp]],axis=0))*100.
        rad_50CI_i = ((rad75-rad25)/np.nanmean(PAD_profiles_rad['20m'][pkeys[pp]],axis=0))*100.
        rad_95CI = np.nansum(rad_95CI_i,axis=0)/np.sum(np.isfinite(rad_95CI_i),axis=0)
        rad_50CI = np.nansum(rad_50CI_i,axis=0)/np.sum(np.isfinite(rad_50CI_i),axis=0)
        if pp==2:
            sp[pp].plot(rad_50CI[2:],heights[2:],'-',c=colour[1],linewidth=1,label='M2 50% CI')
            sp[pp].plot(rad_95CI[2:],heights[2:],':',c=colour[1],linewidth=1,label='M2 95% CI')
            lgd = sp[pp].legend(loc=9, bbox_to_anchor=(0.5, -0.3), ncol=3)
        else:
            sp[pp].plot(rad_95CI[2:],heights[2:],':',c=colour[1],linewidth=1)
            sp[pp].plot(rad_50CI[2:],heights[2:],'-',c=colour[1],linewidth=1)

    yticklabels = axb.get_yticklabels() + axc.get_yticklabels() + axd.get_yticklabels() + axe.get_yticklabels()
    plt.setp(yticklabels,visible=False)

    axa.set_ylim(0,89)
    axa.set_xlim(0,520)

    for ax in sp:
        ax.locator_params(axis='x',nbins=5)
        ax.yaxis.set_ticks_position('both')

    plt.subplots_adjust(hspace=0.4, wspace = 0.2, bottom = 0.3)
    plt.savefig(figure_name)
    #plt.show()
    return 0

"""
Plot sensitivity of PAI estimates
"""
def plot_PAI_sensitivity(figure_number,figure_name,
                    res_profiles_MH_Belian,res_profiles_MH_BNorth,res_profiles_MH_E,
                    res_profiles_rad_Belian,res_profiles_rad_BNorth,res_profiles_rad_E,
                    res_profiles_wt_Belian,res_profiles_wt_BNorth,res_profiles_wt_E,
                    dens_profiles_MH_Belian,dens_profiles_MH_BNorth,dens_profiles_MH_E,
                    dens_profiles_rad_Belian,dens_profiles_rad_BNorth,dens_profiles_rad_E,
                    dens_profiles_wt_Belian,dens_profiles_wt_BNorth,dens_profiles_wt_E):

    res_keys = ['5m','10m','20m','50m','100m']
    N_res = len(res_keys)
    dens_keys = ['25']
    N_dens = len(dens_keys)
    dens = []
    res = []
    PAI_MH = np.arange(0)
    PAI_rad = np.arange(0)
    PAI_wt = np.arange(0)
    plot = []
    plot_name = ['MLA01','SAF03','SAF02']
    for rr in range(0,N_res):
        skip_rad=False
        if res_keys[rr]=='5m':
            skip_rad=True
        for dd in range(0,N_dens):
            for pp in range(0,len(plot_name)):
                res_profiles_wt_Belian[res_keys[rr]][dens_keys[dd]][~np.isfinite(res_profiles_wt_Belian[res_keys[rr]][dens_keys[dd]])]=np.nan
                res_profiles_wt_BNorth[res_keys[rr]][dens_keys[dd]][~np.isfinite(res_profiles_wt_BNorth[res_keys[rr]][dens_keys[dd]])]=np.nan
                res_profiles_wt_E[res_keys[rr]][dens_keys[dd]][~np.isfinite(res_profiles_wt_E[res_keys[rr]][dens_keys[dd]])]=np.nan
                if pp==0:
                    N_obs = np.sum(np.nanmean(res_profiles_MH_Belian[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1).size
                    PAI_wt = np.append(PAI_wt,np.sum(np.nanmean(res_profiles_wt_Belian[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    PAI_MH = np.append(PAI_MH,np.sum(np.nanmean(res_profiles_MH_Belian[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    if skip_rad:
                        PAI_rad = np.append(PAI_rad,np.zeros(N_obs)*np.nan)
                    else:
                        PAI_rad = np.append(PAI_rad,np.sum(np.nanmean(res_profiles_rad_Belian[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    for ii in range(0,N_obs):
                        dens.append(dens_keys[dd])
                        res.append(res_keys[rr])
                        plot.append(plot_name[pp])
                if pp==1:
                    N_obs = np.sum(np.nanmean(res_profiles_MH_E[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1).size
                    PAI_wt = np.append(PAI_wt,np.sum(np.nanmean(res_profiles_wt_E[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    PAI_MH = np.append(PAI_MH,np.sum(np.nanmean(res_profiles_MH_E[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    if skip_rad:
                        PAI_rad = np.append(PAI_rad,np.zeros(N_obs)*np.nan)
                    else:
                        PAI_rad = np.append(PAI_rad,np.sum(np.nanmean(res_profiles_rad_E[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    for ii in range(0,N_obs):
                        dens.append(dens_keys[dd])
                        res.append(res_keys[rr])
                        plot.append(plot_name[pp])
                if pp==2:
                    N_obs = np.sum(np.nanmean(res_profiles_MH_BNorth[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1).size
                    PAI_wt = np.append(PAI_wt,np.sum(np.nanmean(res_profiles_wt_BNorth[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    PAI_MH = np.append(PAI_MH,np.sum(np.nanmean(res_profiles_MH_BNorth[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    if skip_rad:
                        PAI_rad = np.append(PAI_rad,np.zeros(N_obs)*np.nan)
                    else:
                        PAI_rad = np.append(PAI_rad,np.sum(np.nanmean(res_profiles_rad_BNorth[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))

                    for ii in range(0,N_obs):
                        dens.append(dens_keys[dd])
                        res.append(res_keys[rr])
                        plot.append(plot_name[pp])

    df1 = pd.DataFrame({'plot' : plot,'point density' : dens,'resolution' : res,'PAI_MH':PAI_MH,'PAI_wt':PAI_wt,'PAI_rad':PAI_rad})

    pp=0
    fig = plt.figure(figure_number, facecolor='White',figsize=(8,6))
    axa = plt.subplot2grid((2,3),(0,0))
    axa.set_title('a - M1a', fontsize=10)
    axa.set_ylabel('PAI',fontsize=axis_size)
    sns.violinplot(x='resolution',y='PAI_MH',data=df1,linewidth=0.5,inner='box',
                    scale='width',hue="plot",palette = colour,dodge=False,saturation=0.7,cut=0)
    axa.set_xlabel('grid resolution')
    axa.set_ylabel('PAI',fontsize=axis_size)

    axb = plt.subplot2grid((2,3),(0,1),sharex=axa,sharey=axa)
    axb.set_title('b - M1b', fontsize=10)
    axb.set_ylabel('PAI',fontsize=axis_size)
    sns.violinplot(x='resolution',y='PAI_wt',data=df1,linewidth=0.5,inner='box',
                    scale='width',hue="plot",palette = colour,dodge=False,saturation=0.7,cut=0)
    axb.set_xlabel('grid resolution')
    axb.set_ylabel('')

    axc = plt.subplot2grid((2,3),(0,2),sharex=axa,sharey=axa)
    axc.set_title('c - M2', fontsize=10)
    sns.violinplot(x='resolution',y='PAI_rad',data=df1,linewidth=0.5,inner='box',
                    scale='width',hue="plot",palette = colour,dodge=False,saturation=0.7,cut=0)
    axc.set_ylabel('')
    axc.set_xlabel('grid resolution')

    # --------------------------------------------------------------------------
    # Now plot pulse density sensitivity
    res_keys = ['20m']; N_res = len(res_keys)
    dens_keys = ['5','10','15','20','25']; N_dens = len(dens_keys)
    dens = []; res = []; plot = []
    PAI_MH = np.arange(0); PAI_rad = np.arange(0); PAI_wt = np.arange(0)
    plot_name = ['MLA01','SAF03','SAF02']
    for rr in range(0,N_res):
        for dd in range(0,N_dens):
            for pp in range(0,len(plot_name)):
                dens_profiles_wt_Belian[res_keys[rr]][dens_keys[dd]][~np.isfinite(dens_profiles_wt_Belian[res_keys[rr]][dens_keys[dd]])]=np.nan
                dens_profiles_wt_BNorth[res_keys[rr]][dens_keys[dd]][~np.isfinite(dens_profiles_wt_BNorth[res_keys[rr]][dens_keys[dd]])]=np.nan
                dens_profiles_wt_E[res_keys[rr]][dens_keys[dd]][~np.isfinite(dens_profiles_wt_E[res_keys[rr]][dens_keys[dd]])]=np.nan
                if pp==0:
                    PAI_MH = np.append(PAI_MH,np.sum(np.nanmean(dens_profiles_MH_Belian[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    PAI_wt = np.append(PAI_wt,np.sum(np.nanmean(dens_profiles_wt_Belian[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    PAI_rad = np.append(PAI_rad,np.sum(np.nanmean(dens_profiles_rad_Belian[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    N_obs = np.sum(np.nanmean(dens_profiles_rad_Belian[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1).size
                    for ii in range(0,N_obs):
                        dens.append(dens_keys[dd])
                        res.append(res_keys[rr])
                        plot.append(plot_name[pp])
                if pp==1:
                    PAI_wt = np.append(PAI_wt,np.sum(np.nanmean(dens_profiles_wt_E[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    PAI_MH = np.append(PAI_MH,np.sum(np.nanmean(dens_profiles_MH_E[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    PAI_rad = np.append(PAI_rad,np.sum(np.nanmean(dens_profiles_rad_E[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    N_obs = np.sum(np.nanmean(dens_profiles_rad_E[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1).size
                    for ii in range(0,N_obs):
                        dens.append(dens_keys[dd])
                        res.append(res_keys[rr])
                        plot.append(plot_name[pp])
                if pp==2:
                    PAI_wt = np.append(PAI_wt,np.sum(np.nanmean(dens_profiles_wt_BNorth[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    PAI_MH = np.append(PAI_MH,np.sum(np.nanmean(dens_profiles_MH_BNorth[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    PAI_rad = np.append(PAI_rad,np.sum(np.nanmean(dens_profiles_rad_BNorth[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    N_obs = np.sum(np.nanmean(dens_profiles_rad_BNorth[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1).size
                    for ii in range(0,N_obs):
                        dens.append(dens_keys[dd])
                        res.append(res_keys[rr])
                        plot.append(plot_name[pp])

    df2 = pd.DataFrame({'plot' : plot,'density' : dens,'resolution' : res,'PAI_MH':PAI_MH,'PAI_wt':PAI_wt,'PAI_rad':PAI_rad})
    df2['density'][df2['density']=='5']='05'
    pp=0
    axd = plt.subplot2grid((2,3),(1,0),sharey=axa)
    axd.set_title('d - M1a', fontsize=10)
    axd.set_ylabel('PAI',fontsize=axis_size)
    sns.violinplot(x='density',y='PAI_MH',data=df2,linewidth=0.5,inner='box',
                    scale='width',hue="plot",palette = colour,dodge=False,saturation=0.7,cut=0)
    axd.set_xlabel('pulse density / pts m$^{-2}$')
    axd.set_ylabel('PAI',fontsize=axis_size)

    axe = plt.subplot2grid((2,3),(1,1),sharex=axd,sharey=axa)
    axe.set_title('e - M1b', fontsize=10)
    axe.set_ylabel('PAI',fontsize=axis_size)
    sns.violinplot(x='density',y='PAI_wt',data=df2,linewidth=0.5,inner='box',
                    scale='width',hue="plot",palette = colour,dodge=False,saturation=0.7,cut=0)
    axe.set_ylabel('')

    axe.set_xlabel('pulse density / pts m$^{-2}$')

    axf = plt.subplot2grid((2,3),(1,2),sharex=axd,sharey=axa)
    axf.set_title('f - M2', fontsize=10)
    sns.violinplot(x='density',y='PAI_rad',data=df2,linewidth=0.5,inner='box',
                    scale='width',hue="plot",palette = colour,dodge=False,saturation=0.7,cut=0)
    axf.set_ylabel('')
    axf.set_xlabel('pulse density / pts m$^{-2}$')

    axa.set_ylim(top=13)

    axes =[axa,axb,axc,axd,axe,axf]
    plt.setp(axb.get_yticklabels(),visible=False);plt.setp(axc.get_yticklabels(),visible=False)
    plt.setp(axe.get_yticklabels(),visible=False);plt.setp(axf.get_yticklabels(),visible=False)

    for ax in axes:
        ax.yaxis.set_ticks_position('both')
        if ax!=axf:
            ax.legend_.remove()
        else:
            axf.legend(loc =9, bbox_to_anchor=(-0.7, -0.3),ncol=3)
    plt.subplots_adjust(wspace = 0.2,hspace=0.4,bottom=0.2)
    plt.savefig(figure_name)

    return 0


def plot_PAI_sensitivity_resolution(figure_number,figure_name,PAD_profiles_MH_Belian,PAD_profiles_MH_BNorth,PAD_profiles_MH_E,
                    PAD_profiles_rad_Belian,PAD_profiles_rad_BNorth,PAD_profiles_rad_E,
                    PAD_profiles_wt_Belian,PAD_profiles_wt_BNorth,PAD_profiles_wt_E):

    res_keys = ['5m','10m','20m','50m','100m']
    N_res = len(res_keys)
    dens_keys = ['25']
    N_dens = len(dens_keys)
    dens = []
    res = []
    PAI_MH = np.arange(0)
    PAI_rad = np.arange(0)
    PAI_wt = np.arange(0)
    plot = []
    plot_name = ['MLA01','SAF03','SAF02']
    for rr in range(0,N_res):
        skip_rad=False
        if res_keys[rr]=='5m':
            skip_rad=True
        for dd in range(0,N_dens):
            for pp in range(0,len(plot_name)):
                PAD_profiles_wt_Belian[res_keys[rr]][dens_keys[dd]][~np.isfinite(PAD_profiles_wt_Belian[res_keys[rr]][dens_keys[dd]])]=np.nan
                PAD_profiles_wt_BNorth[res_keys[rr]][dens_keys[dd]][~np.isfinite(PAD_profiles_wt_BNorth[res_keys[rr]][dens_keys[dd]])]=np.nan
                PAD_profiles_wt_E[res_keys[rr]][dens_keys[dd]][~np.isfinite(PAD_profiles_wt_E[res_keys[rr]][dens_keys[dd]])]=np.nan
                if pp==0:
                    N_obs = np.sum(np.nanmean(PAD_profiles_MH_Belian[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1).size
                    PAI_wt = np.append(PAI_wt,np.sum(np.nanmean(PAD_profiles_wt_Belian[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    PAI_MH = np.append(PAI_MH,np.sum(np.nanmean(PAD_profiles_MH_Belian[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    if skip_rad:
                        PAI_rad = np.append(PAI_rad,np.zeros(N_obs)*np.nan)
                    else:
                        PAI_rad = np.append(PAI_rad,np.sum(np.nanmean(PAD_profiles_rad_Belian[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    for ii in range(0,N_obs):
                        dens.append(dens_keys[dd])
                        res.append(res_keys[rr])
                        plot.append(plot_name[pp])
                if pp==1:
                    N_obs = np.sum(np.nanmean(PAD_profiles_MH_E[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1).size
                    PAI_wt = np.append(PAI_wt,np.sum(np.nanmean(PAD_profiles_wt_E[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    PAI_MH = np.append(PAI_MH,np.sum(np.nanmean(PAD_profiles_MH_E[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    if skip_rad:
                        PAI_rad = np.append(PAI_rad,np.zeros(N_obs)*np.nan)
                    else:
                        PAI_rad = np.append(PAI_rad,np.sum(np.nanmean(PAD_profiles_rad_E[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    for ii in range(0,N_obs):
                        dens.append(dens_keys[dd])
                        res.append(res_keys[rr])
                        plot.append(plot_name[pp])
                if pp==2:
                    N_obs = np.sum(np.nanmean(PAD_profiles_MH_BNorth[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1).size
                    PAI_wt = np.append(PAI_wt,np.sum(np.nanmean(PAD_profiles_wt_BNorth[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    PAI_MH = np.append(PAI_MH,np.sum(np.nanmean(PAD_profiles_MH_BNorth[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    if skip_rad:
                        PAI_rad = np.append(PAI_rad,np.zeros(N_obs)*np.nan)
                    else:
                        PAI_rad = np.append(PAI_rad,np.sum(np.nanmean(PAD_profiles_rad_BNorth[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))

                    for ii in range(0,N_obs):
                        dens.append(dens_keys[dd])
                        res.append(res_keys[rr])
                        plot.append(plot_name[pp])

    df = pd.DataFrame({'plot' : plot,'point density' : dens,'resolution' : res,'PAI_MH':PAI_MH,'PAI_wt':PAI_wt,'PAI_rad':PAI_rad})

    pp=0
    fig = plt.figure(figure_number, facecolor='White',figsize=(8,3))
    axa = plt.subplot2grid((1,3),(0,0))
    axa.set_title('a - M1a', fontsize=10)
    axa.set_ylabel('PAI',fontsize=axis_size)
    sns.violinplot(x='resolution',y='PAI_MH',data=df,linewidth=0.5,inner='box',
                    scale='width',hue="plot",palette = colour,dodge=False,saturation=0.7,cut=0)
    axa.set_xlabel('grid resolution')
    axa.set_ylabel('PAI',fontsize=axis_size)

    axb = plt.subplot2grid((1,3),(0,1),sharex=axa,sharey=axa)
    axb.set_title('b - M1b', fontsize=10)
    axb.set_ylabel('PAI',fontsize=axis_size)
    sns.violinplot(x='resolution',y='PAI_wt',data=df,linewidth=0.5,inner='box',
                    scale='width',hue="plot",palette = colour,dodge=False,saturation=0.7,cut=0)
    axb.set_xlabel('grid resolution')
    axb.set_ylabel('')

    axc = plt.subplot2grid((1,3),(0,2),sharex=axa,sharey=axa)
    axc.set_title('c - M2', fontsize=10)
    sns.violinplot(x='resolution',y='PAI_rad',data=df,linewidth=0.5,inner='box',
                    scale='width',hue="plot",palette = colour,dodge=False,saturation=0.7,cut=0)
    axc.set_ylabel('')
    axc.set_xlabel('grid resolution')
    axa.set_ylim(top=13)

    axes =[axa,axb,axc]
    plt.setp(axb.get_yticklabels(),visible=False);plt.setp(axc.get_yticklabels(),visible=False)

    for ax in axes:
        ax.yaxis.set_ticks_position('both')
        if ax!=axc:
            ax.legend_.remove()
        else:
            axc.legend(loc =9, bbox_to_anchor=(-0.7, -0.3),ncol=3)
    plt.subplots_adjust(wspace = 0.2,hspace=0.3,bottom=0.3)
    plt.savefig(figure_name)

    #plt.show()

def plot_PAI_sensitivity_density(figure_number,figure_name,PAD_profiles_MH_Belian,PAD_profiles_MH_BNorth,PAD_profiles_MH_E,
                    PAD_profiles_rad_Belian,PAD_profiles_rad_BNorth,PAD_profiles_rad_E,
                    PAD_profiles_wt_Belian,PAD_profiles_wt_BNorth,PAD_profiles_wt_E):

    res_keys = ['20m']; N_res = len(res_keys)
    dens_keys = ['5','10','15','20','25']; N_dens = len(dens_keys)
    dens = []; res = []; plot = []
    PAI_MH = np.arange(0); PAI_rad = np.arange(0); PAI_wt = np.arange(0)
    plot_name = ['MLA01','SAF03','SAF02']
    for rr in range(0,N_res):
        for dd in range(0,N_dens):
            for pp in range(0,len(plot_name)):
                PAD_profiles_wt_Belian[res_keys[rr]][dens_keys[dd]][~np.isfinite(PAD_profiles_wt_Belian[res_keys[rr]][dens_keys[dd]])]=np.nan
                PAD_profiles_wt_BNorth[res_keys[rr]][dens_keys[dd]][~np.isfinite(PAD_profiles_wt_BNorth[res_keys[rr]][dens_keys[dd]])]=np.nan
                PAD_profiles_wt_E[res_keys[rr]][dens_keys[dd]][~np.isfinite(PAD_profiles_wt_E[res_keys[rr]][dens_keys[dd]])]=np.nan
                if pp==0:
                    PAI_MH = np.append(PAI_MH,np.sum(np.nanmean(PAD_profiles_MH_Belian[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    PAI_wt = np.append(PAI_wt,np.sum(np.nanmean(PAD_profiles_wt_Belian[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    PAI_rad = np.append(PAI_rad,np.sum(np.nanmean(PAD_profiles_rad_Belian[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    N_obs = np.sum(np.nanmean(PAD_profiles_rad_Belian[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1).size
                    for ii in range(0,N_obs):
                        dens.append(dens_keys[dd])
                        res.append(res_keys[rr])
                        plot.append(plot_name[pp])
                if pp==1:
                    PAI_wt = np.append(PAI_wt,np.sum(np.nanmean(PAD_profiles_wt_E[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    PAI_MH = np.append(PAI_MH,np.sum(np.nanmean(PAD_profiles_MH_E[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    PAI_rad = np.append(PAI_rad,np.sum(np.nanmean(PAD_profiles_rad_E[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    N_obs = np.sum(np.nanmean(PAD_profiles_rad_E[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1).size
                    for ii in range(0,N_obs):
                        dens.append(dens_keys[dd])
                        res.append(res_keys[rr])
                        plot.append(plot_name[pp])
                if pp==2:
                    PAI_wt = np.append(PAI_wt,np.sum(np.nanmean(PAD_profiles_wt_BNorth[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    PAI_MH = np.append(PAI_MH,np.sum(np.nanmean(PAD_profiles_MH_BNorth[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    PAI_rad = np.append(PAI_rad,np.sum(np.nanmean(PAD_profiles_rad_BNorth[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                    N_obs = np.sum(np.nanmean(PAD_profiles_rad_BNorth[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1).size
                    for ii in range(0,N_obs):
                        dens.append(dens_keys[dd])
                        res.append(res_keys[rr])
                        plot.append(plot_name[pp])

    df = pd.DataFrame({'plot' : plot,'density' : dens,'resolution' : res,'PAI_MH':PAI_MH,'PAI_wt':PAI_wt,'PAI_rad':PAI_rad})
    df['density'][df['density']=='5']='05'
    pp=0
    fig = plt.figure(figure_number, facecolor='White',figsize=(8,3))
    axa = plt.subplot2grid((1,3),(0,0))
    axa.set_title('a - M1a', fontsize=10)
    axa.set_ylabel('PAI',fontsize=axis_size)
    sns.violinplot(x='density',y='PAI_MH',data=df,linewidth=0.5,inner='box',
                    scale='width',hue="plot",palette = colour,dodge=False,saturation=0.7,cut=0)
    axa.set_xlabel('pulse density / pts m$^{-2}$')
    axa.set_ylabel('PAI',fontsize=axis_size)

    axb = plt.subplot2grid((1,3),(0,1),sharex=axa,sharey=axa)
    axb.set_title('a - M1b', fontsize=10)
    axb.set_ylabel('PAI',fontsize=axis_size)
    sns.violinplot(x='density',y='PAI_wt',data=df,linewidth=0.5,inner='box',
                    scale='width',hue="plot",palette = colour,dodge=False,saturation=0.7,cut=0)
    axb.set_ylabel('')

    axb.set_xlabel('pulse density / pts m$^{-2}$')

    axc = plt.subplot2grid((1,3),(0,2),sharex=axa,sharey=axa)
    axc.set_title('b - M2', fontsize=10)
    sns.violinplot(x='density',y='PAI_rad',data=df,linewidth=0.5,inner='box',
                    scale='width',hue="plot",palette = colour,dodge=False,saturation=0.7,cut=0)
    axc.set_ylabel('')
    axc.set_xlabel('pulse density / pts m$^{-2}$')
    axa.set_ylim(top=13)

    axes =[axa,axb,axc]
    plt.setp(axb.get_yticklabels(),visible=False);plt.setp(axc.get_yticklabels(),visible=False)

    for ax in axes:
        ax.yaxis.set_ticks_position('both')
        if ax!=axc:
            ax.legend_.remove()
        else:
            axc.legend(loc =9, bbox_to_anchor=(-0.7, -0.3),ncol=3)
    plt.subplots_adjust(wspace = 0.2,hspace=0.3,bottom=0.3)
    plt.savefig(figure_name)

    #plt.show()
