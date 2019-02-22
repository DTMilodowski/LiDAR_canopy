import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rcParams
import matplotlib.ticker as mticker
import matplotlib.ticker as plticker
from matplotlib import colorbar
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

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
scale = np.arange(0.,6.)
scale /=5.
cmap_colour = cmap(scale)



"""
Plot sensitivity of 1ha profiles to spatial resolution
"""
def plot_profile_sensitivity_resolution(figure_number,figure_name,heights,PAD_profiles_MH_Belian,
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
    plt.show()
    return 0


"""
Plot sensitivity of 1ha profiles to point density
"""
def plot_profile_sensitivity_point_density(figure_number,figure_name,heights,PAD_profiles_MH_Belian,
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
    plt.show()
    return 0


"""
Plot fraction of unsampled layers in profile
"""
def plot_penetration_limits(figure_number,figure_name,heights,penetration_lim_Belian,
                            penetration_lim_E,penetration_lim_BNorth):
    ## Belian
    fig =plt.figure(figure_number, facecolor='White',figsize=[6,3])
    ax1a = plt.subplot2grid((1,3),(0,0))
    ax1a.set_ylabel('height / m',fontsize=axis_size)
    ax1a.annotate('a - MLA01\n(old growth)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1a.set_xlabel('occluded voxels / %',fontsize=axis_size)

    ax1b = plt.subplot2grid((1,3),(0,1),sharex=ax1a,sharey=ax1a)
    ax1b.annotate('b - SAF03\n(moderately logged)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1b.set_xlabel('occluded voxels / %',fontsize=axis_size)

    ax1c = plt.subplot2grid((1,3),(0,2),sharex=ax1a,sharey=ax1a)
    ax1c.annotate('c - SAF02\n(heavily logged)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
    ax1c.set_xlabel('occluded voxels / %',fontsize=axis_size)

    # loop through the subplots
    sp = [ax1a,ax1b,ax1c]
    pkeys = ['2m','5m','10m','20m','50m']
    for pp in range(0,len(pkeys)):

        lim_OG = np.nanmean(penetration_lim_Belian[pkeys[pp]]['40'],axis=1)*100.
        lim_ML = np.nanmean(penetration_lim_E[pkeys[pp]]['40'],axis=1)*100.
        lim_HL = np.nanmean(penetration_lim_BNorth[pkeys[pp]]['40'],axis=1)*100.

        limperc_OG = np.nanpercentile(lim_OG,[2.5,50,97.5],axis=0)
        limperc_ML = np.nanpercentile(lim_ML,[2.5,50,97.5],axis=0)
        limperc_HL = np.nanpercentile(lim_HL,[2.5,50,97.5],axis=0)

        ax1a.fill_betweenx(heights[2:],limperc_OG[0][2:],limperc_OG[2][2:],color=cmap_colour[pp],alpha=0.3)
        ax1a.plot(limperc_OG[1][2:],heights[2:],'-',c=cmap_colour[pp],linewidth=1.5)
        ax1b.fill_betweenx(heights[2:],limperc_ML[0][2:],limperc_ML[2][2:],color=cmap_colour[pp],alpha=0.3)
        ax1b.plot(limperc_ML[1][2:],heights[2:],'-',c=cmap_colour[pp],linewidth=1.5)
        ax1c.fill_betweenx(heights[2:],limperc_HL[0][2:],limperc_HL[2][2:],color=cmap_colour[pp],alpha=0.3)
        ax1c.plot(limperc_HL[1][2:],heights[2:],'-',c=cmap_colour[pp],linewidth=1.5,label=pkeys[pp])

    ax1a.set_ylim(0,89)
    ax1a.set_xlim(0,100)

    yticklabels = ax1b.get_yticklabels() + ax1c.get_yticklabels()
    plt.setp(yticklabels,visible=False)

    for ax in sp:
        ax.locator_params(axis='x',nbins=5)
        ax.yaxis.set_ticks_position('both')

    plt.subplots_adjust(hspace=0.4, wspace = 0.4,bottom=0.2)

    ax1c.legend(loc=5)
    plt.savefig(figure_name)
    plt.show()
    return 0
