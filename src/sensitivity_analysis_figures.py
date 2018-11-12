import numpy as np
import sys
from matplotlib import pyplot as plt
from matplotlib import rcParams
from scipy import stats

# Set up some basiic parameters for the plots
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = 8
rcParams['legend.numpoints'] = 1
axis_size = rcParams['font.size']+2
colour = ['#46E900','#1A2BCE','#E0007F']

# Load in the sensitivity analysis results
max_height=80
heights = np.arange(0.,max_height)+1

PAD_profiles_MH_Belian = np.load('MH_sensitivity_Belian.npy')[()]
PAD_profiles_MH_BNorth = np.load('MH_sensitivity_BNorth.npy')[()]
PAD_profiles_MH_E = np.load('MH_sensitivity_E.npy')[()]

PAD_profiles_rad1_Belian = np.load('rad1_sensitivity_Belian.npy')[()]
PAD_profiles_rad1_BNorth = np.load('rad1_sensitivity_BNorth.npy')[()]
PAD_profiles_rad1_E = np.load('rad1_sensitivity_E.npy')[()]

PAD_profiles_rad2_Belian = np.load('rad2_sensitivity_Belian.npy')[()]
PAD_profiles_rad2_BNorth = np.load('rad2_sensitivity_BNorth.npy')[()]
PAD_profiles_rad2_E = np.load('rad2_sensitivity_E.npy')[()]

penetration_lim_Belian = np.load('penetration_limit_Belian.npy')[()]
penetration_lim_BNorth = np.load('penetration_limit_BNorth.npy')[()]
penetration_lim_E = np.load('penetration_limit_E.npy')[()]

# Figure 1 - sensitivity of 1 ha profiles to point density
## Belian
plt.figure(1, facecolor='White',figsize=[8,9])
ax1a = plt.subplot2grid((3,6),(0,0))
ax1a.set_title('MLA01 - old growth', fontsize=10)
ax1a.set_ylabel('height / m',fontsize=axis_size)
ax1a.annotate('a - 2 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax1b = plt.subplot2grid((3,6),(0,1),sharex=ax1a,sharey=ax1a)
ax1b.annotate('b - 5 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax1c = plt.subplot2grid((3,6),(0,2),sharex=ax1a,sharey=ax1a)
ax1c.annotate('c - 10 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax1d = plt.subplot2grid((3,6),(0,3),sharex=ax1a,sharey=ax1a)
ax1d.annotate('d - 20 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax1e = plt.subplot2grid((3,6),(0,4),sharex=ax1a,sharey=ax1a)
ax1e.annotate('e - 50 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax1f = plt.subplot2grid((3,6),(0,5),sharex=ax1a,sharey=ax1a)
ax1f.annotate('f - 100 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

## E
ax1g = plt.subplot2grid((3,6),(1,0),sharex=ax1a,sharey=ax1a)
ax1g.set_title('SAF03 - moderately logged', fontsize=10)
ax1g.set_ylabel('height / m',fontsize=axis_size)
ax1g.annotate('g - 2 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax1h = plt.subplot2grid((3,6),(1,1),sharex=ax1a,sharey=ax1a)
ax1h.annotate('h - 5 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax1i = plt.subplot2grid((3,6),(1,2),sharex=ax1a,sharey=ax1a)
ax1i.annotate('i - 10 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax1j = plt.subplot2grid((3,6),(1,3),sharex=ax1a,sharey=ax1a)
ax1j.annotate('j - 20 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax1k = plt.subplot2grid((3,6),(1,4),sharex=ax1a,sharey=ax1a)
ax1k.annotate('k - 50 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax1l = plt.subplot2grid((3,6),(1,5),sharex=ax1a,sharey=ax1a)
ax1l.annotate('l - 100 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

## B North
ax1m = plt.subplot2grid((3,6),(2,0),sharex=ax1a,sharey=ax1a)
ax1m.set_title('SAF02 - heavily logged', fontsize=10)
ax1m.set_xlabel('PAD / m$^2$m$^{-3}$',fontsize=axis_size)
ax1m.set_ylabel('height / m',fontsize=axis_size)
ax1m.annotate('m - 2 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax1n = plt.subplot2grid((3,6),(2,1),sharex=ax1a,sharey=ax1a)
ax1n.set_xlabel('PAD / m$^2$m$^{-3}$',fontsize=axis_size)
ax1n.annotate('n - 5 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax1o = plt.subplot2grid((3,6),(2,2),sharex=ax1a,sharey=ax1a)
ax1o.set_xlabel('PAD / m$^2$m$^{-3}$',fontsize=axis_size)
ax1o.annotate('o - 10 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax1p = plt.subplot2grid((3,6),(2,3),sharex=ax1a,sharey=ax1a)
ax1p.set_xlabel('PAD / m$^2$m$^{-3}$',fontsize=axis_size)
ax1p.annotate('p - 20 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax1q = plt.subplot2grid((3,6),(2,4),sharex=ax1a,sharey=ax1a)
ax1q.set_xlabel('PAD / m$^2$m$^{-3}$',fontsize=axis_size)
ax1q.annotate('q - 50 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax1r = plt.subplot2grid((3,6),(2,5),sharex=ax1a,sharey=ax1a)
ax1r.set_xlabel('PAD / m$^2$m$^{-3}$',fontsize=axis_size)
ax1r.annotate('r - 100 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

# loop through the subplots
sp = [ax1a,ax1b,ax1c,ax1d,ax1e,ax1f,ax1g,ax1h,ax1i,ax1j,ax1k,ax1l,ax1m,ax1n,ax1o,ax1p,ax1q,ax1r]
pkeys = ['2m','5m','10m','20m','50m','100m']
for pp in range(0,len(sp)):

    if pp<6:
        rad1_sp = np.nanmean(PAD_profiles_rad1_Belian[pkeys[pp]]['40'],axis=1)
        rad2_sp = np.nanmean(PAD_profiles_rad2_Belian[pkeys[pp]]['40'],axis=1)
        MH_sp = np.nanmean(PAD_profiles_MH_Belian[pkeys[pp]]['40'],axis=1)
    elif pp<12:
        rad1_sp = np.nanmean(PAD_profiles_rad1_E[pkeys[pp-6]]['40'],axis=1)
        rad2_sp = np.nanmean(PAD_profiles_rad2_E[pkeys[pp-6]]['40'],axis=1)
        MH_sp = np.nanmean(PAD_profiles_MH_E[pkeys[pp-6]]['40'],axis=1)
    else:
        rad1_sp = np.nanmean(PAD_profiles_rad1_BNorth[pkeys[pp-12]]['40'],axis=1)
        rad2_sp = np.nanmean(PAD_profiles_rad2_BNorth[pkeys[pp-12]]['40'],axis=1)
        MH_sp = np.nanmean(PAD_profiles_MH_BNorth[pkeys[pp-12]]['40'],axis=1)


    rad1 = np.nanmedian(rad1_sp,axis=0)
    rad1perc = np.nanpercentile(rad1_sp,[2.5,97.5],axis=0)
    sp[pp].fill_betweenx(heights[2:],rad1perc[0][2:],rad1perc[1][2:],color=colour[1],alpha=0.3)

    rad2 = np.nanmedian(rad2_sp,axis=0)
    rad2perc = np.nanpercentile(rad2_sp,[2.5,97.5],axis=0)
    sp[pp].fill_betweenx(heights[2:],rad2perc[0][2:],rad2perc[1][2:],color=colour[2],alpha=0.3)

    MH = np.nanmedian(MH_sp,axis=0)
    MHperc = np.nanpercentile(MH_sp,[2.5,97.5],axis=0)
    sp[pp].fill_betweenx(heights[2:],MHperc[0][2:],MHperc[1][2:],color=colour[0],alpha=0.3)

    if pp==14:
        sp[pp].plot(MH[2:],heights[2:],'-',c=colour[0],linewidth=1.5, label = 'MacArthur-Horn')
        sp[pp].plot(rad1[2:],heights[2:],'-',c=colour[1],linewidth=1.5, label = 'rad trans (Detto)')
        sp[pp].plot(rad2[2:],heights[2:],'-',c=colour[2],linewidth=1.5, label = 'rad trans (modified)')
        sp[pp].legend(loc=9, bbox_to_anchor=(0.5, -0.3), ncol=3)
    else:
        sp[pp].plot(MH[2:],heights[2:],'-',c=colour[0],linewidth=1.5)
        sp[pp].plot(rad1[2:],heights[2:],'-',c=colour[1],linewidth=1.5)
        sp[pp].plot(rad2[2:],heights[2:],'-',c=colour[2],linewidth=1.5)

ax1a.set_ylim(0,80)
ax1a.set_xlim(0,0.59)

ax1a.locator_params(axis='x',nbins=5)
ax1b.locator_params(axis='x',nbins=5)
ax1c.locator_params(axis='x',nbins=5)
ax1d.locator_params(axis='x',nbins=5)
ax1e.locator_params(axis='x',nbins=5)
ax1f.locator_params(axis='x',nbins=5)
ax1g.locator_params(axis='x',nbins=5)
ax1h.locator_params(axis='x',nbins=5)
ax1i.locator_params(axis='x',nbins=5)
ax1j.locator_params(axis='x',nbins=5)
ax1k.locator_params(axis='x',nbins=5)
ax1l.locator_params(axis='x',nbins=5)
ax1m.locator_params(axis='x',nbins=5)
ax1n.locator_params(axis='x',nbins=5)
ax1o.locator_params(axis='x',nbins=5)
ax1p.locator_params(axis='x',nbins=5)
ax1q.locator_params(axis='x',nbins=5)
ax1r.locator_params(axis='x',nbins=5)

yticklabels = ax1b.get_yticklabels() + ax1c.get_yticklabels() + ax1d.get_yticklabels() + ax1e.get_yticklabels() + ax1f.get_yticklabels() + \
              ax1h.get_yticklabels() + ax1i.get_yticklabels() + ax1j.get_yticklabels() + ax1k.get_yticklabels() + ax1l.get_yticklabels() + \
              ax1n.get_yticklabels() + ax1o.get_yticklabels() + ax1p.get_yticklabels() + ax1q.get_yticklabels() + ax1r.get_yticklabels()
xticklabels = ax1a.get_xticklabels() + ax1b.get_xticklabels() + ax1c.get_xticklabels() + ax1d.get_xticklabels() + ax1e.get_xticklabels() + ax1f.get_xticklabels() + \
              ax1g.get_xticklabels() + ax1h.get_xticklabels() + ax1i.get_xticklabels() + ax1j.get_xticklabels() + ax1k.get_xticklabels() + ax1l.get_xticklabels()
plt.setp(yticklabels,visible=False)
#plt.setp(xticklabels,visible=False)
plt.subplots_adjust(hspace=0.2, wspace = 0.2, bottom = 0.2)
plt.savefig('PAD_resolution_sensitivity_1ha_average.png')
plt.show()





#==============================================================================================
# 2) Equivalent plots but this time using different point densities at target resolution
pkeys=['5','10','15','20','25','30','40']
plt.figure(2, facecolor='White',figsize=[7,9])
ax2a = plt.subplot2grid((3,5),(0,0))
ax2a.set_title('MLA01 - old growth', fontsize=10)
ax2a.set_ylabel('height / m',fontsize=axis_size)
ax2a.annotate('a - 5 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax2b = plt.subplot2grid((3,5),(0,1),sharex=ax2a,sharey=ax2a)
ax2b.annotate('b - 10 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax2c = plt.subplot2grid((3,5),(0,2),sharex=ax2a,sharey=ax2a)
ax2c.annotate('c - 20 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax2d = plt.subplot2grid((3,5),(0,3),sharex=ax2a,sharey=ax2a)
ax2d.annotate('d - 30 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax2e = plt.subplot2grid((3,5),(0,4),sharex=ax2a,sharey=ax2a)
ax2e.annotate('e - 40 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax2f = plt.subplot2grid((3,5),(1,0),sharex=ax2a,sharey=ax2a)
ax2f.set_title('SAF03 - moderately logged', fontsize=10)
ax2f.set_ylabel('height / m',fontsize=axis_size)
ax2f.annotate('f - 5 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax2g = plt.subplot2grid((3,5),(1,1),sharex=ax2a,sharey=ax2a)
ax2g.annotate('g - 10 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax2h = plt.subplot2grid((3,5),(1,2),sharex=ax2a,sharey=ax2a)
ax2h.annotate('h - 20 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax2i = plt.subplot2grid((3,5),(1,3),sharex=ax2a,sharey=ax2a)
ax2i.annotate('i - 30 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax2j = plt.subplot2grid((3,5),(1,4),sharex=ax2a,sharey=ax2a)
ax2j.annotate('j - 40 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax2k = plt.subplot2grid((3,5),(2,0),sharex=ax2a,sharey=ax2a)
ax2k.set_title('SAF02 - heavily logged', fontsize=10)
ax2k.set_xlabel('PAD / m$^2$m$^{-3}$',fontsize=axis_size)
ax2k.set_ylabel('height / m',fontsize=axis_size)
ax2k.annotate('k - 5 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax2l = plt.subplot2grid((3,5),(2,1),sharex=ax2a,sharey=ax2a)
ax2l.set_xlabel('PAD / m$^2$m$^{-3}$',fontsize=axis_size)
ax2l.annotate('l - 10 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax2m = plt.subplot2grid((3,5),(2,2),sharex=ax2a,sharey=ax2a)
ax2m.set_xlabel('PAD / m$^2$m$^{-3}$',fontsize=axis_size)
ax2m.annotate('m - 20 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax2n = plt.subplot2grid((3,5),(2,3),sharex=ax2a,sharey=ax2a)
ax2n.set_xlabel('PAD / m$^2$m$^{-3}$',fontsize=axis_size)
ax2n.annotate('n - 30 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax2o = plt.subplot2grid((3,5),(2,4),sharex=ax2a,sharey=ax2a)
ax2o.set_xlabel('PAD / m$^2$m$^{-3}$',fontsize=axis_size)
ax2o.annotate('o - 40 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

# loop through the subplots
sp = [ax2a,ax2b,ax2c,ax2d,ax2e,ax2f,ax2g,ax2h,ax2i,ax2j,ax2k,ax2l,ax2m,ax2n,ax2o]
pkeys = ['5','10','20','30','40']
for pp in range(0,len(sp)):
    if pp<5:
        rad1_sp = np.nanmean(PAD_profiles_rad1_Belian['20m'][pkeys[pp]],axis=1)
        rad2_sp = np.nanmean(PAD_profiles_rad2_Belian['20m'][pkeys[pp]],axis=1)
        MH_sp = np.nanmean(PAD_profiles_MH_Belian['20m'][pkeys[pp]],axis=1)
    elif pp<10:
        rad1_sp = np.nanmean(PAD_profiles_rad1_E['20m'][pkeys[pp-5]],axis=1)
        rad2_sp = np.nanmean(PAD_profiles_rad2_E['20m'][pkeys[pp-5]],axis=1)
        MH_sp = np.nanmean(PAD_profiles_MH_E['20m'][pkeys[pp-5]],axis=1)
    else:
        rad1_sp = np.nanmean(PAD_profiles_rad1_BNorth['20m'][pkeys[pp-10]],axis=1)
        rad2_sp = np.nanmean(PAD_profiles_rad2_BNorth['20m'][pkeys[pp-10]],axis=1)
        MH_sp = np.nanmean(PAD_profiles_MH_BNorth['20m'][pkeys[pp-10]],axis=1)

    rad2 = np.nanmedian(rad2_sp,axis=0)
    rad2perc = np.nanpercentile(rad2_sp,[2.5,97.5],axis=0)
    sp[pp].fill_betweenx(heights[2:],rad2perc[0][2:],rad2perc[1][2:],color=colour[2],alpha=0.3)

    rad1 = np.nanmedian(rad1_sp,axis=0)
    rad1perc = np.percentile(rad1_sp,[2.5,97.5],axis=0)
    sp[pp].fill_betweenx(heights[2:],rad1perc[0][2:],rad1perc[1][2:],color=colour[1],alpha=0.3)

    MH = np.median(MH_sp,axis=0)
    MHperc = np.percentile(MH_sp,[2.5,97.5],axis=0)
    sp[pp].fill_betweenx(heights[2:],MHperc[0][2:],MHperc[1][2:],color=colour[0],alpha=0.3)

    if pp==12:
        sp[pp].plot(MH[2:],heights[2:],'-',c=colour[0],linewidth=1.5, label = 'MacArthur-Horn')
        sp[pp].plot(rad1[2:],heights[2:],'-',c=colour[1],linewidth=1.5, label = 'rad trans (Detto)')
        sp[pp].plot(rad2[2:],heights[2:],'-',c=colour[2],linewidth=1.5, label = 'rad trans (modified)')
        sp[pp].legend(loc=9, bbox_to_anchor=(0.5, -0.3), ncol=3)#(loc='upper right')
    else:
        sp[pp].plot(MH[2:],heights[2:],'-',c=colour[0],linewidth=1.5)
        sp[pp].plot(rad1[2:],heights[2:],'-',c=colour[1],linewidth=1.5)
        sp[pp].plot(rad2[2:],heights[2:],'-',c=colour[2],linewidth=1.5)

ax2a.set_ylim(0,80)
ax2a.set_xlim(0,0.59)
ax2a.locator_params(axis='x',nbins=5)
ax2b.locator_params(axis='x',nbins=5)
ax2c.locator_params(axis='x',nbins=5)
ax2d.locator_params(axis='x',nbins=5)
ax2e.locator_params(axis='x',nbins=5)
ax2f.locator_params(axis='x',nbins=5)
ax2g.locator_params(axis='x',nbins=5)
ax2h.locator_params(axis='x',nbins=5)
ax2i.locator_params(axis='x',nbins=5)
ax2j.locator_params(axis='x',nbins=5)
ax2k.locator_params(axis='x',nbins=5)
ax2l.locator_params(axis='x',nbins=5)
ax2m.locator_params(axis='x',nbins=5)
ax2n.locator_params(axis='x',nbins=5)
ax2o.locator_params(axis='x',nbins=5)

yticklabels = ax2b.get_yticklabels() + ax2c.get_yticklabels() + ax2d.get_yticklabels() + ax2e.get_yticklabels() + \
              ax2g.get_yticklabels() + ax2h.get_yticklabels() + ax2i.get_yticklabels() + ax2j.get_yticklabels() + \
              ax2l.get_yticklabels() + ax2m.get_yticklabels() + ax2n.get_yticklabels() + ax2o.get_yticklabels()
xticklabels = ax1a.get_xticklabels() + ax1b.get_xticklabels() + ax1c.get_xticklabels() + ax1d.get_xticklabels() + ax1e.get_xticklabels() + \
              ax1f.get_xticklabels() + ax1g.get_xticklabels() + ax1h.get_xticklabels() + ax1i.get_xticklabels() + ax1j.get_xticklabels()
plt.setp(yticklabels,visible=False)
plt.setp(xticklabels,visible=False)
plt.subplots_adjust(hspace=0.2, wspace = 0.2, bottom = 0.2)
plt.savefig('PAD_point_density_sensitivity_1ha_averages.png')
plt.show()

#===============================================================================
# Sensitivity of individual profiles to resolution - only underaken at one plot
# as these trends are common across the three.
## Belian
plt.figure(3, facecolor='White',figsize=[8,4])
ax3a = plt.subplot2grid((1,6),(0,0))
ax3a.set_ylabel('height / m',fontsize=axis_size)
ax3a.annotate('a - 2 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax3b = plt.subplot2grid((1,6),(0,1),sharex=ax3a,sharey=ax3a)
ax3b.annotate('b - 5 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax3c = plt.subplot2grid((1,6),(0,2),sharex=ax3a,sharey=ax3a)
ax3c.annotate('c - 10 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax3d = plt.subplot2grid((1,6),(0,3),sharex=ax3a,sharey=ax3a)
ax3d.annotate('d - 20 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax3e = plt.subplot2grid((1,6),(0,4),sharex=ax3a,sharey=ax3a)
ax3e.annotate('e - 50 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax3f = plt.subplot2grid((1,6),(0,5),sharex=ax3a,sharey=ax3a)
ax3f.annotate('f - 100 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

# loop through the subplots
sp = [ax3a,ax3b,ax3c,ax3d,ax3e,ax3f]
pkeys = ['2m','5m','10m','20m','50m','100m']
for pp in range(0,len(sp)):
    MHl = np.nanpercentile(PAD_profiles_MH_Belian[pkeys[pp]]['40'],2.5,axis=0)
    MH25 = np.nanpercentile(PAD_profiles_MH_Belian[pkeys[pp]]['40'],25,axis=0)
    MH75 = np.nanpercentile(PAD_profiles_MH_Belian[pkeys[pp]]['40'],75,axis=0)
    MHu = np.nanpercentile(PAD_profiles_MH_Belian[pkeys[pp]]['40'],97.5,axis=0)
    MH_95CI_i = ((MHu-MHl)/np.nanmean(PAD_profiles_MH_Belian[pkeys[pp]]['40'],axis=0))*100.
    MH_50CI_i = ((MH75-MH25)/np.nanmean(PAD_profiles_MH_Belian[pkeys[pp]]['40'],axis=0))*100.
    MH_95CI = np.nansum(MH_95CI_i,axis=0)/np.sum(np.isfinite(MH_95CI_i),axis=0)
    MH_50CI = np.nansum(MH_50CI_i,axis=0)/np.sum(np.isfinite(MH_50CI_i),axis=0)
    sp[pp].plot(MH_95CI[2:],heights[2:],':',c=colour[0],linewidth=1)
    sp[pp].plot(MH_50CI[2:],heights[2:],'-',c=colour[0],linewidth=1)
    if pp==2:
        sp[pp].plot(MH_95CI[2:],heights[2:],':',c=colour[0],linewidth=1,label='MacArthur-Horn 95% CI')
        sp[pp].plot(MH_50CI[2:],heights[2:],'-',c=colour[0],linewidth=1,label='MacArthur-Horn 50% CI')
    else:
        sp[pp].plot(MH_95CI[2:],heights[2:],':',c=colour[0],linewidth=1)
        sp[pp].plot(MH_50CI[2:],heights[2:],'-',c=colour[0],linewidth=1)

    rad1l = np.nanpercentile(PAD_profiles_rad1_Belian[pkeys[pp]]['40'],2.5,axis=0)
    rad125 = np.nanpercentile(PAD_profiles_rad1_Belian[pkeys[pp]]['40'],25,axis=0)
    rad175 = np.nanpercentile(PAD_profiles_rad1_Belian[pkeys[pp]]['40'],75,axis=0)
    rad1u = np.nanpercentile(PAD_profiles_rad1_Belian[pkeys[pp]]['40'],97.5,axis=0)
    rad1_95CI_i = ((rad1u-rad1l)/np.nanmean(PAD_profiles_rad1_Belian[pkeys[pp]]['40'],axis=0))*100.
    rad1_50CI_i = ((rad175-rad125)/np.nanmean(PAD_profiles_rad1_Belian[pkeys[pp]]['40'],axis=0))*100.
    rad1_95CI = np.nansum(rad1_95CI_i,axis=0)/np.sum(np.isfinite(rad1_95CI_i),axis=0)
    rad1_50CI = np.nansum(rad1_50CI_i,axis=0)/np.sum(np.isfinite(rad1_50CI_i),axis=0)
    if pp==2:
        sp[pp].plot(rad1_95CI[2:],heights[2:],':',c=colour[1],linewidth=1,label='rad. trans. (Detto) 95% CI')
        sp[pp].plot(rad1_50CI[2:],heights[2:],'-',c=colour[1],linewidth=1,label='rad. trans. (Detto) 50% CI')
    else:
        sp[pp].plot(rad1_95CI[2:],heights[2:],':',c=colour[1],linewidth=1)
        sp[pp].plot(rad1_50CI[2:],heights[2:],'-',c=colour[1],linewidth=1)

    rad2l = np.nanpercentile(PAD_profiles_rad2_Belian[pkeys[pp]]['40'],2.5,axis=0)
    rad225 = np.nanpercentile(PAD_profiles_rad2_Belian[pkeys[pp]]['40'],25,axis=0)
    rad275 = np.nanpercentile(PAD_profiles_rad2_Belian[pkeys[pp]]['40'],75,axis=0)
    rad2u = np.nanpercentile(PAD_profiles_rad2_Belian[pkeys[pp]]['40'],97.5,axis=0)
    rad2_95CI_i = ((rad2u-rad2l)/np.nanmean(PAD_profiles_rad2_Belian[pkeys[pp]]['40'],axis=0))*100.
    rad2_50CI_i = ((rad275-rad225)/np.nanmean(PAD_profiles_rad2_Belian[pkeys[pp]]['40'],axis=0))*100.
    rad2_95CI = np.nansum(rad2_95CI_i,axis=0)/np.sum(np.isfinite(rad2_95CI_i),axis=0)
    rad2_50CI = np.nansum(rad2_50CI_i,axis=0)/np.sum(np.isfinite(rad2_50CI_i),axis=0)
    if pp==2:
        sp[pp].plot(rad2_95CI[2:],heights[2:],':',c=colour[2],linewidth=1,label='rad. trans. (mod) 95% CI')
        sp[pp].plot(rad2_50CI[2:],heights[2:],'-',c=colour[2],linewidth=1,label='rad. trans. (mod) 50% CI')
        lgd = sp[pp].legend(loc=9, bbox_to_anchor=(0.5, -0.2), ncol=3)#(loc='upper right')
    else:
        sp[pp].plot(rad2_95CI[2:],heights[2:],':',c=colour[2],linewidth=1)
        sp[pp].plot(rad2_50CI[2:],heights[2:],'-',c=colour[2],linewidth=1)

ax3a.set_ylim(0,89)
ax3a.set_xlim(0,520)
ax3a.locator_params(axis='x',nbins=5)
ax3b.locator_params(axis='x',nbins=5)
ax3c.locator_params(axis='x',nbins=5)
ax3d.locator_params(axis='x',nbins=5)
ax3e.locator_params(axis='x',nbins=5)
ax3f.locator_params(axis='x',nbins=5)

yticklabels = ax3b.get_yticklabels() + ax3c.get_yticklabels() + ax3d.get_yticklabels() + ax3e.get_yticklabels() + ax3f.get_yticklabels()
plt.setp(yticklabels,visible=False)
plt.savefig('PAD_resolution_sensitivity_individual_CIs.png')
plt.show()


#===============================================================================
# As above, but now looking at point density
## Belian
plt.figure(4, facecolor='White',figsize=[7,3])
ax4a = plt.subplot2grid((1,5),(0,0))
ax4a.set_ylabel('height / m',fontsize=axis_size)
ax4a.annotate('a - 5 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax4b = plt.subplot2grid((1,5),(0,1),sharex=ax4a,sharey=ax4a)
ax4b.annotate('b - 10 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax4c = plt.subplot2grid((1,5),(0,2),sharex=ax4a,sharey=ax4a)
ax4c.annotate('c - 20 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax4d = plt.subplot2grid((1,5),(0,3),sharex=ax4a,sharey=ax4a)
ax4d.annotate('d - 30 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax4e = plt.subplot2grid((1,5),(0,4),sharex=ax4a,sharey=ax4a)
ax4e.annotate('e - 40 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

# loop through the subplots
sp = [ax4a,ax4b,ax4c,ax4d,ax4e]
pkeys = ['5','10','20','30','40']
for pp in range(0,len(sp)):
    MHl = np.nanpercentile(PAD_profiles_MH_Belian['20m'][pkeys[pp]],2.5,axis=0)
    MH25 = np.nanpercentile(PAD_profiles_MH_Belian['20m'][pkeys[pp]],25,axis=0)
    MH75 = np.nanpercentile(PAD_profiles_MH_Belian['20m'][pkeys[pp]],75,axis=0)
    MHu = np.nanpercentile(PAD_profiles_MH_Belian['20m'][pkeys[pp]],97.5,axis=0)
    MH_95CI_i = ((MHu-MHl)/np.nanmean(PAD_profiles_MH_Belian['20m'][pkeys[pp]],axis=0))*100.
    MH_50CI_i = ((MH75-MH25)/np.nanmean(PAD_profiles_MH_Belian['20m'][pkeys[pp]],axis=0))*100.
    MH_95CI = np.nansum(MH_95CI_i,axis=0)/np.sum(np.isfinite(MH_95CI_i),axis=0)
    MH_50CI = np.nansum(MH_50CI_i,axis=0)/np.sum(np.isfinite(MH_50CI_i),axis=0)
    sp[pp].plot(MH_95CI[2:],heights[2:],':',c=colour[0],linewidth=1)
    sp[pp].plot(MH_50CI[2:],heights[2:],'-',c=colour[0],linewidth=1)
    if pp==2:
        sp[pp].plot(MH_95CI[2:],heights[2:],':',c=colour[0],linewidth=1,label='MacArthur-Horn 95% CI')
        sp[pp].plot(MH_50CI[2:],heights[2:],'-',c=colour[0],linewidth=1,label='MacArthur-Horn 50% CI')
    else:
        sp[pp].plot(MH_95CI[2:],heights[2:],':',c=colour[0],linewidth=1)
        sp[pp].plot(MH_50CI[2:],heights[2:],'-',c=colour[0],linewidth=1)

    rad1l = np.nanpercentile(PAD_profiles_rad1_Belian['20m'][pkeys[pp]],2.5,axis=0)
    rad125 = np.nanpercentile(PAD_profiles_rad1_Belian['20m'][pkeys[pp]],25,axis=0)
    rad175 = np.nanpercentile(PAD_profiles_rad1_Belian['20m'][pkeys[pp]],75,axis=0)
    rad1u = np.nanpercentile(PAD_profiles_rad1_Belian['20m'][pkeys[pp]],97.5,axis=0)
    rad1_95CI_i = ((rad1u-rad1l)/np.nanmean(PAD_profiles_rad1_Belian['20m'][pkeys[pp]],axis=0))*100.
    rad1_50CI_i = ((rad175-rad125)/np.nanmean(PAD_profiles_rad1_Belian['20m'][pkeys[pp]],axis=0))*100.
    rad1_95CI = np.nansum(rad1_95CI_i,axis=0)/np.sum(np.isfinite(rad1_95CI_i),axis=0)
    rad1_50CI = np.nansum(rad1_50CI_i,axis=0)/np.sum(np.isfinite(rad1_50CI_i),axis=0)
    if pp==2:
        sp[pp].plot(rad1_95CI[2:],heights[2:],':',c=colour[1],linewidth=1,label='rad. trans. (Detto) 95% CI')
        sp[pp].plot(rad1_50CI[2:],heights[2:],'-',c=colour[1],linewidth=1,label='rad. trans. (Detto) 50% CI')
    else:
        sp[pp].plot(rad1_95CI[2:],heights[2:],':',c=colour[1],linewidth=1)
        sp[pp].plot(rad1_50CI[2:],heights[2:],'-',c=colour[1],linewidth=1)

    rad2l = np.nanpercentile(PAD_profiles_rad2_Belian['20m'][pkeys[pp]],2.5,axis=0)
    rad225 = np.nanpercentile(PAD_profiles_rad2_Belian['20m'][pkeys[pp]],25,axis=0)
    rad275 = np.nanpercentile(PAD_profiles_rad2_Belian['20m'][pkeys[pp]],75,axis=0)
    rad2u = np.nanpercentile(PAD_profiles_rad2_Belian['20m'][pkeys[pp]],97.5,axis=0)
    rad2_95CI_i = ((rad2u-rad2l)/np.nanmean(PAD_profiles_rad2_Belian['20m'][pkeys[pp]],axis=0))*100.
    rad2_50CI_i = ((rad275-rad225)/np.nanmean(PAD_profiles_rad2_Belian['20m'][pkeys[pp]],axis=0))*100.
    rad2_95CI = np.nansum(rad2_95CI_i,axis=0)/np.sum(np.isfinite(rad2_95CI_i),axis=0)
    rad2_50CI = np.nansum(rad2_50CI_i,axis=0)/np.sum(np.isfinite(rad2_50CI_i),axis=0)
    if pp==2:
        sp[pp].plot(rad2_95CI[2:],heights[2:],':',c=colour[2],linewidth=1,label='rad. trans. (mod) 95% CI')
        sp[pp].plot(rad2_50CI[2:],heights[2:],'-',c=colour[2],linewidth=1,label='rad. trans. (mod) 50% CI')
        lgd = sp[pp].legend(loc=9, bbox_to_anchor=(0.5, -0.2), ncol=3)#(loc='upper right')
    else:
        sp[pp].plot(rad2_95CI[2:],heights[2:],':',c=colour[2],linewidth=1)
        sp[pp].plot(rad2_50CI[2:],heights[2:],'-',c=colour[2],linewidth=1)

ax4a.set_ylim(0,89)
ax4a.set_xlim(0,520)
ax4a.locator_params(axis='x',nbins=5)
ax4b.locator_params(axis='x',nbins=5)
ax4c.locator_params(axis='x',nbins=5)
ax4d.locator_params(axis='x',nbins=5)
ax4e.locator_params(axis='x',nbins=5)

yticklabels = ax4b.get_yticklabels() + ax4c.get_yticklabels() + ax4d.get_yticklabels() + ax4e.get_yticklabels()
plt.setp(yticklabels,visible=False)
plt.savefig('PAD_point_density_sensitivity_individual_CIs.png')
plt.show()




#===============================================================================
# 5) Penetration limits
# A plot showing the penetration limits for the LiDAR surveys at different
# resolutions
plt.figure(5, facecolor='White',figsize=[8,9])
ax1a = plt.subplot2grid((3,6),(0,0))
ax1a.set_ylabel('height / m',fontsize=axis_size)
ax1a.set_title('MLA01', fontsize=10)
ax1a.annotate('a - 2 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax1b = plt.subplot2grid((3,6),(0,1),sharex=ax1a,sharey=ax1a)
ax1b.annotate('b - 5 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax1c = plt.subplot2grid((3,6),(0,2),sharex=ax1a,sharey=ax1a)
ax1c.annotate('c - 10 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax1d = plt.subplot2grid((3,6),(0,3),sharex=ax1a,sharey=ax1a)
ax1d.annotate('d - 20 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax1e = plt.subplot2grid((3,6),(0,4),sharex=ax1a,sharey=ax1a)
ax1e.annotate('e - 50 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax1f = plt.subplot2grid((3,6),(0,5),sharex=ax1a,sharey=ax1a)
ax1f.annotate('f - 100 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

## E
ax1g = plt.subplot2grid((3,6),(1,0),sharex=ax1a,sharey=ax1a)
ax1g.set_ylabel('height / m',fontsize=axis_size)
ax1g.set_title('SAF03', fontsize=10)
ax1g.annotate('g - 2 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax1h = plt.subplot2grid((3,6),(1,1),sharex=ax1a,sharey=ax1a)
ax1h.annotate('h - 5 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax1i = plt.subplot2grid((3,6),(1,2),sharex=ax1a,sharey=ax1a)
ax1i.annotate('i - 10 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax1j = plt.subplot2grid((3,6),(1,3),sharex=ax1a,sharey=ax1a)
ax1j.annotate('j - 20 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax1k = plt.subplot2grid((3,6),(1,4),sharex=ax1a,sharey=ax1a)
ax1k.annotate('k - 50 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax1l = plt.subplot2grid((3,6),(1,5),sharex=ax1a,sharey=ax1a)
ax1l.annotate('l - 100 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

## B North
ax1m = plt.subplot2grid((3,6),(2,0),sharex=ax1a,sharey=ax1a)
ax1m.set_title('SAF02', fontsize=10)
ax1m.set_xlabel('fraction occluded (%)',fontsize=axis_size)
ax1m.set_ylabel('height / m',fontsize=axis_size)
ax1m.annotate('m - 2 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax1n = plt.subplot2grid((3,6),(2,1),sharex=ax1a,sharey=ax1a)
ax1n.set_xlabel('fraction occluded (%)',fontsize=axis_size)
ax1n.annotate('n - 5 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax1o = plt.subplot2grid((3,6),(2,2),sharex=ax1a,sharey=ax1a)
ax1o.set_xlabel('fraction occluded (%)',fontsize=axis_size)
ax1o.annotate('o - 10 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax1p = plt.subplot2grid((3,6),(2,3),sharex=ax1a,sharey=ax1a)
ax1p.set_xlabel('fraction occluded (%)',fontsize=axis_size)
ax1p.annotate('p - 20 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax1q = plt.subplot2grid((3,6),(2,4),sharex=ax1a,sharey=ax1a)
ax1q.set_xlabel('fraction occluded (%)',fontsize=axis_size)
ax1q.annotate('q - 50 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax1r = plt.subplot2grid((3,6),(2,5),sharex=ax1a,sharey=ax1a)
ax1r.set_xlabel('fraction occluded (%)',fontsize=axis_size)
ax1r.annotate('r - 100 m', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

# loop through the subplots
sp = [ax1a,ax1b,ax1c,ax1d,ax1e,ax1f,ax1g,ax1h,ax1i,ax1j,ax1k,ax1l,ax1m,ax1n,ax1o,ax1p,ax1q,ax1r]
pkeys = ['2m','5m','10m','20m','50m','100m']
for pp in range(0,len(sp)):

    if pp<6:
        lim_sp = np.nanmean(penetration_lim_Belian[pkeys[pp]]['40'],axis=1)*100.
    elif pp<12:
        lim_sp = np.nanmean(penetration_lim_E[pkeys[pp-6]]['40'],axis=1)*100.
    else:
        lim_sp = np.nanmean(penetration_lim_BNorth[pkeys[pp-12]]['40'],axis=1)*100.


    lim = np.nanmedian(lim_sp,axis=0)
    limperc = np.nanpercentile(lim_sp,[2.5,97.5],axis=0)
    sp[pp].fill_betweenx(heights[2:],limperc[0][2:],limperc[1][2:],color=colour[1],alpha=0.3)
    sp[pp].plot(lim[2:],heights[2:],'-',c=colour[0],linewidth=1.5)

ax1a.set_ylim(0,80)
ax1a.set_xlim(0,110)

ax1a.locator_params(axis='x',nbins=2)
ax1b.locator_params(axis='x',nbins=2)
ax1c.locator_params(axis='x',nbins=2)
ax1d.locator_params(axis='x',nbins=2)
ax1e.locator_params(axis='x',nbins=2)
ax1f.locator_params(axis='x',nbins=2)
ax1g.locator_params(axis='x',nbins=2)
ax1h.locator_params(axis='x',nbins=2)
ax1i.locator_params(axis='x',nbins=2)
ax1j.locator_params(axis='x',nbins=2)
ax1k.locator_params(axis='x',nbins=2)
ax1l.locator_params(axis='x',nbins=2)
ax1m.locator_params(axis='x',nbins=2)
ax1n.locator_params(axis='x',nbins=2)
ax1o.locator_params(axis='x',nbins=2)
ax1p.locator_params(axis='x',nbins=2)
ax1q.locator_params(axis='x',nbins=2)
ax1r.locator_params(axis='x',nbins=2)

yticklabels = ax1b.get_yticklabels() + ax1c.get_yticklabels() + ax1d.get_yticklabels() + ax1e.get_yticklabels() + ax1f.get_yticklabels() + \
              ax1h.get_yticklabels() + ax1i.get_yticklabels() + ax1j.get_yticklabels() + ax1k.get_yticklabels() + ax1l.get_yticklabels() + \
              ax1n.get_yticklabels() + ax1o.get_yticklabels() + ax1p.get_yticklabels() + ax1q.get_yticklabels() + ax1r.get_yticklabels()
xticklabels = ax1a.get_xticklabels() + ax1b.get_xticklabels() + ax1c.get_xticklabels() + ax1d.get_xticklabels() + ax1e.get_xticklabels() + ax1f.get_xticklabels() + \
              ax1g.get_xticklabels() + ax1h.get_xticklabels() + ax1i.get_xticklabels() + ax1j.get_xticklabels() + ax1k.get_xticklabels() + ax1l.get_xticklabels()
plt.setp(yticklabels,visible=False)
plt.setp(xticklabels,visible=False)
plt.subplots_adjust(hspace=0.2, wspace = 0.)
plt.savefig('penetration_limits_sensitivity_1ha_average.png')
plt.show()


#===============================================================================
# 6) PAI vs resolution at different shot spacings for a specified plot
import seaborn as sns
import pandas as pd

res_keys = ['2m','5m','10m','20m','50m','100m']
N_res = len(res_keys)
dens_keys = ['5','20','40']
N_dens = len(dens_keys)
dens = []
res = []
PAI_MH = np.arange(0)
PAI_rad1 = np.arange(0)
PAI_rad2 = np.arange(0)
plot = []
plot_name = ['MLA01','SAF03','SAF02']
for rr in range(0,N_res):
    for dd in range(0,N_dens):
        for pp in range(0,len(plot_name)):
            if pp==0:
                PAI_MH = np.append(PAI_MH,np.sum(np.nanmean(PAD_profiles_MH_Belian[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                PAI_rad1 = np.append(PAI_rad1,np.sum(np.nanmean(PAD_profiles_rad1_Belian[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                PAI_rad2 = np.append(PAI_rad2,np.sum(np.nanmean(PAD_profiles_rad2_Belian[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                N_obs = np.sum(np.nanmean(PAD_profiles_rad1_Belian[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1).size
                for ii in range(0,N_obs):
                    dens.append(dens_keys[dd])
                    res.append(res_keys[rr])
                    plot.append(plot_name[pp])
            if pp==1:
                PAI_MH = np.append(PAI_MH,np.sum(np.nanmean(PAD_profiles_MH_E[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                PAI_rad1 = np.append(PAI_rad1,np.sum(np.nanmean(PAD_profiles_rad1_E[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                PAI_rad2 = np.append(PAI_rad2,np.sum(np.nanmean(PAD_profiles_rad2_E[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                N_obs = np.sum(np.nanmean(PAD_profiles_rad1_E[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1).size
                for ii in range(0,N_obs):
                    dens.append(dens_keys[dd])
                    res.append(res_keys[rr])
                    plot.append(plot_name[pp])
            if pp==2:
                PAI_MH = np.append(PAI_MH,np.sum(np.nanmean(PAD_profiles_MH_BNorth[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                PAI_rad1 = np.append(PAI_rad1,np.sum(np.nanmean(PAD_profiles_rad1_BNorth[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                PAI_rad2 = np.append(PAI_rad2,np.sum(np.nanmean(PAD_profiles_rad2_BNorth[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1))
                N_obs = np.sum(np.nanmean(PAD_profiles_rad1_BNorth[res_keys[rr]][dens_keys[dd]],axis=1)[:,2:],axis=1).size
                for ii in range(0,N_obs):
                    dens.append(dens_keys[dd])
                    res.append(res_keys[rr])
                    plot.append(plot_name[pp])
df = pd.DataFrame({'plot' : plot,'point density' : dens,'resolution' : res,'PAI_MH':PAI_MH,'PAI_rad1':PAI_rad1,'PAI_rad2':PAI_rad2})

pp=0
plt.figure(6, facecolor='White',figsize=[9,9])
# row one of matrix - point density is 5 pts m-2
mask = np.all(((df['plot']==plot_name[pp]),(df['point density']=='5')),axis=0)

ax6a = plt.subplot2grid((3,3),(0,0))
ax6a.set_ylabel('PAI',fontsize=axis_size)
ax6a.annotate('a - MacArthur-Horn', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
sns.violinplot(x='resolution',y='PAI_MH',data=df[mask],inner=None,linewidth=0.5,scale='width',color=colour[0])
ax6a.set_ylabel('')
ax6a.set_xlabel('')

ax6b = plt.subplot2grid((3,3),(0,1),sharex=ax6a,sharey=ax6a)
ax6b.set_title('Point density = 5 pts m$^{-2}$', fontsize=10)
ax6b.annotate('b - rad. trans. (Detto)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
sns.violinplot(x='resolution',y='PAI_rad1',data=df[mask],inner=None,linewidth=0.5,scale='width',color=colour[1])
ax6b.set_ylabel('')
ax6b.set_xlabel('')
ax6e.set_ylabel('PAI',fontsize=axis_size)

ax6c = plt.subplot2grid((3,3),(0,2),sharex=ax6a,sharey=ax6a)
ax6c.annotate('c - rad. trans (modified)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
sns.violinplot(x='resolution',y='PAI_rad2',data=df[mask],inner=None,linewidth=0.5,scale='width',color=colour[2])
ax6c.set_ylabel('')
ax6c.set_xlabel('')

# row two of matrix - point density is 20 pts m-2
mask = np.all(((df['plot']==plot_name[pp]),(df['point density']=='20')),axis=0)

ax6d = plt.subplot2grid((3,3),(1,0),sharex=ax6a,sharey=ax6a)
ax6d.annotate('d', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
sns.violinplot(x='resolution',y='PAI_MH',data=df[mask],inner=None,linewidth=0.5,scale='width',color=colour[0])
ax6d.set_ylabel('')
ax6d.set_xlabel('')
ax6d.set_ylabel('PAI',fontsize=axis_size)

ax6e = plt.subplot2grid((3,3),(1,1),sharex=ax6a,sharey=ax6a)
ax6e.set_title('Point density = 20 pts m$^{-2}$', fontsize=10)
ax6e.annotate('e', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
sns.violinplot(x='resolution',y='PAI_rad1',data=df[mask],inner=None,linewidth=0.5,scale='width',color=colour[1])
ax6e.set_ylabel('')
ax6e.set_xlabel('')

ax6f = plt.subplot2grid((3,3),(1,2),sharex=ax6a,sharey=ax6a)
ax6f.annotate('f', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
sns.violinplot(x='resolution',y='PAI_rad2',data=df[mask],inner=None,linewidth=0.5,scale='width',color=colour[2])
ax6f.set_ylabel('')
ax6f.set_xlabel('')


# row three of matrix - point density is 40 pts m-2
mask = np.all(((df['plot']==plot_name[pp]),(df['point density']=='40')),axis=0)

ax6g = plt.subplot2grid((3,3),(2,0),sharex=ax6a,sharey=ax6a)
ax6g.annotate('g', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
sns.violinplot(x='resolution',y='PAI_MH',data=df[mask],inner=None,linewidth=0.5,scale='width',color=colour[0])
ax6g.set_ylabel('')
ax6g.set_xlabel('')
ax6g.set_ylabel('PAI',fontsize=axis_size)

ax6h = plt.subplot2grid((3,3),(2,1),sharex=ax6a,sharey=ax6a)
ax6h.set_title('Point density = 40 pts m$^{-2}$', fontsize=10)
ax6h.annotate('h', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
sns.violinplot(x='resolution',y='PAI_rad1',data=df[mask],inner=None,linewidth=0.5,scale='width',color=colour[1])
ax6h.set_ylabel('')
ax6h.set_xlabel('')

ax6i = plt.subplot2grid((3,3),(2,2),sharex=ax6a,sharey=ax6a)
ax6i.annotate('i', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
sns.violinplot(x='resolution',y='PAI_rad2',data=df[mask],inner=None,linewidth=0.5,scale='width',color=colour[2])
ax6i.set_ylabel('')
ax6i.set_xlabel('')

ax6a.set_xlim(xmax=14)
plt.subplots_adjust(wspace = 0.2,hspace=0.4)
plt.savefig('PAI_sensitivity_revised.png')

plt.show()
