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
plt.figure(1, facecolor='White',figsize=[9,9])
ax1a = plt.subplot2grid((3,6),(0,0))
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
kk =0
for pp in range(0,len(sp)):

    rad1_sp = np.nanmean(PAD_profiles_rad1_Belian[pkeys[kk]]['40'],axis=1)
    rad1 = np.nanmedian(rad1_sp,axis=0)
    rad1perc = np.nanpercentile(rad1_sp,[2.5,25,75,95.5],axis=0)
    sp[pp].fill_betweenx(heights[2:],rad1perc[0][2:],rad1perc[3][2:],color=colour[1],alpha=0.3)
    sp[pp].fill_betweenx(heights[2:],rad1perc[1][2:],rad1perc[2][2:],color=colour[1],alpha=0.3)

    rad2_sp = np.nanmean(PAD_profiles_rad2_Belian[pkeys[kk]]['40'],axis=1)
    rad2 = np.nanmedian(rad2_sp,axis=0)
    rad2perc = np.nanpercentile(rad2_sp,[2.5,25,75,95.5],axis=0)
    sp[pp].fill_betweenx(heights[2:],rad2perc[0][2:],rad2perc[3][2:],color=colour[2],alpha=0.3)
    sp[pp].fill_betweenx(heights[2:],rad2perc[1][2:],rad2perc[2][2:],color=colour[2],alpha=0.3)

    MH_sp = np.nanmean(PAD_profiles_MH_Belian[pkeys[kk]]['40'],axis=1)
    MH = np.nanmedian(MH_sp,axis=0)
    MHperc = np.nanpercentile(MH_sp,[2.5,25,75,95.5],axis=0)
    sp[pp].fill_betweenx(heights[2:],MHperc[0][2:],MHperc[3][2:],color=colour[0],alpha=0.3)
    sp[pp].fill_betweenx(heights[2:],MHperc[1][2:],MHperc[2][2:],color=colour[0],alpha=0.3)

    if pp==2:
        sp[pp].plot(MH[2:],heights[2:],'-',c=colour[0],linewidth=1.5, label = 'MacArthur-Horn')
        sp[pp].plot(rad1[2:],heights[2:],'-',c=colour[1],linewidth=1.5, label = 'rad trans (Detto)')
        sp[pp].plot(rad2[2:],heights[2:],'-',c=colour[2],linewidth=1.5, label = 'rad trans (modified)')
        sp[pp].plot(rad1[2:],heights[2:],'-',c=colour[1],linewidth=1.5)
        sp[pp].plot(MH[2:],heights[2:],'-',c=colour[0],linewidth=1.5)
        sp[pp].legend(loc=9, bbox_to_anchor=(0.5, -0.2), ncol=3)
    else:
        sp[pp].plot(MH[2:],heights[2:],'-',c=colour[0],linewidth=1.5)
        sp[pp].plot(rad1[2:],heights[2:],'-',c=colour[1],linewidth=1.5)
        sp[pp].plot(rad2[2:],heights[2:],'-',c=colour[2],linewidth=1.5)

    kk+=1
    if kk==length(pkeys):
        k=0




"""
# 1b) the widths of the 50% (solid) and 95% (dashed) confidence intervals of each vertical layer
# this gives an indication of the reliability of the profile at a given canopy depth
ax1f = plt.subplot2grid((2,5),(1,0),sharey=ax1a)
ax1f.set_xlabel('CI / %',fontsize=axis_size)
ax1f.set_ylabel('height / m',fontsize=axis_size)
ax1f.annotate('f - 5 m', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)

ax1g = plt.subplot2grid((2,5),(1,1),sharex=ax1f,sharey=ax1a)
ax1g.set_xlabel('CI / %',fontsize=axis_size)
ax1g.annotate('g - 10 m', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)

ax1h = plt.subplot2grid((2,5),(1,2),sharex=ax1f,sharey=ax1a)
ax1h.set_xlabel('CI / %',fontsize=axis_size)
ax1h.annotate('h - 20 m', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)

ax1i = plt.subplot2grid((2,5),(1,3),sharex=ax1f,sharey=ax1a)
ax1i.set_xlabel('CI / %',fontsize=axis_size)
ax1i.annotate('i - 50 m', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)

ax1j = plt.subplot2grid((2,5),(1,4),sharex=ax1f,sharey=ax1a)
ax1j.set_xlabel('CI / %',fontsize=axis_size)
ax1j.annotate('j - 100 m', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)

sp = [ax1f,ax1g,ax1h,ax1i,ax1j]
for pp in range(0,len(sp)):
  MHl = np.percentile(PAD_profiles_MH[pkeys[pp]]['40'],2.5,axis=0)
  MH25 = np.percentile(PAD_profiles_MH[pkeys[pp]]['40'],25,axis=0)
  MH75 = np.percentile(PAD_profiles_MH[pkeys[pp]]['40'],75,axis=0)
  MHu = np.percentile(PAD_profiles_MH[pkeys[pp]]['40'],97.5,axis=0)
  MH_95CI_i = ((MHu-MHl)/np.mean(PAD_profiles_MH[pkeys[pp]]['40'],axis=0))*100.
  MH_50CI_i = ((MH75-MH25)/np.mean(PAD_profiles_MH[pkeys[pp]]['40'],axis=0))*100.
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

  rad1l = np.percentile(PAD_profiles_rad1[pkeys[pp]]['40'],2.5,axis=0)
  rad125 = np.percentile(PAD_profiles_rad1[pkeys[pp]]['40'],25,axis=0)
  rad175 = np.percentile(PAD_profiles_rad1[pkeys[pp]]['40'],75,axis=0)
  rad1u = np.percentile(PAD_profiles_rad1[pkeys[pp]]['40'],97.5,axis=0)
  rad1_95CI_i = ((rad1u-rad1l)/np.mean(PAD_profiles_rad1[pkeys[pp]]['40'],axis=0))*100.
  rad1_50CI_i = ((rad175-rad125)/np.mean(PAD_profiles_rad1[pkeys[pp]]['40'],axis=0))*100.
  rad1_95CI = np.nansum(rad1_95CI_i,axis=0)/np.sum(np.isfinite(rad1_95CI_i),axis=0)
  rad1_50CI = np.nansum(rad1_50CI_i,axis=0)/np.sum(np.isfinite(rad1_50CI_i),axis=0)
  if pp==2:
      sp[pp].plot(rad1_95CI[2:],heights[2:],':',c=colour[1],linewidth=1,label='rad. trans. (Detto) 95% CI')
      sp[pp].plot(rad1_50CI[2:],heights[2:],'-',c=colour[1],linewidth=1,label='rad. trans. (Detto) 50% CI')
  else:
      sp[pp].plot(rad1_95CI[2:],heights[2:],':',c=colour[1],linewidth=1)
      sp[pp].plot(rad1_50CI[2:],heights[2:],'-',c=colour[1],linewidth=1)

  rad2l = np.percentile(PAD_profiles_rad2[pkeys[pp]]['40'],2.5,axis=0)
  rad225 = np.percentile(PAD_profiles_rad2[pkeys[pp]]['40'],25,axis=0)
  rad275 = np.percentile(PAD_profiles_rad2[pkeys[pp]]['40'],75,axis=0)
  rad2u = np.percentile(PAD_profiles_rad2[pkeys[pp]]['40'],97.5,axis=0)
  rad2_95CI_i = ((rad2u-rad2l)/np.mean(PAD_profiles_rad2[pkeys[pp]]['40'],axis=0))*100.
  rad2_50CI_i = ((rad275-rad225)/np.mean(PAD_profiles_rad2[pkeys[pp]]['40'],axis=0))*100.
  rad2_95CI = np.nansum(rad2_95CI_i,axis=0)/np.sum(np.isfinite(rad2_95CI_i),axis=0)
  rad2_50CI = np.nansum(rad2_50CI_i,axis=0)/np.sum(np.isfinite(rad2_50CI_i),axis=0)
  if pp==2:
      sp[pp].plot(rad2_95CI[2:],heights[2:],':',c=colour[2],linewidth=1,label='rad. trans. (mod) 95% CI')
      sp[pp].plot(rad2_50CI[2:],heights[2:],'-',c=colour[2],linewidth=1,label='rad. trans. (mod) 50% CI')
      lgd = sp[pp].legend(loc=9, bbox_to_anchor=(0.5, -0.2), ncol=3)#(loc='upper right')
  else:
      sp[pp].plot(rad2_95CI[2:],heights[2:],':',c=colour[2],linewidth=1)
      sp[pp].plot(rad2_50CI[2:],heights[2:],'-',c=colour[2],linewidth=1)

ax1a.set_ylim(0,80)
ax1a.set_xlim(0,0.7)
ax1f.set_xlim(0,520)
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

yticklabels = ax1b.get_yticklabels() + ax1c.get_yticklabels() + ax1d.get_yticklabels() + ax1e.get_yticklabels() + ax1g.get_yticklabels() + ax1h.get_yticklabels() + ax1i.get_yticklabels() + ax1j.get_yticklabels()
plt.setp(yticklabels,visible=False)
plt.subplots_adjust(hspace=0.4, wspace = 0.1, bottom = 0.2)
plt.savefig('PAD_resolution_sensitivity.png')
plt.show()
"""

#==============================================================================================
# 2) Equivalent plots but this time using different point densities at target resolution
pkeys=['5','10','15','20','25','30','40']
plt.figure(2, facecolor='White',figsize=[9,9])
ax2a = plt.subplot2grid((2,5),(0,0))
ax2a.set_xlabel('PAD / m$^2$m$^{-3}$',fontsize=axis_size)
ax2a.set_ylabel('height / m',fontsize=axis_size)
ax2a.annotate('a - 5 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax2b = plt.subplot2grid((2,5),(0,1),sharex=ax2a,sharey=ax2a)
ax2b.set_xlabel('PAD / m$^2$m$^{-3}$',fontsize=axis_size)
ax2b.annotate('b - 10 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax2c = plt.subplot2grid((2,5),(0,2),sharex=ax2a,sharey=ax2a)
ax2c.set_xlabel('PAD / m$^2$m$^{-3}$',fontsize=axis_size)
ax2c.annotate('c - 20 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax2d = plt.subplot2grid((2,5),(0,3),sharex=ax2a,sharey=ax2a)
ax2d.set_xlabel('PAD / m$^2$m$^{-3}$',fontsize=axis_size)
ax2d.annotate('d - 30 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax2e = plt.subplot2grid((2,5),(0,4),sharex=ax2a,sharey=ax2a)
ax2e.set_xlabel('PAD / m$^2$m$^{-3}$',fontsize=axis_size)
ax2e.annotate('e - 40 pts m$^{-2}$', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

# loop through the subplots
sp = [ax2a,ax2b,ax2c,ax2d,ax2e]
pkeys = ['5','10','20','30','40']
for pp in range(0,len(sp)):
  rad2_sp = np.mean(PAD_profiles_rad2['20m'][pkeys[pp]],axis=1)
  rad2 = np.median(rad2_sp,axis=0)
  rad2perc = np.percentile(rad2_sp,[2.5,25,75,95.5],axis=0)
  sp[pp].fill_betweenx(heights[2:],rad2perc[0][2:],rad2perc[3][2:],color=colour[2],alpha=0.3)
  sp[pp].fill_betweenx(heights[2:],rad2perc[1][2:],rad2perc[2][2:],color=colour[2],alpha=0.3)

  rad1_sp = np.mean(PAD_profiles_rad1['20m'][pkeys[pp]],axis=1)
  rad1 = np.median(rad1_sp,axis=0)
  rad1perc = np.percentile(rad1_sp,[2.5,25,75,95.5],axis=0)
  sp[pp].fill_betweenx(heights[2:],rad1perc[0][2:],rad1perc[3][2:],color=colour[1],alpha=0.3)
  sp[pp].fill_betweenx(heights[2:],rad1perc[1][2:],rad1perc[2][2:],color=colour[1],alpha=0.3)

  MH_sp = np.mean(PAD_profiles_MH['20m'][pkeys[pp]],axis=1)
  MH = np.median(MH_sp,axis=0)
  MHperc = np.percentile(MH_sp,[2.5,25,75,95.5],axis=0)
  sp[pp].fill_betweenx(heights[2:],MHperc[0][2:],MHperc[3][2:],color=colour[0],alpha=0.3)
  sp[pp].fill_betweenx(heights[2:],MHperc[1][2:],MHperc[2][2:],color=colour[0],alpha=0.3)

  if pp==2:
    sp[pp].plot(MH[2:],heights[2:],'-',c=colour[0],linewidth=1.5, label = 'MacArthur-Horn')
    sp[pp].plot(rad1[2:],heights[2:],'-',c=colour[1],linewidth=1.5, label = 'rad trans (Detto)')
    sp[pp].plot(rad2[2:],heights[2:],'-',c=colour[2],linewidth=1.5, label = 'rad trans (modified)')
    sp[pp].plot(rad1[2:],heights[2:],'-',c=colour[1],linewidth=1.5)
    sp[pp].plot(MH[2:],heights[2:],'-',c=colour[0],linewidth=1.5)
    sp[pp].legend(loc=9, bbox_to_anchor=(0.5, -0.2), ncol=3)#(loc='upper right')
  else:
    sp[pp].plot(MH[2:],heights[2:],'-',c=colour[0],linewidth=1.5)
    sp[pp].plot(rad1[2:],heights[2:],'-',c=colour[1],linewidth=1.5)
    sp[pp].plot(rad2[2:],heights[2:],'-',c=colour[2],linewidth=1.5)

# 1b) the widths of the 50% (solid) and 95% (dashed) confidence intervals of each vertical layer
# this gives an indication of the reliability of the profile at a given canopy depth
ax2f = plt.subplot2grid((2,5),(1,0),sharey=ax2a)
ax2f.set_xlabel('CI / %',fontsize=axis_size)
ax2f.set_ylabel('height / m',fontsize=axis_size)
ax2f.annotate('f - 5 pts m$^{-2}$', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)

ax2g = plt.subplot2grid((2,5),(1,1),sharex=ax2f,sharey=ax2a)
ax2g.set_xlabel('CI / %',fontsize=axis_size)
ax2g.annotate('g - 10 pts m$^{-2}$', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)

ax2h = plt.subplot2grid((2,5),(1,2),sharex=ax2f,sharey=ax2a)
ax2h.set_xlabel('CI / %',fontsize=axis_size)
ax2h.annotate('h - 20 pts m$^{-2}$', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)

ax2i = plt.subplot2grid((2,5),(1,3),sharex=ax2f,sharey=ax2a)
ax2i.set_xlabel('CI / %',fontsize=axis_size)
ax2i.annotate('i - 30 pts m$^{-2}$', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)

ax2j = plt.subplot2grid((2,5),(1,4),sharex=ax2f,sharey=ax2a)
ax2j.set_xlabel('CI / %',fontsize=axis_size)
ax2j.annotate('j - 40 pts m$^{-2}$', xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=10)

sp = [ax2f,ax2g,ax2h,ax2i,ax2j]
for pp in range(0,len(sp)):
    MHl = np.percentile(PAD_profiles_MH['20m'][pkeys[pp]],2.5,axis=0)
    MH25 = np.percentile(PAD_profiles_MH['20m'][pkeys[pp]],25,axis=0)
    MH75 = np.percentile(PAD_profiles_MH['20m'][pkeys[pp]],75,axis=0)
    MHu = np.percentile(PAD_profiles_MH['20m'][pkeys[pp]],97.5,axis=0)
    MH_95CI_i = ((MHu-MHl)/np.mean(PAD_profiles_MH['20m'][pkeys[pp]],axis=0))*100.
    MH_50CI_i = ((MH75-MH25)/np.mean(PAD_profiles_MH['20m'][pkeys[pp]],axis=0))*100.
    MH_95CI = np.nansum(MH_95CI_i,axis=0)/np.sum(np.isfinite(MH_95CI_i),axis=0)
    MH_50CI = np.nansum(MH_50CI_i,axis=0)/np.sum(np.isfinite(MH_50CI_i),axis=0)
    if pp==2:
        sp[pp].plot(MH_95CI[2:],heights[2:],':',c=colour[0],linewidth=1,label='MacArthur-Horn 95% CI')
        sp[pp].plot(MH_50CI[2:],heights[2:],'-',c=colour[0],linewidth=1,label='MacArthur-Horn 50% CI')
    else:
        sp[pp].plot(MH_95CI[2:],heights[2:],':',c=colour[0],linewidth=1)
        sp[pp].plot(MH_50CI[2:],heights[2:],'-',c=colour[0],linewidth=1)

    rad1l = np.percentile(PAD_profiles_rad1['20m'][pkeys[pp]],2.5,axis=0)
    rad125 = np.percentile(PAD_profiles_rad1['20m'][pkeys[pp]],25,axis=0)
    rad175 = np.percentile(PAD_profiles_rad1['20m'][pkeys[pp]],75,axis=0)
    rad1u = np.percentile(PAD_profiles_rad1['20m'][pkeys[pp]],97.5,axis=0)
    rad1_95CI_i = ((rad1u-rad1l)/np.mean(PAD_profiles_rad1['20m'][pkeys[pp]],axis=0))*100.
    rad1_50CI_i = ((rad175-rad125)/np.mean(PAD_profiles_rad1['20m'][pkeys[pp]],axis=0))*100.
    rad1_95CI = np.nansum(rad1_95CI_i,axis=0)/np.sum(np.isfinite(rad1_95CI_i),axis=0)
    rad1_50CI = np.nansum(rad1_50CI_i,axis=0)/np.sum(np.isfinite(rad1_50CI_i),axis=0)
    if pp==2:
        sp[pp].plot(rad1_95CI[2:],heights[2:],':',c=colour[1],linewidth=1,label='rad. trans. (Detto) 95% CI')
        sp[pp].plot(rad1_50CI[2:],heights[2:],'-',c=colour[1],linewidth=1,label='rad. trans. (Detto) 50% CI')
    else:
        sp[pp].plot(rad1_95CI[2:],heights[2:],':',c=colour[1],linewidth=1)
        sp[pp].plot(rad1_50CI[2:],heights[2:],'-',c=colour[1],linewidth=1)

    rad2l = np.percentile(PAD_profiles_rad2['20m'][pkeys[pp]],2.5,axis=0)
    rad225 = np.percentile(PAD_profiles_rad2['20m'][pkeys[pp]],25,axis=0)
    rad275 = np.percentile(PAD_profiles_rad2['20m'][pkeys[pp]],75,axis=0)
    rad2u = np.percentile(PAD_profiles_rad2['20m'][pkeys[pp]],97.5,axis=0)
    rad2_95CI_i = ((rad2u-rad2l)/np.mean(PAD_profiles_rad2['20m'][pkeys[pp]],axis=0))*100.
    rad2_50CI_i = ((rad275-rad225)/np.mean(PAD_profiles_rad2['20m'][pkeys[pp]],axis=0))*100.
    rad2_95CI = np.nansum(rad2_95CI_i,axis=0)/np.sum(np.isfinite(rad2_95CI_i),axis=0)
    rad2_50CI = np.nansum(rad2_50CI_i,axis=0)/np.sum(np.isfinite(rad2_50CI_i),axis=0)
    if pp==2:
        sp[pp].plot(rad2_95CI[2:],heights[2:],':',c=colour[2],linewidth=1,label='rad. trans. (mod) 95% CI')
        sp[pp].plot(rad2_50CI[2:],heights[2:],'-',c=colour[2],linewidth=1,label='rad. trans. (mod) 50% CI')
        sp[pp].legend(loc=9, bbox_to_anchor=(0.5, -0.2), ncol=3)#(loc='upper right')
    else:
        sp[pp].plot(rad2_95CI[2:],heights[2:],':',c=colour[2],linewidth=1)
        sp[pp].plot(rad2_50CI[2:],heights[2:],'-',c=colour[2],linewidth=1)

ax2a.set_ylim(0,80)
ax2a.set_xlim(0,0.7)
ax2f.set_xlim(0,520)
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

yticklabels = ax2b.get_yticklabels() + ax2c.get_yticklabels() + ax2d.get_yticklabels() + ax2e.get_yticklabels() + ax2g.get_yticklabels() + ax2h.get_yticklabels() + ax2i.get_yticklabels() + ax2j.get_yticklabels()
plt.setp(yticklabels,visible=False)
plt.subplots_adjust(hspace=0.4, wspace = 0.1, bottom = 0.2)
plt.savefig('PAD_point_density_sensitivity.png')
plt.show()


# 3) PAI vs resolution at different shot spacings
N_res = len(keys)
plt.figure(3, facecolor='White',figsize=[9,9])
ax3a = plt.subplot2grid((3,3),(0,0))
ax3a.set_ylabel('PAI',fontsize=axis_size)
ax3a.annotate('a - MacArthur-Horn', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax3b = plt.subplot2grid((3,3),(0,1),sharex=ax3a)
ax3b.set_title('Point density = 5 pts m$^{-2}$', fontsize=10)
ax3b.annotate('b - rad. trans. (Detto)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax3c = plt.subplot2grid((3,3),(0,2),sharex=ax3a)
ax3c.annotate('c - rad. trans (modified)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

boxplot_a=[]
boxplot_b=[]
boxplot_c=[]
for kk in range(0,N_res):
    PAI_MH= np.sum(np.mean(PAD_profiles_MH[keys[kk]]['5'],axis=1)[:,2:],axis=1)
    PAI_rad1= np.sum(np.mean(PAD_profiles_rad1[keys[kk]]['5'],axis=1)[:,2:],axis=1)
    PAI_rad2= np.sum(np.mean(PAD_profiles_rad2[keys[kk]]['5'],axis=1)[:,2:],axis=1)
    boxplot_a.append(PAI_MH.copy())
    boxplot_b.append(PAI_rad1.copy())
    boxplot_c.append(PAI_rad2.copy())

bpa = ax3a.boxplot(boxplot_a, patch_artist=True)
bpb = ax3b.boxplot(boxplot_b, patch_artist=True)
bpc = ax3c.boxplot(boxplot_c, patch_artist=True)

bps=[bpa,bpb,bpc]
edge_colour = ['#30A000','#0B1681','#960055']
flier_colour = ['#85EE58','#636FDA','#E756A8']
for bb in range(0,3):
    bp = bps[bb]
    for box in bp['boxes']:
        box.set( color=edge_colour[bb])
        box.set( facecolor = colour[bb] )
    for whisker in bp['whiskers']:
        whisker.set(color=edge_colour[bb])
    for cap in bp['caps']:
        cap.set(color=edge_colour[bb])
    for median in bp['medians']:
        median.set(color='white')
    for flier in bp['fliers']:
        flier.set(marker='o', color=flier_colour[bb], alpha=0.5)



ax3d = plt.subplot2grid((3,3),(1,0),sharex=ax3a,sharey=ax3a)
ax3d.set_ylabel('PAI',fontsize=axis_size)
ax3d.annotate('d - MacArthur-Horn', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax3e = plt.subplot2grid((3,3),(1,1),sharex=ax3a,sharey=ax3b)
ax3e.set_title('Point density = 20 pts m$^{-2}$', fontsize=10)
ax3e.annotate('e - rad. trans. (Detto)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

ax3f = plt.subplot2grid((3,3),(1,2),sharex=ax3a,sharey=ax3c)
ax3f.annotate('f - rad. trans (modified)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

boxplot_d=[]
boxplot_e=[]
boxplot_f=[]
for kk in range(0,N_res):

    PAI_MH= np.sum(np.mean(PAD_profiles_MH[keys[kk]]['20'],axis=1)[:,2:],axis=1)
    PAI_rad1= np.sum(np.mean(PAD_profiles_rad1[keys[kk]]['20'],axis=1)[:,2:],axis=1)
    PAI_rad2= np.sum(np.mean(PAD_profiles_rad2[keys[kk]]['20'],axis=1)[:,2:],axis=1)
    boxplot_d.append(PAI_MH.copy())
    boxplot_e.append(PAI_rad1.copy())
    boxplot_f.append(PAI_rad2.copy())

bpd = ax3d.boxplot(boxplot_d, patch_artist=True)
bpe = ax3e.boxplot(boxplot_e, patch_artist=True)
bpf = ax3f.boxplot(boxplot_f, patch_artist=True)

bps=[bpd,bpe,bpf]
edge_colour = ['#30A000','#0B1681','#960055']
flier_colour = ['#85EE58','#636FDA','#E756A8']
for bb in range(0,3):
    bp = bps[bb]
    for box in bp['boxes']:
        box.set( color=edge_colour[bb])
        box.set( facecolor = colour[bb] )
    for whisker in bp['whiskers']:
        whisker.set(color=edge_colour[bb])
    for cap in bp['caps']:
        cap.set(color=edge_colour[bb])
    for median in bp['medians']:
        median.set(color='white')
    for flier in bp['fliers']:
        flier.set(marker='o', color=flier_colour[bb], alpha=0.5)

ax3g = plt.subplot2grid((3,3),(2,0),sharex=ax3a,sharey=ax3a)
ax3g.set_ylabel('PAI',fontsize=axis_size)
ax3g.annotate('g - MacArthur-Horn', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax3g.set_xlabel('Spatial resolution / m$^2$',fontsize=axis_size)

ax3h = plt.subplot2grid((3,3),(2,1),sharex=ax3a,sharey=ax3b)
ax3h.set_title('Point density = 40 pts m$^{-2}$', fontsize=10)
ax3h.annotate('h - rad. trans. (Detto)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax3h.set_xlabel('Spatial resolution / m$^2$',fontsize=axis_size)

ax3i = plt.subplot2grid((3,3),(2,2),sharex=ax3a,sharey=ax3c)
ax3i.annotate('i - rad. trans (modified)', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
ax3i.set_xlabel('Spatial resolution / m$^2$',fontsize=axis_size)

boxplot_g=[]
boxplot_h=[]
boxplot_i=[]
for kk in range(0,N_res):

    PAI_MH= np.sum(np.mean(PAD_profiles_MH[keys[kk]]['40'],axis=1)[:,2:],axis=1)
    PAI_rad1= np.sum(np.mean(PAD_profiles_rad1[keys[kk]]['40'],axis=1)[:,2:],axis=1)
    PAI_rad2= np.sum(np.mean(PAD_profiles_rad2[keys[kk]]['40'],axis=1)[:,2:],axis=1)
    boxplot_g.append(PAI_MH.copy())
    boxplot_h.append(PAI_rad1.copy())
    boxplot_i.append(PAI_rad2.copy())

bpg = ax3g.boxplot(boxplot_g, patch_artist=True)
bph = ax3h.boxplot(boxplot_h, patch_artist=True)
bpi = ax3i.boxplot(boxplot_i, patch_artist=True)

bps=[bpg,bph,bpi]
edge_colour = ['#30A000','#0B1681','#960055']
flier_colour = ['#85EE58','#636FDA','#E756A8']
for bb in range(0,3):
    bp = bps[bb]
    for box in bp['boxes']:
        box.set( color=edge_colour[bb])
        box.set( facecolor = colour[bb] )
    for whisker in bp['whiskers']:
        whisker.set(color=edge_colour[bb])
    for cap in bp['caps']:
        cap.set(color=edge_colour[bb])
    for median in bp['medians']:
        median.set(color='white')
    for flier in bp['fliers']:
        flier.set(marker='o', color=flier_colour[bb], alpha=0.5)

x_locs = [1,2,3,4,5,6]
ax3c.set_xticks(x_locs)
xticks=ax3c.get_xticks().tolist()
xticks=keys
ax3c.set_xticklabels(xticks,fontsize=axis_size)

xticklabels = ax3a.get_xticklabels() + ax3b.get_xticklabels() + ax3c.get_xticklabels() + ax3d.get_xticklabels() + ax3e.get_xticklabels() + ax3f.get_xticklabels()

yticklabels = ax3h.get_yticklabels() + ax3b.get_yticklabels() + ax3c.get_yticklabels() + ax3i.get_yticklabels() + ax3e.get_yticklabels() + ax3f.get_yticklabels()
ax3a.set_ylim(5,11)
ax3b.set_ylim(2,10)
ax3c.set_ylim(6,20)

plt.setp(xticklabels,visible=False)
plt.subplots_adjust(wspace = 0.4,hspace=0.2)
plt.savefig('PAI_sensitivity.png')

plt.show()
