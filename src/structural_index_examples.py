import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

import structural_metrics as structure

pi=np.pi

# First up, let's create some generic canopy profiles
h = np.arange(0,81.,step=4)
h_hr = np.arange(0,81.,step=1)
lambda_1 = 81/3.
lambda_2 = 40
lambda_3 = 16
profile_1 = (1+np.sin(pi*2*h/lambda_1)+np.sin(pi*2*h/lambda_2))*(80-h)/5.
profile_1[profile_1<0]=0

profile_2 = (1+np.cos(pi*2*h/lambda_1)+0.5*np.cos(pi*2.3*h/lambda_2))*(80-h)/5.
profile_2[profile_2<0]=0

profile_1_hr = (1+np.sin(pi*2*h_hr/lambda_1)+np.sin(pi*2*h_hr/lambda_2))*(80-h_hr)/5.
profile_1_hr[profile_1_hr<0]=0

profile_2_hr = (1+np.cos(pi*2*h_hr/lambda_1)+0.5*np.cos(pi*2.3*h_hr/lambda_2))*(80-h_hr)/5.
profile_2_hr[profile_2_hr<0]=0

#-------------------------------------------------------------------------------
# Figure 1 - Illustration of the Frechet distance
N=h.size
y, x = np.meshgrid(h,h)
z = x.copy()

for i in range(0,N):
    for j in range(0,N):
        z[i,j] = structure.euc_dist([profile_1[i],x[i,j]],[profile_2[j],y[i,j]])
dx = h[1]-h[0]
dy = h[1]-h[0]
z = z[:-1, :-1]
levels = MaxNLocator(nbins=25).tick_values(z.min(), z.max())
cmap = plt.get_cmap('YlGnBu')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

plt.figure(1, facecolor='White',figsize=[5,5])
ax = plt.subplot2grid((3,3),(1,0),rowspan=2,colspan=2)
# contours are *point* based plots, so convert our bound into point
# centers
cf = ax.contourf(x[:-1, :-1] + dx/2., y[:-1, :-1] + dy/2., z, levels=levels,
                  cmap=cmap)



axa = plt.subplot2grid((3,3),(1,2),rowspan=2)
axa.plot(profile_2,h,'-o',color='blue')
axa.fill_between(profile_2,0,h,color='blue',alpha=0.3)

axb = plt.subplot2grid((3,3),(0,0),colspan=2)
axb.plot(h,profile_1,'-o',color='blue')
axb.fill_betweenx(profile_1,0,h,color='blue',alpha=0.7)


ax.set_xticklabels([])
ax.set_yticklabels([])
axa.set_xticklabels([])
axa.set_yticklabels([])
axb.set_xticklabels([])
axb.set_yticklabels([])

ax.set_xlabel('position along profile b')
ax.set_ylabel('position along profile a')
axa.annotate('b', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)
axb.annotate('a', xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=10)

plt.tight_layout()
plt.savefig('Figures/Frechet_example.png')
plt.show()

#-------------------------------------------------------------------------------
# Figure 2 - Illustration of the VSI
h1,p1= structure.find_maxima(h_hr,profile_1_hr)
h2,p2 = structure.find_maxima(h_hr,profile_2_hr)

plt.figure(2, facecolor='White',figsize=[4,4])
ax2 = plt.subplot2grid((1,1),(0,0))
ax2.plot(profile_1_hr,h_hr,'-',color='blue')
ax2.plot(profile_2_hr,h_hr,'-',color='red')
ax2.plot(p1,h1,'o',color='blue')
ax2.plot(p2,h2,'o',color='red')
for i in range(0,p1.size):
    ax2.plot((0,p1[i]),(h1[i],h1[i]),'--',color='blue')
for i in range(0,p2.size):
    ax2.plot((0,p2[i]),(h2[i],h2[i]),'--',color='red')

ax2.set_xticklabels([])
ax2.set_yticklabels([])
ax2.set_ylabel('height')
ax2.set_xlabel('LAD')

plt.tight_layout()
plt.savefig('Figures/VSI_example.png')
plt.show()
