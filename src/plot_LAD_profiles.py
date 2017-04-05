# This library hosts some scripts to produce plots of LAD profiles
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rcParams

# This function produces a figure that provides a detailed look at LAD variations across a GEM ICP.
# The layout of the subplot axes represents the subplot distribution 
# Profiles_in is a 2D array (N_profiles x N_layers) for the 25 subplots in a given plot.
# Heights_in is an array giving the heights of the top of each layer
# color_string is the colour you'd like to use for your plots
# label_string is a string you use to label the plot e.g. the GEM plot name
# figure_name is the filename you'd like to use for the resultant figure
def plot_subplot_LAD_profiles(Profiles_in,Heights_in,color_string,label_string,figure_name):
    n_subplots,n_layers = Profiles_in.shape
    Profiles = np.zeros((n_subplots,n_layers*2))
    Heights = np.zeros(n_layers*2)
    for i in range(0,n_layers):
        Profiles[:,i*2] = Profiles_in[:,i]
        Profiles[:,i*2+1] = Profiles_in[:,i]
        Heights[i*2] = Heights_in[i]
        Heights[i*2+1] = Heights_in[i]
    Heights=np.roll(Heights,1)
    Heights[0] = 0
    dz = Heights_in[1]-Heights_in[0]
    Profiles[np.isfinite(Profiles)==False]=0

    plt.figure(1,facecolor='White',figsize=[15,10])
    
    ax40 = plt.subplot2grid((10,15),(0,0),colspan=2,rowspan=2)
    ax40.plot(Profiles[20,:],Heights,color=color_string)
    ax40.fill_betweenx(Heights,0,Profiles[20,:],color=color_string,alpha=0.5)

    ax41 = plt.subplot2grid((10,15),(0,2),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax41.plot(Profiles[21,:],Heights,color=color_string)
    ax41.fill_betweenx(Heights,0,Profiles[21,:],color=color_string,alpha=0.5)

    ax42 = plt.subplot2grid((10,15),(0,4),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax42.plot(Profiles[22,:],Heights,color=color_string)
    ax42.fill_betweenx(Heights,0,Profiles[22,:],color=color_string,alpha=0.5)

    ax43 = plt.subplot2grid((10,15),(0,6),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax43.plot(Profiles[23,:],Heights,color=color_string)
    ax43.fill_betweenx(Heights,0,Profiles[23,:],color=color_string,alpha=0.5)

    ax44 = plt.subplot2grid((10,15),(0,8),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax44.plot(Profiles[24,:],Heights,color=color_string)
    ax44.fill_betweenx(Heights,0,Profiles[24,:],color=color_string,alpha=0.5)

    ax30 = plt.subplot2grid((10,15),(2,0),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax30.plot(Profiles[19,:],Heights,color=color_string)
    ax30.fill_betweenx(Heights,0,Profiles[19,:],color=color_string,alpha=0.5)

    ax31 = plt.subplot2grid((10,15),(2,2),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax31.plot(Profiles[18,:],Heights,color=color_string)
    ax31.fill_betweenx(Heights,0,Profiles[18,:],color=color_string,alpha=0.5)

    ax32 = plt.subplot2grid((10,15),(2,4),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax32.plot(Profiles[17,:],Heights,color=color_string)
    ax32.fill_betweenx(Heights,0,Profiles[17,:],color=color_string,alpha=0.5)

    ax33 = plt.subplot2grid((10,15),(2,6),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax33.plot(Profiles[16,:],Heights,color=color_string)
    ax33.fill_betweenx(Heights,0,Profiles[16,:],color=color_string,alpha=0.5)

    ax34 = plt.subplot2grid((10,15),(2,8),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax34.plot(Profiles[15,:],Heights,color=color_string)
    ax34.fill_betweenx(Heights,0,Profiles[15,:],color=color_string,alpha=0.5)

    ax20 = plt.subplot2grid((10,15),(4,0),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax20.plot(Profiles[10,:],Heights,color=color_string)
    ax20.fill_betweenx(Heights,0,Profiles[10,:],color=color_string,alpha=0.5)

    ax21 = plt.subplot2grid((10,15),(4,2),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax21.plot(Profiles[11,:],Heights,color=color_string)
    ax21.fill_betweenx(Heights,0,Profiles[11,:],color=color_string,alpha=0.5)

    ax22 = plt.subplot2grid((10,15),(4,4),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax22.plot(Profiles[12,:],Heights,color=color_string)
    ax22.fill_betweenx(Heights,0,Profiles[12,:],color=color_string,alpha=0.5)

    ax23 = plt.subplot2grid((10,15),(4,6),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax23.plot(Profiles[13,:],Heights,color=color_string)
    ax23.fill_betweenx(Heights,0,Profiles[13,:],color=color_string,alpha=0.5)

    ax24 = plt.subplot2grid((10,15),(4,8),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax24.plot(Profiles[14,:],Heights,color=color_string)
    ax24.fill_betweenx(Heights,0,Profiles[14,:],color=color_string,alpha=0.5)

    ax10 = plt.subplot2grid((10,15),(6,0),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax10.plot(Profiles[9,:],Heights,color=color_string)
    ax10.fill_betweenx(Heights,0,Profiles[9,:],color=color_string,alpha=0.5)

    ax11 = plt.subplot2grid((10,15),(6,2),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax11.plot(Profiles[8,:],Heights,color=color_string)
    ax11.fill_betweenx(Heights,0,Profiles[8,:],color=color_string,alpha=0.5)

    ax12 = plt.subplot2grid((10,15),(6,4),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax12.plot(Profiles[7,:],Heights,color=color_string)
    ax12.fill_betweenx(Heights,0,Profiles[7,:],color=color_string,alpha=0.5)

    ax13 = plt.subplot2grid((10,15),(6,6),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax13.plot(Profiles[6,:],Heights,color=color_string)
    ax13.fill_betweenx(Heights,0,Profiles[6,:],color=color_string,alpha=0.5)

    ax14 = plt.subplot2grid((10,15),(6,8),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax14.plot(Profiles[5,:],Heights,color=color_string)
    ax14.fill_betweenx(Heights,0,Profiles[5,:],color=color_string,alpha=0.5)

    ax00 = plt.subplot2grid((10,15),(8,0),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax00.plot(Profiles[0,:],Heights,color=color_string)
    ax00.fill_betweenx(Heights,0,Profiles[0,:],color=color_string,alpha=0.5)

    ax01 = plt.subplot2grid((10,15),(8,2),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax01.plot(Profiles[1,:],Heights,color=color_string)
    ax01.fill_betweenx(Heights,0,Profiles[1,:],color=color_string,alpha=0.5)

    ax02 = plt.subplot2grid((10,15),(8,4),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax02.plot(Profiles[2,:],Heights,color=color_string)
    ax02.fill_betweenx(Heights,0,Profiles[2,:],color=color_string,alpha=0.5)

    ax03 = plt.subplot2grid((10,15),(8,6),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax03.plot(Profiles[3,:],Heights,color=color_string)
    ax03.fill_betweenx(Heights,0,Profiles[3,:],color=color_string,alpha=0.5)

    ax04 = plt.subplot2grid((10,15),(8,8),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax04.plot(Profiles[4,:],Heights,color=color_string)
    ax04.fill_betweenx(Heights,0,Profiles[4,:],color=color_string,alpha=0.5)

    ax40.set_ylim(ymin=0,ymax=80)
    ax40.set_xlim(xmin=0,xmax=0.8)

    xticklabels = ax00.get_xticklabels()+ax01.get_xticklabels()+ax02.get_xticklabels()+ax03.get_xticklabels()+ax04.get_xticklabels() + ax10.get_xticklabels()+ax11.get_xticklabels()+ax12.get_xticklabels()+ax13.get_xticklabels()+ax14.get_xticklabels() + ax20.get_xticklabels()+ax21.get_xticklabels()+ax22.get_xticklabels()+ax23.get_xticklabels()+ax24.get_xticklabels() + ax30.get_xticklabels()+ax31.get_xticklabels()+ax32.get_xticklabels()+ax33.get_xticklabels()+ax34.get_xticklabels() + ax40.get_xticklabels()+ax41.get_xticklabels()+ax42.get_xticklabels()+ax43.get_xticklabels()+ax44.get_xticklabels()

    yticklabels = ax00.get_yticklabels()+ax01.get_yticklabels()+ax02.get_yticklabels()+ax03.get_yticklabels()+ax04.get_yticklabels() + ax10.get_yticklabels()+ax11.get_yticklabels()+ax12.get_yticklabels()+ax13.get_yticklabels()+ax14.get_yticklabels() + ax20.get_yticklabels()+ax21.get_yticklabels()+ax22.get_yticklabels()+ax23.get_yticklabels()+ax24.get_yticklabels() + ax30.get_yticklabels()+ax31.get_yticklabels()+ax32.get_yticklabels()+ax33.get_yticklabels()+ax34.get_yticklabels() + ax40.get_yticklabels()+ax41.get_yticklabels()+ax42.get_yticklabels()+ax43.get_yticklabels()+ax44.get_yticklabels()


    plt.setp(yticklabels, visible=False)
    plt.setp(xticklabels, visible=False)
    plt.subplots_adjust(hspace=0.001,wspace=0.001)

    ax_main = plt.subplot2grid((10,15),(0,10),colspan=5,rowspan=10)
    ax_main.plot(np.mean(Profiles,axis=0),Heights,color=color_string)
    ax_main.fill_betweenx(Heights,0,np.mean(Profiles,axis=0),color=color_string,alpha=0.5)

    ax_main.set_ylim(ymin=0,ymax=80)
    ax_main.set_xlim(xmin=0,xmax=0.8)

    LAI ='%.2f' % np.sum(np.mean(Profiles,axis=0)*dz/2.)
    ax_main.annotate('LAI = ' + LAI, xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=rcParams['font.size']+2)
    ax_main.annotate(label_string, xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2)

    ax_main.set_xlabel('Leaf Area Density')
    ax_main.set_ylabel('Height / m')
    ax_main.yaxis.tick_right()
    ax_main.yaxis.set_label_position("right")
    #plt.tight_layout()
    plt.savefig(figure_name+".png",format="png")
    #plt.show()

# Equivalent function that instead plots profiles of transmittance or absorption through the canopy at both the subplot level and plot averages
def plot_subplot_transmittance_profiles(Profiles,Heights,color_string,label_string,figure_name):
    n_subplots,n_layers = Profiles.shape
    dz = Heights[1]-Heights[0]
    Profiles[np.isfinite(Profiles)==False]=0
    plt.figure(1,facecolor='White',figsize=[15,10])
    
    ax40 = plt.subplot2grid((10,15),(0,0),colspan=2,rowspan=2)
    ax40.plot(Profiles[20,:],Heights,color=color_string)

    ax41 = plt.subplot2grid((10,15),(0,2),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax41.plot(Profiles[21,:],Heights,color=color_string)

    ax42 = plt.subplot2grid((10,15),(0,4),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax42.plot(Profiles[22,:],Heights,color=color_string)

    ax43 = plt.subplot2grid((10,15),(0,6),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax43.plot(Profiles[23,:],Heights,color=color_string)

    ax44 = plt.subplot2grid((10,15),(0,8),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax44.plot(Profiles[24,:],Heights,color=color_string)

    ax30 = plt.subplot2grid((10,15),(2,0),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax30.plot(Profiles[19,:],Heights,color=color_string)

    ax31 = plt.subplot2grid((10,15),(2,2),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax31.plot(Profiles[18,:],Heights,color=color_string)

    ax32 = plt.subplot2grid((10,15),(2,4),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax32.plot(Profiles[17,:],Heights,color=color_string)

    ax33 = plt.subplot2grid((10,15),(2,6),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax33.plot(Profiles[16,:],Heights,color=color_string)

    ax34 = plt.subplot2grid((10,15),(2,8),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax34.plot(Profiles[15,:],Heights,color=color_string)

    ax20 = plt.subplot2grid((10,15),(4,0),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax20.plot(Profiles[10,:],Heights,color=color_string)

    ax21 = plt.subplot2grid((10,15),(4,2),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax21.plot(Profiles[11,:],Heights,color=color_string)

    ax22 = plt.subplot2grid((10,15),(4,4),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax22.plot(Profiles[12,:],Heights,color=color_string)

    ax23 = plt.subplot2grid((10,15),(4,6),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax23.plot(Profiles[13,:],Heights,color=color_string)

    ax24 = plt.subplot2grid((10,15),(4,8),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax24.plot(Profiles[14,:],Heights,color=color_string)

    ax10 = plt.subplot2grid((10,15),(6,0),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax10.plot(Profiles[9,:],Heights,color=color_string)

    ax11 = plt.subplot2grid((10,15),(6,2),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax11.plot(Profiles[8,:],Heights,color=color_string)

    ax12 = plt.subplot2grid((10,15),(6,4),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax12.plot(Profiles[7,:],Heights,color=color_string)

    ax13 = plt.subplot2grid((10,15),(6,6),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax13.plot(Profiles[6,:],Heights,color=color_string)

    ax14 = plt.subplot2grid((10,15),(6,8),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax14.plot(Profiles[5,:],Heights,color=color_string)

    ax00 = plt.subplot2grid((10,15),(8,0),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax00.plot(Profiles[0,:],Heights,color=color_string)

    ax01 = plt.subplot2grid((10,15),(8,2),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax01.plot(Profiles[1,:],Heights,color=color_string)

    ax02 = plt.subplot2grid((10,15),(8,4),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax02.plot(Profiles[2,:],Heights,color=color_string)

    ax03 = plt.subplot2grid((10,15),(8,6),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax03.plot(Profiles[3,:],Heights,color=color_string)

    ax04 = plt.subplot2grid((10,15),(8,8),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax04.plot(Profiles[4,:],Heights,color=color_string)

    ax40.set_ylim(ymin=0,ymax=80)

    xticklabels = ax00.get_xticklabels()+ax01.get_xticklabels()+ax02.get_xticklabels()+ax03.get_xticklabels()+ax04.get_xticklabels() + ax10.get_xticklabels()+ax11.get_xticklabels()+ax12.get_xticklabels()+ax13.get_xticklabels()+ax14.get_xticklabels() + ax20.get_xticklabels()+ax21.get_xticklabels()+ax22.get_xticklabels()+ax23.get_xticklabels()+ax24.get_xticklabels() + ax30.get_xticklabels()+ax31.get_xticklabels()+ax32.get_xticklabels()+ax33.get_xticklabels()+ax34.get_xticklabels() + ax40.get_xticklabels()+ax41.get_xticklabels()+ax42.get_xticklabels()+ax43.get_xticklabels()+ax44.get_xticklabels()

    yticklabels = ax00.get_yticklabels()+ax01.get_yticklabels()+ax02.get_yticklabels()+ax03.get_yticklabels()+ax04.get_yticklabels() + ax10.get_yticklabels()+ax11.get_yticklabels()+ax12.get_yticklabels()+ax13.get_yticklabels()+ax14.get_yticklabels() + ax20.get_yticklabels()+ax21.get_yticklabels()+ax22.get_yticklabels()+ax23.get_yticklabels()+ax24.get_yticklabels() + ax30.get_yticklabels()+ax31.get_yticklabels()+ax32.get_yticklabels()+ax33.get_yticklabels()+ax34.get_yticklabels() + ax40.get_yticklabels()+ax41.get_yticklabels()+ax42.get_yticklabels()+ax43.get_yticklabels()+ax44.get_yticklabels()


    plt.setp(yticklabels, visible=False)
    plt.setp(xticklabels, visible=False)
    plt.subplots_adjust(hspace=0.001,wspace=0.001)

    ax_main = plt.subplot2grid((10,15),(0,10),colspan=5,rowspan=10)
    ax_main.plot(np.mean(Profiles,axis=0),Heights,color=color_string)

    ax_main.set_ylim(ymin=0,ymax=80)
    
    Transmission ='%.2f' % np.mean(Profiles[:,Heights==2])*dz
    ax_main.annotate('Transmission to 2m height = ' + Transmission, xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=rcParams['font.size']+2)
    ax_main.annotate(label_string, xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2)

    ax_main.set_xlabel('Fractional Light Transmittance')
    ax_main.set_ylabel('Height / m')
    ax_main.yaxis.tick_right()
    ax_main.yaxis.set_label_position("right")
    #plt.tight_layout()
    plt.savefig(figure_name+".png",format="png")
    #plt.show()

# An equivalent function that instead plots profiles of transmittance through the canopy at both the subplot level and plot averages
def plot_subplot_absorption_profiles(Profiles,Heights,color_string,label_string,figure_name):
    print Profiles.shape
    print Heights.shape
    n_subplots,n_layers = Profiles.shape
    dz = Heights[1]-Heights[0]
    Profiles[np.isfinite(Profiles)==False]=0
    plt.figure(1,facecolor='White',figsize=[15,10])
    
    ax40 = plt.subplot2grid((10,15),(0,0),colspan=2,rowspan=2)
    ax40.plot(Profiles[20,:],Heights,color=color_string)

    ax41 = plt.subplot2grid((10,15),(0,2),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax41.plot(Profiles[21,:],Heights,color=color_string)

    ax42 = plt.subplot2grid((10,15),(0,4),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax42.plot(Profiles[22,:],Heights,color=color_string)

    ax43 = plt.subplot2grid((10,15),(0,6),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax43.plot(Profiles[23,:],Heights,color=color_string)

    ax44 = plt.subplot2grid((10,15),(0,8),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax44.plot(Profiles[24,:],Heights,color=color_string)

    ax30 = plt.subplot2grid((10,15),(2,0),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax30.plot(Profiles[19,:],Heights,color=color_string)

    ax31 = plt.subplot2grid((10,15),(2,2),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax31.plot(Profiles[18,:],Heights,color=color_string)

    ax32 = plt.subplot2grid((10,15),(2,4),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax32.plot(Profiles[17,:],Heights,color=color_string)

    ax33 = plt.subplot2grid((10,15),(2,6),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax33.plot(Profiles[16,:],Heights,color=color_string)

    ax34 = plt.subplot2grid((10,15),(2,8),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax34.plot(Profiles[15,:],Heights,color=color_string)

    ax20 = plt.subplot2grid((10,15),(4,0),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax20.plot(Profiles[10,:],Heights,color=color_string)

    ax21 = plt.subplot2grid((10,15),(4,2),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax21.plot(Profiles[11,:],Heights,color=color_string)

    ax22 = plt.subplot2grid((10,15),(4,4),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax22.plot(Profiles[12,:],Heights,color=color_string)

    ax23 = plt.subplot2grid((10,15),(4,6),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax23.plot(Profiles[13,:],Heights,color=color_string)

    ax24 = plt.subplot2grid((10,15),(4,8),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax24.plot(Profiles[14,:],Heights,color=color_string)

    ax10 = plt.subplot2grid((10,15),(6,0),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax10.plot(Profiles[9,:],Heights,color=color_string)

    ax11 = plt.subplot2grid((10,15),(6,2),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax11.plot(Profiles[8,:],Heights,color=color_string)

    ax12 = plt.subplot2grid((10,15),(6,4),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax12.plot(Profiles[7,:],Heights,color=color_string)

    ax13 = plt.subplot2grid((10,15),(6,6),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax13.plot(Profiles[6,:],Heights,color=color_string)

    ax14 = plt.subplot2grid((10,15),(6,8),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax14.plot(Profiles[5,:],Heights,color=color_string)

    ax00 = plt.subplot2grid((10,15),(8,0),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax00.plot(Profiles[0,:],Heights,color=color_string)

    ax01 = plt.subplot2grid((10,15),(8,2),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax01.plot(Profiles[1,:],Heights,color=color_string)

    ax02 = plt.subplot2grid((10,15),(8,4),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax02.plot(Profiles[2,:],Heights,color=color_string)

    ax03 = plt.subplot2grid((10,15),(8,6),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax03.plot(Profiles[3,:],Heights,color=color_string)

    ax04 = plt.subplot2grid((10,15),(8,8),colspan=2,rowspan=2,sharex=ax40,sharey=ax40)
    ax04.plot(Profiles[4,:],Heights,color=color_string)

    ax40.set_ylim(ymin=0,ymax=80)

    xticklabels = ax00.get_xticklabels()+ax01.get_xticklabels()+ax02.get_xticklabels()+ax03.get_xticklabels()+ax04.get_xticklabels() + ax10.get_xticklabels()+ax11.get_xticklabels()+ax12.get_xticklabels()+ax13.get_xticklabels()+ax14.get_xticklabels() + ax20.get_xticklabels()+ax21.get_xticklabels()+ax22.get_xticklabels()+ax23.get_xticklabels()+ax24.get_xticklabels() + ax30.get_xticklabels()+ax31.get_xticklabels()+ax32.get_xticklabels()+ax33.get_xticklabels()+ax34.get_xticklabels() + ax40.get_xticklabels()+ax41.get_xticklabels()+ax42.get_xticklabels()+ax43.get_xticklabels()+ax44.get_xticklabels()

    yticklabels = ax00.get_yticklabels()+ax01.get_yticklabels()+ax02.get_yticklabels()+ax03.get_yticklabels()+ax04.get_yticklabels() + ax10.get_yticklabels()+ax11.get_yticklabels()+ax12.get_yticklabels()+ax13.get_yticklabels()+ax14.get_yticklabels() + ax20.get_yticklabels()+ax21.get_yticklabels()+ax22.get_yticklabels()+ax23.get_yticklabels()+ax24.get_yticklabels() + ax30.get_yticklabels()+ax31.get_yticklabels()+ax32.get_yticklabels()+ax33.get_yticklabels()+ax34.get_yticklabels() + ax40.get_yticklabels()+ax41.get_yticklabels()+ax42.get_yticklabels()+ax43.get_yticklabels()+ax44.get_yticklabels()


    plt.setp(yticklabels, visible=False)
    plt.setp(xticklabels, visible=False)
    plt.subplots_adjust(hspace=0.001,wspace=0.001)

    ax_main = plt.subplot2grid((10,15),(0,10),colspan=5,rowspan=10)
    ax_main.plot(np.mean(Profiles,axis=0),Heights,color=color_string)

    ax_main.set_ylim(ymin=0,ymax=80)

    Absorption='%.2f' % np.sum(np.mean(Profiles[:,Heights>=2],axis=0)*dz)
    ax_main.annotate('Absorption to 2m height = ' + Absorption, xy=(0.95,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='right', verticalalignment='top', fontsize=rcParams['font.size']+2)
    ax_main.annotate(label_string, xy=(0.05,0.95), xycoords='axes fraction',backgroundcolor='none',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+2)

    ax_main.set_xlabel('Fractional Light Absortption')
    ax_main.set_ylabel('Height / m')
    ax_main.yaxis.tick_right()
    ax_main.yaxis.set_label_position("right")
    #plt.tight_layout()
    plt.savefig(figure_name+".png",format="png")
    #plt.show()
