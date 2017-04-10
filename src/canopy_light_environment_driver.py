import numpy as np
import sys
from matplotlib import pyplot as plt
import LiDAR_tools as lidar
import auxilliary_functions as aux
import LiDAR_MacHorn_LAD_profiles as LAD1
import LiDAR_radiative_transfer_LAD_profiles as LAD2
import structural_metrics as structure
import plot_LAD_profiles as plot_LAD
import canopy_microclimate as clim

sys.path.append('/home/dmilodow/DataStore_DTM/BALI/SPA_BALI_data_and_analysis/scripts/field_data/')
import load_field_data as field

# define the profiles' location
LAD_file = './output/BALI_subplot_LAD_profiles_MacHorn_1m.npz'
LAD_profiles = np.load(LAD_file)
Plots = LAD_profiles.keys()
heights = np.arange(1.,81.)
light_absorption = {}
light_transmittance = {}
k = 0.7
n_plots = len(Plots)
color_string = 'blue'
label_string = '-'
for pp in range(0,n_plots):
    plot = Plots[pp]
    print plot
    I = np.zeros(LAD_profiles[plot].shape)
    A = np.zeros(LAD_profiles[plot].shape)
    n_sub = I.shape[0]
    for ss in range(0,n_sub):
        print '======'
        I[ss,:]=clim.estimate_canopy_light_transmittance(LAD_profiles[plot][ss],heights,k)
        A[ss,:]=clim.estimate_canopy_light_absorption(I[ss,:],k)

    light_transmittance[plot] = I.copy()
    light_absorption[plot] = A.copy()

    figure_name = output_dir + '/light_environment/'+Plots[pp]+'_subplot_transmitance'
    plot_LAD.plot_subplot_transmittance_profiles(light_transmittance[plot],heights,color_string,label_string,figure_name)
    figure_name = output_dir + '/light_environment/'+Plots[pp]+'_subplot_absorption'
    plot_LAD.plot_subplot_absorption_profiles(light_absorption[plot],heights,color_string,label_string,figure_name)

OutFile = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/src/output/BALI_subplot_lighttransmittance'
np.savez(OutFile+'.npz', **light_transmittance)
OutFile = '/home/dmilodow/DataStore_DTM/BALI/LiDAR/src/output/BALI_subplot_light_absorption'
np.savez(OutFile+'.npz', **light_absorption)



