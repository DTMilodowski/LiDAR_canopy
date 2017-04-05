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
heights = np.arange(1.,81.)
light_absorption = {}
light_transmittance = {}
k = 0.034
n_plots = len(LAD_profiles.keys())
for pp in range(0,n_plots):
    plot = LAD_profiles.keys()[pp]
    print plot
    I = np.zeros(LAD_profiles[plot].shape)
    A = np.zeros(LAD_profiles[plot].shape)
    n_sub = I.shape[0]
    for ss in range(0,n_sub):
        I[ss,:]=clim.estimate_canopy_light_transmittance(LAD_profiles[plot][ss],heights,k)
        A[ss,:]=clim.estimate_canopy_light_absorption(I,k)

    light_transmittance[plot] = I.copy()
    light_absorption[plot] = A.copy()
