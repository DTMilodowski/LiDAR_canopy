# This package hosts analyses to investigate sub-canopy microclimate based on the
# vertical distribution of leaf area within the canopy.

import numpy as np
import sys
from matplotlib import pyplot as plt
import LiDAR_tools as lidar
import auxilliary_functions as aux
import LiDAR_MacHorn_LAD_profiles as LAD1
import LiDAR_radiative_transfer_LAD_profiles as LAD2
import inventory_based_LAD_profiles as field

# Estimate distribution of light within canopy with a given LAD, assuming the
# assumption of horizontal homogeneity.
# Light transmittance at any point in a continuous vertical column of voxels to
# a canopy level, j=v, is given by:
#     I_v = I_0*exp(-k*sum(LAD[j=v:-1]))
# k is the extinction coefficient defining Beer-Lambert-like decay of light
# irradiance through the canopy
# I_0 is the incident radiation at the top of the canopy 
def estimate_canopy_light_transmittance(LAD,heights,k,I_0):
    N_levels = heights.size
    dz = np.abs(heights[1]-heights[0])
    I = I_0*np.exp( -k*np.cumsum(LAD[::-1])[::-1]*dz )

    return I
