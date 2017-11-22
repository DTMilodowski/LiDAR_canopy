#===============================================================================
# PAI_limitations_figures.py
# D.T.Milodowski, November 2017
#-------------------------------------------------------------------------------
# This function contains the scripts used to produce the figures in the paper:
# "Point density imposes limitations on PAI estimation from discrete return
# LiDAR"
#-------------------------------------------------------------------------------
import numpy as np

# import plotting libraries
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap, shiftgrid
import sys
sys.path.append('/home/dmilodow/DataStore_DTM/FOREST2020/EOdata/EO_data_processing/src/plot_EO_data/colormap/')

import colormaps as cmaps
plt.register_cmap(name='viridis', cmap=cmaps.viridis)
plt.register_cmap(name='plasma', cmap=cmaps.plasma)
plt.set_cmap(cmaps.viridis)

# import required LiDAR libaries
import LiDAR_MacHorn_LAD_profiles.py as MH

#-------------------------------------------------------------------------------
# Figure 1 - this figure illustrates the analytical solution presented in this
# paper describing the threshold PAI that can be detected using discrete return
# LiDAR

