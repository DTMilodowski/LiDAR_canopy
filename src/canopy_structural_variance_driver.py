import numpy as np
from scipy import stats
import sys
from matplotlib import pyplot as plt
import fiona
import LiDAR_io as io
import LiDAR_tools as lidar
import auxilliary_functions as aux
import LiDAR_MacHorn_LAD_profiles as LAD1

from matplotlib import rcParams


# Set up some basiic parameters for the plots
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = 8
rcParams['legend.numpoints'] = 1
axis_size = rcParams['font.size']+2

# Some files
las_list = 'lastile_list.csv' ## CHANGE AS REQUIRED
pts_shp = 'sample_points.shp' ## CHANGE AS REQUIRED

# Some parameters
radius = 10.

# Open shapefile and get points & bounding box for ROI
shapefile= fiona.open(pts_shp)
N_pts = len(list(shapefile))
pts = np.zeros((N_pts,2))
for pp in range(0,N_pts):
    pts[pp,:]=np.asarray(list(shapefile)[pp]['geometry']['coordinates'])

W,S,E,N = shapefile.bounds
W-=radius
S-=radius
E+=radius
N+=radius

# Read in LiDAR points for region of interest
# - first find the tiles
