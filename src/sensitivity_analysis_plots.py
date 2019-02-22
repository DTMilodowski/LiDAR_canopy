import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rcParams
import matplotlib.ticker as mticker
import matplotlib.ticker as plticker
from matplotlib import colorbar
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = 8
rcParams['legend.numpoints'] = 1
axis_size = rcParams['font.size']+2
colour = ['#46E900','#1A2BCE','#E0007F']
n_subplots=25
