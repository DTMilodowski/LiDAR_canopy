# LiDAR_canopy
**David T. Milodowski**, School of GeoSciences, University of Edinburgh
This repository hosts a python library to extract quantitative information about forest canopies from discrete return LiDAR point clouds.

## What is in the repository?
This repository hosts code currently under development.  Consequently, there are many different modules, many of which are surplus to requirements for calculating basic metrics.  For now, I won't give a comprehensive list of the enclosed routines, but rather provide pointers for intrepid explorers to start navigating their way through the core modules.

- LiDAR_io.py : Interfacing with laspy for i/o and spatial filtering of LiDAR returns based on neighbourhood or polygon boundaries; also constructs kdtrees
- LiDAR_tools.py : a bunch of useful tools to interact with LiDAR in python includes spatial filtering of LiDAR returns based on neighbourhood or polygon boundaries.
- auxilliary_functions.py : some generic useful stuff
- LiDAR_MacHorn_LAD_profiles.py : this code calculates the vertical leaf area distribution (_sensu stricto_ plant area distribution) using the method proposed by Stark et al. [Ecology Letters, 2012], which utilises the MacArthur-Horn-based approach to invert the vertical distribution of first returns
- LiDAR_radiative_transfer_LAD_profiles.py : this code calculates the vertical leaf area distribution (_sensu stricto_ plant area distribution) using the method proposed by Detto et al. [Journal of Geophysical Research, 2015], which inverts LiDAR distributions based on a stochastic radiative transfer model.  This includes modifications to the model by me to account for variable absorption of LiDAR pulses by canopy material.
- inventory_based_LAD_profiles.py : some code to produce crown volume distribution estimates based on detailed field observations and/or allometric relationships
- structural_metrics.py : some functions to extract quantitative metrics of canopy structure and structural complexity
- least_squares_fitting.py : routines to perform least squares fitting of various types, including polynomial and surface fitting, and affine transformations
- plot_LAD_profiles.py : makes some plots based on the GEM plot layout (probably specific to my purposes)
- canopy_structure_paper_figs.py : driver function to undertake all analysis in paper draft
- canopy_microclimate.py : some routines to estimate light environments within the canopy
- create_subplot_grids_least_squares_affine_transformation.py : a driver function to fit subplot grids to GEM plot vertices using a least squares affine (rotation & translation) transformation

There are plenty of other bits and bobs around, although many things are driver functions that are no longer used or are for more exploratory analyses, so I won't detail them here.

