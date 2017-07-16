## LiDAR_io.py
## This library contains the code to read LiDAR point clouds from .las files.
## It then creates kd-trees associated with each tile. This will subsequently
## allow v. efficient spatial searches to be done.  Moreover, the when dealing
## with >10^6 points (build time on our linux server ~0.5 s), we split the data
## into multiple trees due to the non-linear increase in build time with the
## number of points. 
