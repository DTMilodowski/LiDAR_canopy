# This is a quick script to convert laz file to las files using las2las, which
# is one of the open source components of lastools.
# The basic idea is to loop through a list of laz files and convert to las files
# using the identical naming structure.  It's a bit of a hash, but it works ok.
import numpy as np
import os

# Declare paths
laz_dir = '/disk/scratch/local.2/dmilodow/BALI/LiDAR/Danum/Danum_CHM_laz_tiles'
las_dir = '/disk/scratch/local.2/dmilodow/BALI/LiDAR/Danum'

# create file list
os.system("ls %s/*.laz > %s/laz_list.txt" % (laz_dir,laz_dir))
laz_files = np.genfromtxt("%s/laz_list.txt" % laz_dir,delimiter=',',dtype='S256')
N = len(laz_files)
for ii in range(0,N):
    filename = laz_files[ii].split('/')[-1]
    prefix = filename.split('.')[0]
    print "tile %i of %i: %s" % (ii+1,N,prefix)
    las_file = '%s/%s.las' % (las_dir,prefix)
    os.system("las2las %s %s" % (laz_file, las_file))
    
print "Done" 
