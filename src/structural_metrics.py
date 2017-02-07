## This library hosts functions to quantify aspects of the canopy structure, for example canopy heterogeneity in the horizontal and vertical dimensions.
import numpy as np

#-------------------------------------------------------------------------------
# Frechet number calculation - this algorithm was coded up by Max Bareiss and 
# can be found here:
#     https://www.snip2code.com/Snippet/76076/Fr-chet-Distance-in-Python
# The Frechet number is a metric to describe the similarity of two curves.  In
# this instance the curves are simplified to polygonal curves of an arbitrary
# number of points, giving the discrete Frechet difference.  This approximation # gives great advantages with respect to performance.
# The original algorithm was developed here:
#     Eiter, T. and Mannila, H., 1994. Computing discrete Fréchet distance.
#     Tech. Report CD-TR 94/64, Information Systems Department, Technical
#      University of Vienna.


# Calculate Euclidean distance between two points.
def euc_dist(pt1,pt2):
    return np.sqrt((pt2[0]-pt1[0])*(pt2[0]-pt1[0])+(pt2[1]-pt1[1])*(pt2[1]-pt1[1]))

def _c(ca,i,j,P,Q):
    if ca[i,j] > -1:
        return ca[i,j]
    elif i == 0 and j == 0:
        ca[i,j] = euc_dist(P[0],Q[0])
    elif i > 0 and j == 0:
        ca[i,j] = max(_c(ca,i-1,0,P,Q),euc_dist(P[i],Q[0]))
    elif i == 0 and j > 0:
        ca[i,j] = max(_c(ca,0,j-1,P,Q),euc_dist(P[0],Q[j]))
    elif i > 0 and j > 0:
        ca[i,j] = max(min(_c(ca,i-1,j,P,Q),_c(ca,i-1,j-1,P,Q),_c(ca,i,j-1,P,Q)),euc_dist(P[i],Q[j]))
    else:
        ca[i,j] = float("inf")
    return ca[i,j]

""" Computes the discrete frechet distance between two polygonal lines
Algorithm: http://www.kr.tuwien.ac.at/staff/eiter/et-archive/cdtr9464.pdf
P and Q are arrays of 2-element arrays (points)
"""
def frechetDist(P,Q):
    ca = np.ones((len(P),len(Q)))
    ca = np.multiply(ca,-1)
    return _c(ca,len(P)-1,len(Q)-1,P,Q)

#------------------------------------------------------------------------------
