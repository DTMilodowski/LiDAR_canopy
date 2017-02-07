## This library hosts functions to quantify aspects of the canopy structure, for example canopy heterogeneity in the horizontal and vertical dimensions.
import numpy as np

#--------------------------------------------------------------------------------------------------------------
# Frechet number calculation - this algorithm was coded up by Max Bareiss and can be found here:
#     https://www.snip2code.com/Snippet/76076/Fr-chet-Distance-in-Python
# The Frechet number is a metric to describe the similarity of two curves.  In this instance the curves are
# simplified to polygonal curves of an arbitrary number of points, giving the discrete Frechet difference. 
# This approximation gives great advantages with respect to performance.
# The original algorithm was developed here:
#     Eiter, T. and Mannila, H., 1994. Computing discrete Fréchet distance.
#     Tech. Report CD-TR 94/64, Information Systems Department, Technical
#     University of Vienna.

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

#--------------------------------------------------------------------------------------------------------------

# calculate mean Frechet distance -> horizontal structural index suggested in:
#    Tello, M., Cazcarra-Bes, V., Pardini, M. and Papathanassiou, K., 2015, July. Structural classification of
#    forest by means of L-band tomographic SAR. In Geoscience and Remote Sensing Symposium (IGARSS), 2015 IEEE 
#    International (pp. 5288-5291). IEEE.
def calculate_mean_Frechet_distance(vertical_profiles):

    return mean_Fr

#--------------------------------------------------------------------------------------------------------------
# calculate vertical forest structural index (VSI)
# This calculates a measure of forest structural heterogeneity based on the number and vertical distribution of 
# canopy layers picked out by the remote sensing product.  Follows method suggested in:
#    Tello, M., Cazcarra-Bes, V., Pardini, M. and Papathanassiou, K., 2015, July. Structural classification of
#    forest by means of L-band tomographic SAR. In Geoscience and Remote Sensing Symposium (IGARSS), 2015 IEEE 
#    International (pp. 5288-5291). IEEE.
def calculate_VSI(vertical_profiles):

    return VSI


def find_maxima(signal_x, signal_y):

    pks = []
    N=signal_x.size
    for i in range(1,N-1):
        if np.all((signal_y[i]>signal_y[i-1],signal_y[i]>signal_y[i+1])):
            pks.append(i)
    Npks = len(pks)
    peak_x = np.zeros(Npks)
    peak_y=peak_x.copy()
    for i in range(0,Npks):
        peak_x = signal_x[pks[i]]
        peak_y = signal_y[pks[i]]
    return peak_x, peak_y
