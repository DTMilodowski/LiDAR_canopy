## This library hosts functions to quantify aspects of the canopy structure, for example canopy heterogeneity 
## in the horizontal and vertical dimensions.
import numpy as np
from scipy import stats
import least_squares_fitting as lstsq
from matplotlib import pyplot as plt
#--------------------------------------------------------------------------------------------------------------
# Frechet number calculation - this algorithm was coded up by Max Bareiss and can be found here:
#     https://www.snip2code.com/Snippet/76076/Fr-chet-Distance-in-Python
# The Frechet number is a metric to describe the similarity of two curves.  In this instance the curves are
# simplified to polygonal curves of an arbitrary number of points, giving the discrete Frechet difference. 
# This approximation gives great advantages with respect to performance.
# The original algorithm was developed here:
#     Eiter, T. and Mannila, H., 1994. Computing discrete Frechet distance. Tech. Report CD-TR 94/64,
#     Information Systems Department, Technical University of Vienna.


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
def calculate_mean_Frechet_distance(vertical_profiles,heights):
    N_profiles, N_heights = vertical_profiles.shape
    N_pairs = N_profiles*(N_profiles-1)/2
    Frechet = np.zeros(N_pairs)
    index = 0
    line1=np.zeros((N_heights,2))
    line1[:,0] = heights.copy()
    line2=line1.copy()
    for i in range(0,N_profiles):
        line1[:,1] = vertical_profiles[i,:]
        for j in range(i+1,N_profiles):
            line2[:,1] = vertical_profiles[j,:]
            Frechet[index] = frechetDist(line1,line2)
            index+=1
    mean_Fr = np.mean(Frechet)
    return mean_Fr

#--------------------------------------------------------------------------------------------------------------
# calculate vertical forest structural index (VSI)
# This calculates a measure of forest structural heterogeneity based on the number and vertical distribution of 
# canopy layers picked out by the remote sensing product.  Follows method suggested in:
#    Tello, M., Cazcarra-Bes, V., Pardini, M. and Papathanassiou, K., 2015, July. Structural classification of
#    forest by means of L-band tomographic SAR. In Geoscience and Remote Sensing Symposium (IGARSS), 2015 IEEE 
#    International (pp. 5288-5291). IEEE.


def retrieve_peaks(vertical_profiles,heights):
    N_profiles,N_heights = vertical_profiles.shape
    peaks, peak_amplitude = find_maxima(heights,vertical_profiles[0])
    for i in range(1,N_profiles):
        peaks_this_iter, peak_amplitude = find_maxima(heights,vertical_profiles[i])
        peaks = np.concatenate((peaks,peaks_this_iter),axis=0)
    return peaks

# filter signal using savitzky-golay filter before retrieving peaks
def retrieve_peaks_with_savitzky_golay_filter(vertical_profiles,heights,filter_window,filter_order=3,threshold = 0):
    N_profiles,N_heights = vertical_profiles.shape
    profile = moving_polynomial_filter(vertical_profiles[0],filter_window,filter_order)
    
    peaks, peak_amplitude = find_maxima(heights,profile,threshold)
    plt.plot(vertical_profiles[0],heights,'-')
    plt.plot(profile,heights,'-')
    plt.plot(peak_amplitude,peaks,'o')
    plt.xlim(xmin=0)
    plt.show()
    for i in range(1,N_profiles):
        profile = moving_polynomial_filter(vertical_profiles[i],filter_window,filter_order)
        peaks_this_iter, peak_amplitude = find_maxima(heights,profile,threshold)
        peaks = np.concatenate((peaks,peaks_this_iter),axis=0)
        plt.plot(vertical_profiles[i],heights,'-')
        plt.plot(profile,heights,'-')
        plt.plot(peak_amplitude,peaks_this_iter,'o')
        plt.xlim(xmin=0)
        plt.show()
    return peaks


# filter signal using gaussian filter before retrieving peaks
def retrieve_peaks_gaussian_convolution(vertical_profiles,heights,sigma=1,threshold=0,plot_profiles=False):
    N_profiles,N_heights = vertical_profiles.shape
    profile = moving_gaussian_filter(vertical_profiles[0],sigma)
    
    peaks, peak_amplitude = find_maxima(heights,profile,threshold)
    index = np.zeros(peaks.size)
    if plot_profiles:
        plt.figure(1,facecolor='White',figsize=[4,4])
        plt.plot(vertical_profiles[0],heights,'-')
        plt.plot(profile,heights,'-')
        plt.plot(peak_amplitude,peaks,'o')
        plt.xlim(xmin=0)
        plt.show()
    for i in range(1,N_profiles):
        profile = moving_gaussian_filter(vertical_profiles[i],sigma)
        peaks_this_iter, peak_amplitude = find_maxima(heights,profile,threshold)
        peaks = np.concatenate((peaks,peaks_this_iter),axis=0)
        index = np.concatenate((index,np.ones(peaks_this_iter.size)*i),axis=0)
        if plot_profiles:
            plt.figure(1,facecolor='White',figsize=[4,4])
            plt.plot(vertical_profiles[i],heights,'-')
            plt.plot(profile,heights,'-')
            plt.plot(peak_amplitude,peaks_this_iter,'o')
            plt.xlim(xmin=0)
            plt.show()
    return peaks, index

def calculate_VSI(peaks):
    # convert number of peaks & their locations in canopy into VSI    
    N_peaks = peaks.size
    
    N_pairs = N_peaks*(N_peaks-1)/2
    separation = np.zeros(N_pairs)
    index = 0
    for i in range(0,N_peaks):
        for j in range(i+1,N_peaks):
            separation[index] = np.sqrt((peaks[i]-peaks[j])*(peaks[i]-peaks[j]))

    # CHECK THIS WITH MARIVI TELLO
    VSI = float(N_peaks)*np.mean(separation)
    return VSI

# alternative metric suggested by Marivi Tello is just to use vertical variance of peaks
def calculate_vertical_structural_variance(peaks):
    vertical_structural_variance = np.var(peaks)
    return vertical_structural_variance
    
# Likewise for the horizontal structural heterogeneity
def calculate_horizontal_structural_heterogeneity_alt(profiles,heights,stand_area):
    pk_hts = retrieve_peaks(vertical_profiles,heights)
    mean_profile = np.mean(vertical_profiles,axis=1)
    stand_ht = np.max(heights[mean_profile>0])
    Npks_upper_canopy = float(np.sum(pk_hts>=0.6*stand_ht))
    HSI = N_peaks_upper_canopy*stand_ht/stand_area
    return HSI



# Simple function to find local maxima based on immediate neighbourhood
def find_maxima(signal_x, signal_y, threshold=0):
    pks = []
    N=signal_x.size
    for i in range(1,N-1):
        if np.all((signal_y[i]>=threshold,signal_y[i]>signal_y[i-1],signal_y[i]>signal_y[i+1])):
            pks.append(i)
    Npks = len(pks)
    peak_x = np.zeros(Npks)
    peak_y=peak_x.copy()
    for i in range(0,Npks):
        peak_x[i] = signal_x[pks[i]]
        peak_y[i] = signal_y[pks[i]]
    return peak_x, peak_y



# function to smooth using moving window with polynomial fit (default is second order).
# Specify window half width (in pixels) for fitting polynomial
# i.e. a Savitzky-Golay filter.  Currently supports up to 4th order fitting
def moving_polynomial_filter(signal_y,window_width,order=2):
    N = signal_y.size
    window_half_width = window_width//2
    y_filt = np.zeros(N)

    firstvals = signal_y[0] - np.abs(signal_y[1:window_half_width+1][::-1] - signal_y[0] )
    lastvals = signal_y[-1] + np.abs(signal_y[-window_half_width-1:-1][::-1] - signal_y[-1])
    y_temp = np.concatenate((firstvals, signal_y, lastvals))

    for i in range(0,N):
        ii = i + window_half_width
        x = np.arange(-window_half_width,window_half_width+1)
        y = y_temp[ii-window_half_width:ii+window_half_width+1]
        coeffs = lstsq.oneD_least_squares_polynomial(x,y,order)
        # at x=0 - the centre of the moving window - polynomial reduces down to constant
        y_filt[i] = coeffs[-1]
    
    return y_filt


# function to convolve a signal with a gaussian curve (comprising three standard deviations)
# to provide a signal filter.
def get_gaussian_kernel(sigma):
    x = np.arange(2*3*sigma+1)-3*sigma
    kernel = np.exp(-(x**2)/(2*sigma**2))
    return kernel/kernel.sum()

# Specify i) input signal, ii) standard deviation width, sigma.  sigma must be an integer, and represents a number of cells
# the total width of the gaussian filter will be 2*3*sigma+1 cells
# Boundary conditions are reflected
def moving_gaussian_filter(signal_y,sigma):
    N = signal_y.size
    window_width = 3*2*sigma+1
    window_half_width = window_width//2
    y_filt = np.zeros(N)
    kernel = get_gaussian_kernel(float(sigma))
    firstvals = signal_y[0] - np.abs(signal_y[1:window_half_width+1][::-1] - signal_y[0] )
    lastvals = signal_y[-1] + np.abs(signal_y[-window_half_width-1:-1][::-1] - signal_y[-1])
    y_temp = np.concatenate((firstvals, signal_y, lastvals))
    y_filt = np.convolve(y_temp,kernel,'valid')
    return y_filt




# Calculate Shannon Index
# for use of this metric to assess canopy structural diversity see MacArthur & MacArthur, 1961; Stark et al., 2012
def calculate_Shannon_index(P):
    p=P/P.sum()
    S = -np.sum(p[p>0]*np.log(p[p>0]))
    return S


# Other useful structural metrics
# Canopy Shape [Asner et al., Biogeosciences, 2014]
# P is the height at which there is a maximum "canopy volume"
# H is the 99th percentile of the height distribution
# Asner and following papers use the LiDAR point cloud, but this doesn't really tell you about actual vegetation distribution
# Could also consider using the derived plant area profiles, although even then, unclear whether more useful than distribution
# moments
def calculate_canopy_shape(heights,density):
    # make sure profiles are orientated correctly
    if heights[0]>heights[1]:
        heights = heights[::-1]
        density = density[::-1]
    
    cumulative_density = np.cumsum(density)
    
    perc99 = 0.99*cumulative_density[-1]
    H = heights[cumulative_density>=perc99][0]
    P = heights[density==density.max()][0]
    return H, P, P/H

# Calculate moments
# Descriptive statistics to characterise distributions
# Mean is obvious
# Std deviation tells you about spread around the mean
# Skew tells you about distribution within the profile (e.g. heavy tailed etc.)
# Kurtosis tells you about "peakiness of distribution"
def calculate_moments_of_distribution(x,p):
    N = p.sum()
    u = np.sum(x*p)/N # mean
    v = np.sum(p*(x-u)**2)/N  # variance
    s = (np.sum(p*(x-u)**3)/N) / ((np.sum(p*(x-u)**2)/N)**(3/2.))  # skew
    k = (np.sum(p*(x-u)**4)/N) / ((np.sum(p*(x-u)**2)/N)**(2.)) # kurtosis

    return u,np.sqrt(v),s,k

# Calculate number of canopy layers based on regions of continuous PAD [Clark et al., Ecology Letters, 2008]
# This is trying to replicate the field process of counting canopy layers
# simpler and easier to interpret than picking peaks from smoothed distributions
def calculate_number_of_contiguous_layers(heights,density,minimum_density):
    is_layer = np.zeros(heights.size+1)
    is_layer[:-1] = density>=minimum_density
    n_layers = np.sum(is_layer[:-1]-is_layer[1:]==1)
    return n_layers
