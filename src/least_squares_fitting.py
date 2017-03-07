import numpy as np

def oneD_least_squares_polynomial(x,y,order):
    return 0


def twoD_least_squares_polynomial(x,y,z,order=1):
    # first up - mask nodata values
    mask = np.isfinite(z)
    X = x[mask]
    Y = y[mask]
    Z = z[mask]
    ones = np.ones(Z.size)
    A=np.array([])
    
    if order == 1:
        # coefficients returned should be in order a-f for:
        # z = ax + by + c
        A = np.array([X,Y,ones]).T
    elif order == 2:
        print "order = 2"
        # coefficients returned should be in order a-f for:
        # z = ax^2 + by^2 + cxy + dx + ey + f
        A = np.array([X**2, Y**2, X*Y, X, Y, ones]).T
    b = Z.copy()
    coeff, r, rank, s = np.linalg.lstsq(A, b)
    return coeff
