import numpy as np

def 1D_least_squares_polynomial(x,y,order):



def 2D_least_squares_polynomial(x,y,z,order):
    # first up - mask nodata values
    mask = np.isfinite(z)
    X = z[mask]
    Y = y[mask]
    Z = z[mask]
    ones = np.zeros(Z.size)
    A=np.array([])
    
    if order == 1:
        A = np.array([X,Y,ones])
    elif order == 2:
        # coefficients returned should be in order a-f for:
        # z = ax^2 + by^2 + cxy + dx + ey + f
        A = np.array([X**2, Y**2, X*Y, X, Y, ones])

    
    b = Z.copy()

    coeff, r, rank, s = np.linalg.lstsq(A, b)
