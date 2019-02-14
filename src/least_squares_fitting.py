import numpy as np
# This function calculates best fitting polynomial line (up to order 3) for a
# given set of input points x and y, using linear least squares inversion.
# Default is for 1st order.
# x and y are one dimensional arrays
def oneD_least_squares_polynomial(x,y,order=1):
    # first up - mask nodata values
    mask = np.isfinite(y)
    X = x[mask]
    Y = y[mask]
    ones = np.ones(Y.size)
    A=np.array([])

    if order == 1:
        # coefficients returned should be in order a-b for:
        # z = ax + b
        A = np.array([X,ones]).T
    elif order == 2:
        # coefficients returned should be in order a-c for:
        # z = ax^2 + bx + c
        A = np.array([X**2, X, ones]).T
    elif order == 3:
        # coefficients returned should be in order a-d for:
        # z = ax^3 + bx^2 + cx + d
        A = np.array([X**3, X**2, X, ones]).T
    elif order == 4:
        # coefficients returned should be in order a-e for:
        # z = ax^4 + bx^3 + cx^2 + dx + e
        A = np.array([X**4, X**3, X**2, X, ones]).T
    b = Y.copy()
    coeff, r, rank, s = np.linalg.lstsq(A, b)
    return coeff

# This function calculates best fitting polynomial surface (up to order 3) for a
# given set of input points x,y and z, using linear least squares inversion.
# Default is for 1st order.
# x,y and z are one dimensional arrays
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
        # coefficients returned should be in order a-f for:
        # z = ax^2 + by^2 + cxy + dx + ey + f
        A = np.array([X**2, Y**2, X*Y, X, Y, ones]).T
    elif order == 3:
        # coefficients returned should be in order a-j for:
        # z = ax^3 + by^3 + cyx^2 + dxy^2 + ex^2 + fy^2 + gxy + hx + iy + j
        A = np.array([X**3, Y**3, X**2*Y, Y**2*X, X**2, Y**2, X*Y, X, Y, ones]).T
    b = Z.copy()
    coeff, r, rank, s = np.linalg.lstsq(A, b)
    return coeff

# perform a least squares affine transformation to map a set of coordinates x,y to a transformed coordinate set x',y'.
# input arguments:
#   x -       the original x coordinate
#   y -       the original y coordinate
#   x_prime - the target x coordinate
#   y_prime - the target y coordinate
#   shear   - a flag to determine whether points can be sheared in addition
#             to rotation and translation (default is False)
def least_squares_affine_matrix(x,y,x_prime,y_prime,shear=False):
    # first get points for a given plot and build matrices
    n_points = x.size
    A=np.zeros((2*n_points,6))
    b=np.zeros(2*n_points)

    for ii in range(0,n_points):
        A[ii,0]=x[ii]
        A[ii,1]=y[ii]
        A[ii,2]=1

        A[ii+n_points,3]=x[ii]
        A[ii+n_points,4]=y[ii]
        A[ii+n_points,5]=1

        b[ii] = x_prime[ii]
        b[ii+n_points] = y_prime[ii]

    # now perform least squares inversion to find parameters of the optimal 2D rotation-translation affine transformation
    h, res, rank, s = np.linalg.lstsq(A,b)
    # now construct the affine transformation matrix
    affine = np.zeros((3,3))
    affine[0,0]=h[0]
    affine[0,1]=h[1]
    affine[0,2]=h[2]
    affine[1,0]=h[3]
    affine[1,1]=h[4]
    affine[1,2]=h[5]
    affine[2,0]=0
    affine[2,1]=0
    affine[2,2]=1

    return affine


# Apply affine transformation
def apply_affine_transformation(x,y,affine):

    X = np.asarray([x,y,np.ones(x.size)])
    X_prime = np.dot(affine,X)

    return X_prime[0], X_prime[1]
