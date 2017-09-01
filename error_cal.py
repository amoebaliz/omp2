import numpy as np



# Vandermonde matrix of x
V = np.vander(X,2)

# QR decomposition: r = Triangular factor
Q,R = np.linalg.qr(V)

# Inverse of R
Rinv = np.linalg.inv(R)

# normr
# df


#  Estimate of the covariance matrix of p
Rinv * Rinv.transpose() * normr**2/df


