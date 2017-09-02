import numpy as np

x = [1,2,3] 
y = [2,3.9,6]
n = 1

# POLYFIT
p = np.polyfit(x,y,n)
yp = p[0]*np.array(x)+p[1]

# RESIDUALS
r = y-yp

# DEGREES OF FREEDOM
df = len(x)-(n+1) 

# Vandermonde matrix of x
V = np.vander(x,2)

# QR decomposition: r = Triangular factor
Q,R = np.linalg.qr(V)

if n == 0:
   e = np.sqrt(1+Q**2)
else:
   e = np.sqrt(1+np.sum(Q**2,axis=1))

delta = np.linalg.norm(r)/df*e

print delta



