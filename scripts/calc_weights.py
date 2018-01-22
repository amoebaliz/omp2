import numpy as np

G = np.array(( \
    # PTEMP   SALT    OXY    PO4   NO3   SILICATE   mass   pvort
    ( 8.87,  33.59,  3.58,  1.72,  0.0,     23.20,   1.0,   0.0),\
    (11.91,  33.46,  5.44,  0.90,  0.0,      6.91,   1.0,   0.0),\
    ( 7.66,  34.32,  0.70,  2.86,  0.0,     55.10,   1.0,   0.0),\
    (10.26,  34.07,  2.39,  2.00,  0.0,     25.98,   1.0,   0.0),\
    ( 6.25,  34.13,  1.21,  2.81,  0.0,     63.11,   1.0,   0.0),\
    ( 9.37,  33.92,  3.66,  1.65,  0.0,     20.78,   1.0,   0.0)))

nwt = 3

# PTEMP   SALT    OXY    PO4   NO3   SILICATE   mass   pvort
if nwt == 3:
       # WITH NPCW
       #                      NPCW   PSA    NPCW   PSA   (NA)  NPCW    (NA)  (NA) 
       max_var =  np.array(( 0.629, 0.017, 0.255, 0.031,  1,  11.721,  1,    1)) 

elif nwt == 2:
       # WITHOUT NPCW
       #                     PEW   PSA   PSA   PSA  (NA) PEW    (NA) (NA) 
       max_var =  np.array((0.531,0.017,0.195,0.031,  1, 6.916, 1,  1))

var_wts = np.var(G,axis=0)
wts = var_wts/max_var
wts[-2:] = np.max(wts)
print wts 
