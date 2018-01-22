# ######### OMP analysis main program version 2  ###################

# Function calls used: qwt2.m qwt_tst.m nansum.m (Philip Morgan, CSIRO)
#  sw_ptmp sw_dens0.m (Philip Morgan, CSIRO) may be called for some data files
#  sw_dist.m (Philip Morgan, CSIRO) is called through the contour2 call

import scipy.io
import numpy as np
import netCDF4 as nc
import math
import pyroms
import matplotlib.pyplot as plt
from omp2 import omp2 
from collections import Counter
import datetime as dt
import matplotlib.dates as pltd

def qwt2(wm_row,ict):

    # WATER MASS ID VALUES
    wm = ('PSA', 'PSA', 'PEW', 'PEW', 'NPCW', 'NPCW')
    
    # WATER TYPE MATRIX
    # lower PSA
    # upper PSA
    # lower PEW
    # upper PEW
    # lower NPCW
    # upper NPCW

    #  The following lines define the water types. The order of parameters is
    #  Note: potential vorticity is multiplied by 10*8.
    wts=np.array(( \
    # PTEMP   SALT    OXY    PO4   NO3   SILICATE   mass   pvort
    ( 8.87,  33.59,  3.58,  1.72,  0.0,     23.20,   1.0,   0.0),\
    (11.91,  33.46,  5.44,  0.90,  0.0,      6.91,   1.0,   0.0),\
    ( 7.66,  34.32,  0.70,  2.86,  0.0,     55.10,   1.0,   0.0),\
    (10.26,  34.07,  2.39,  2.00,  0.0,     25.98,   1.0,   0.0),\
    ( 6.25,  34.13,  1.21,  2.81,  0.0,     63.11,   1.0,   0.0),\
    ( 9.37,  33.92,  3.66,  1.65,  0.0,     20.78,   1.0,   0.0)))

    G1=np.transpose(wts[wm_row,:])
    allsize = wts.shape

    return G1, wm, allsize


def sw_dist(lat,lon,units):

    # -----------------
    # DEFINE CONSTANTS
    # -----------------
    DEG2RAD = (2*math.pi/360)
    RAD2DEG = 1/DEG2RAD
    DEG2MIN = 60
    DEG2NM  = 60
    NM2KM   = 1.8520    # Defined in Pond & Pickard p303.

    npositions = len(lat)
    dlon = np.diff(lon)
    if any(abs(dlon)>180):
       flag = np.where(abs(dlon)>180)[0]
       for ii in flag:
           dlon[ii]= -np.sign(dlon[ii]) * (360 - abs(dlon[ii]))

    latrad = abs(lat*DEG2RAD)
    temp_vals = (latrad[1:]+latrad[:-1])/2.
    dep    = np.cos(temp_vals.squeeze()) * dlon
    dlat   = np.diff(lat.squeeze())

    dist   = DEG2NM*np.sqrt(dlat**2 + dep**2)  # in n.miles
    if units == 'km':   # defaults to n.miles
       dist = dist * NM2KM

    # CALCUALTE ANGLE TO X AXIS
    #phaseangle  = np.angle(dep+dlat*np.sqrt(-1))*RAD2DEG

    return dist #,phaseangle 



#############################################
# OMP2: Python edition
#############################################

print '  '
print 'OMP Analysis version 2 (March 1999)'
print '===================================  '
print '  '

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Originally from  incontr2:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
OMP = 'cla' # classical OMP analysis

# data location
# MACBOOK FILES:
dataset = '/Users/elizabethdrenkard/TOOLS/omp2/scripts/CalCOFI_LINE_093.3.npy'
dataset = '/Users/elizabethdrenkard/TOOLS/omp2/scripts/CalCOFI_LINE_080.0.npy'
# SWFSC FILES:
# dataset = '/Users/liz.drenkard/TOOLS/omp2/scripts/CalCOFI_LINE_093.3.npy'
# dataset = '/Users/liz.drenkard/TOOLS/omp2/scripts/CalCOFI_LINE_080.0.npy'

# data limitations
selection=  '(pdens>=23) & (pdens<=28)'# & (press>300) & (press<600)' 

# Select/deselect potential vorticity by setting switchpot to 'y' or 'n':
switchpot = 'n'

# WEIGHTS FOR VARIABLES
weightset='testwght.npy'

# number of water masses to be included in the analysis
wm = 3

#  Select the water type numbers (row in the water type matrix)
qwt_pos = [0,1,2,3,5] 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loading .npy files into python 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mat_dat=np.load(dataset).item()
globals().update(mat_dat)
weight_dat=np.load(weightset).item()
globals().update(weight_dat)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Variables available for/ to use in analysis
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# variable indices:
# 1: latitude 
# 2: longitude
# 3: pressure
# 4: salinity
# 5: potential temperature
# 6: oxygen
# 7: phosphate
# 8: nitrate
# 9: silicate
#10: potential vorticity
#11: temperature

# NOTE: For historical reasons the two columns mass conservation and potential vorticity are
# swapped in the program so that mass conservation is always the last column, after potential vorticity.
# The arrangement of the water type matrix and the weight vector thus differs from the description
# in the user manual. This should not be of concern but has to be watched when changing the code.

# ORIGINAL
esx = np.array((1,1,1,1,1,1,1,0,1,0,0)) 
# Read the weight and Redfield ratio file
# Check which weights are needed and reset the diagonal:
print Wx.shape
A    = np.diag(Wx)
A.setflags(write=1)
A1   = A[7]  # change order of weights so that mass conservation is last
A[7] = A[6]
A[6] = A1
ratio = ratio.squeeze()
# MANUALLY SET A
A[4] = 0
ratio[4] = -99999                    # no nitrate weight if no nitrate
A[6] = 0
ratio[6] = -99999                    # no pot. vorticity weight if not needed

statind = np.where(A>0)[0]
Wx = np.diag(A[statind])
statind = np.where(ratio>-99999)[0]
redfrat = ratio[statind]         # Redfield ratio for selected variables only
print '  '
# End of if statements for weights and Redfield ratio
G0,wmnames,i = qwt2(qwt_pos,0)

wm_index = []
wm_ind0  = []
wm_ind1  = []
j = 0
print '  '
tit_index = []
for i in range(len(qwt_pos)):
    wm_ind1 = wmnames[qwt_pos[i]]
    k = (wm_ind0==wm_ind1)
    if not k:
       j = j+1
       tit_index.extend([wmnames[qwt_pos[i]]])
    wm_ind0 = wm_ind1
    wm_index.extend([j])

nr_of_wm = wm_index[len(wm_index)-1]

# PTEMP, SALT, OXYGEN, PHOS, SI, MASS
i = (0,1,2,3,5,6)
G1 = G0[i,:]

# SELECT STATIONS 
stations = [26.7,28,30,35,40,45,50,55,60,70,80,90,100,110,120] #LINE 93.3
stations = [51,55,60,70,80,90,100]                             #LINE 80.0
nsta=len(stations)

surf_frac = np.array([], dtype=np.int64).reshape(0,wm,nsta)
cruise_dates = []
for yr in range(1980,2017+1):
    Iy = np.where(np.array(mat_dat['YEAR'])==yr)
    if len(Iy[0])>30:
       mons = np.array(list(set(np.array(mat_dat['MONTH'])[Iy[0]])))
       for mon in mons[(mons>2) & (mons<7)]:
           print yr, mon
           I = np.where((np.array(mat_dat['YEAR'])==yr) & (np.array(mat_dat['MONTH'])==mon))
           stats = np.array(mat_dat['STATION'])[I[0]]
           lat = np.array(mat_dat['LAT'])[I[0]]
           lon = np.array(mat_dat['LONG'])[I[0]]
           ptemp = np.array(mat_dat['PTEMP'])[I[0]] 
           sal = np.array(mat_dat['SALINITY'])[I[0]]
           pdens = np.array(mat_dat['PDENS'])[I[0]]
           oxy = np.array(mat_dat['OXYGEN'])[I[0]]
           ph = np.array(mat_dat['PHOSPHATE'])[I[0]]
           si = np.array(mat_dat['SILICATE'])[I[0]]
           press = np.array(mat_dat['PRESS'])[I[0]]
           dist = sw_dist(lat.squeeze(),lon.squeeze(),'km')
           # This is the main part of it all: The call to omp2.m which does the analysis
           surf_frac = np.concatenate((surf_frac,\
                       omp2(OMP,nr_of_wm,tit_index,qwt_pos,wmnames,Wx,lat,switchpot,selection,stations,stats,yr,mon,lon,esx,\
                       press,sal,oxy,ptemp,pdens,ph,si,G1,wm_index).reshape(1,wm,nsta)),axis=0) 
           cruise_dates = np.append(cruise_dates,pltd.date2num(dt.datetime(yr,mon,1)))

# SAVE variables as file for plotting figures
np.savez('water_mass_fractions_200m',surf_frac,cruise_dates)


