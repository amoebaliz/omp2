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
    (  9.1,  33.50,  3.91,  1.56,  0.0,     20.30,   1.0,   0.0),\
    (11.34,  33.58,  4.93,  1.14,  0.0,     11.25,   1.0,   0.0),\
    ( 7.37,  34.04,  2.45,  2.32,  0.0,     43.42,   1.0,   0.0),\
    ( 9.75,  33.85,  3.66,  1.68,  0.0,     23.65,   1.0,   0.0),\
    ( 6.94,  34.06,  2.15,  2.43,  0.0,     49.97,   1.0,   0.0),\
    ( 9.14,  33.97,  3.24,  1.80,  0.0,     26.16,   1.0,   0.0))) 

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
#dataset = '/Users/liz.drenkard/TOOLS/omp2/scripts/CalCOFI_LINE_093.3.npy'
dataset = '/Users/elizabethdrenkard/TOOLS/omp2/scripts/CalCOFI_LINE_093.3.npy'
# data limitations
selection=  '(pdens>=23) & (pdens<=28)'# & (press>300) & (press<600)' 
# Select/deselect potential vorticity by setting switchpot to 'y' or 'n':
switchpot = 'n'
# Select/deselect variables by setting corresponding switches to 'y' or 'n':
iox = 1 #'y' # oxygen switch
iph = 1 #'y' # phosphate switch
ini = 0 #'n' # nitrate switch
isi = 1 #'y' # silicate switch
var_switches = [1,1,0,1]
# file which contains the weights 
# weightset='/Users/elizabethdrenkard/TOOLS/omp2mats/testwght.mat'
weightset='testwght.npy'
# number of water masses to be included in the analysis
wm = 3
#  Select the water type numbers (row in the water type matrix)
qwt_pos = [0,2,4] # changed from [1,2,3,4]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loading .mat files into python - soon to be netCDF
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

esx = np.array((1,1,1,1,1,1,1,0,1,0,0)) 

# Read the weight and Redfield ratio file
# Check which weights are needed and reset the diagonal:
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

for yr in range(1995,2014+1):
    Iy = np.where(np.array(mat_dat['YEAR'])==yr)
    if len(Iy[0])>30:
       mons = np.array(list(set(np.array(mat_dat['MONTH'])[Iy[0]])))
       for mon in mons[(mons>2) & (mons<7)]:
           print yr, mon
           I = np.where((np.array(mat_dat['YEAR'])==yr) & (np.array(mat_dat['MONTH'])==mon))
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
           #dist,phaseangle = sw_dist(lat.squeeze(),lon.squeeze(),'km')
           cumdist=np.append(0, np.cumsum(dist))
           # This is the main part of it all: The call to omp2.m which does the analysis
           omp2(OMP,nr_of_wm,tit_index,qwt_pos,wmnames,Wx,lat,switchpot,selection,lon,esx,press,sal,oxy,ptemp,pdens,ph,si,G1,wm_index)


#ORIGINAL esx calculation method
#eex = np.zeros(11) # index of available variables:
#key_vars = ['LAT','LONG','press','SALINITY','PTEMP','OXYGEN','PHOSPHATE','ni','SILICATE','pvort','temp']
#for n in range(len(key_vars)):
#    if key_vars[n] in globals():
#       eex[n]=1
    #elif n == len(key_vars)-1:
    #   temp = sw_temp(sal,ptemp,press,0)
    #   eex[n]=1

# Determine the number of variables used in this run:
#esx=np.copy(eex) # index of selected variables
#nvar = 3 + np.sum(var_switches)

#for nt in range(len(var_switches)):
#    esx[nt+5] = var_switches[nt]
# Read the weight and Redfield ratio file
# Check which weights are needed and reset the diagonal:
#A    = np.diag(Wx)

#A.setflags(write=1)
#A1   = A[7]  # change order of weights so that mass conservation is last
#A[7] = A[6]
#A[6] = A1
#ratio = ratio.squeeze()

#if esx[4] == 0:
#   A[0] = 0
#   ratio[0] = -99999                    # no pot. temperature weight if not needed
#if esx[3] == 0:
#   A[1] = 0
#   ratio[1] = -99999                    # no salinity weight if not needed
#if esx[5] == 0:
#   A[2] = 0
#   ratio[2] = -99999                    # no oxygen weight if no oxygen
#if esx[6] == 0:
#   A[3] = 0
#   ratio[3] = -99999                    # no phosphate weight if no phosphate
#if esx[7] == 0:
#   A[4] = 0
#   ratio[4] = -99999                    # no nitrate weight if no nitrate
#if esx[8] == 0:
#   A[5] = 0
#   ratio[5] = -99999                    # no silicate weight if no silicate
#if esx[9] == 0:
#   A[6] = 0
#   ratio[6] = -99999                    # no pot. vorticity weight if not needed

#statind = np.where(A>0)[0]
#Wx = np.diag(A[statind])
#statind = np.where(ratio>-99999)[0]
#redfrat = ratio[statind]         # Redfield ratio for selected variables only
#print '  '
# End of if statements for weights and Redfield ratio

# Read the water types
#G0,wmnames,i = qwt2(qwt_pos,0)

#wm_index = []
#wm_ind0  = []
#wm_ind1  = []
#j = 0
#print '  '
#tit_index = []
#for i in range(len(qwt_pos)):
#    wm_ind1 = wmnames[qwt_pos[i]]
#    k = (wm_ind0==wm_ind1)
#    if not k:
#       j = j+1
#       tit_index.extend([wmnames[qwt_pos[i]]])
#    wm_ind0 = wm_ind1
#    wm_index.extend([j])

#nr_of_wm = wm_index[len(wm_index)-1]

#i = 2
#del G1

# G1 add rows
#G1 = G0[:2,:] # PTEMP, SALT

#if esx[5] == 1: # OXYGEN
#   G1 = np.concatenate((G1,np.array([G0[2,:]])),axis=0)
#   i = i+1
#if esx[6] == 1: # PHOSPHATE
#   G1 = np.concatenate((G1,np.array([G0[3,:]])),axis=0)   
#   i = i+1
#if esx[7] ==1:  # NITRATE
#   G1 = np.concatenate((G1,np.array([G0[4,:]])),axis=0)
#   i = i+1
#if esx[8] == 1: # SILICATE
#   G1 = np.concatenate((G1,np.array([G0[5,:]])),axis=0)
#   i = i+1
#if esx[9] == 1: # MASS
#   G1 = np.concatenate((G1,np.array([G0[7,:]])),axis=0) # Concatenate pvort  
#   i = i+1
#G1 =  np.concatenate((G1,np.array([G0[6,:]])),axis=0) # Concatenate mass
# This is the main part of it all: The call to omp2.m which does the analysis
#omp2(OMP,nr_of_wm,tit_index,qwt_pos,wmnames,Wx,lat,switchpot,selection,lon,esx,press,sal,oxy,ptemp,temp,pdens,ph,ni,G1,wm_index)
#omp2(OMP,nr_of_wm,tit_index,qwt_pos,wmnames,Wx,lat,switchpot,selection,lon,esx,press,sal,oxy,ptemp,pdens,ph,si,G1,wm_index)
# It's all done. Documentation and display is all in omp2.m.
