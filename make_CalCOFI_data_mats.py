import csv
import scipy
import numpy as np
#import seawater as sw
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

# CalCOFI data codes
# http://calcofi.org/new.data/index.php/reporteddata/2013-09-30-23-23-27/database-tables-description

# ----------------
#   ROW  |  DATA
# ----------------
#   02      Station ID (Sta_ID)             CalCOFI Line Station
#   03      Sample ID  (Depth_ID)  
#   49      Depth      (R_DEPTH)            m
#   50      Temp       (R_TEMP)             oC
#   51      Pot. Temp  (R_POTEMP)           oC
#   52      Salinity   (R_SALINITY)         umol/L
#   53      Pot. Dens  (R_SIGMA)            n.a.
#   56      Oxygen     (R_O2)               ml/L
#   58      Silicate   (R_SIO3)             umol/L
#   59      Phosphate  (R_PO4)              umol/L 
#   60      Nitrate    (R_NO3)              umol/L 
#   65      Pressure   (R_PRES)             db
#   67      Oxygen     (R_Oxy_mumol_per_Kg) umol/kg


def insertIntoDataStruct(name,val,aDict):
    if not name in aDict:
        aDict[name] = [(val)]
    else:
        aDict[name].append((val))
# SELECT STATION FOR SPECIFIC WATER TYPE
wt = 2
pkt = .2  # Buffer for potential density (isopycnal) values 
mkt = .2

# WATER MASS SAMPLES: 
stat_ids = ['080.0 080.0', '093.3 030.0','093.3 110.0']
iso_pyc = [25.8,26.5,26.5]

# OPEN/READ IN THE CalCOFI DATA FILE
#op_fil = open('/Users/liz.drenkard/external_data/CalCOFI/194903-201402_Bottle.csv','rU')
op_fil = open('/Users/elizabethdrenkard/Desktop/OMP/194903-201402_Bottle.csv','rU')
rd_fil = csv.reader(op_fil, delimiter=',')

# VARIABLES USED FOR WATER TYPE DEF
sw_vars = ['PTEMP','SALINITY','OXYGEN','SILICATE','PHOSPHATE']
nrow    = [51,52,56,58,59] 
# DEFINE DICTIONARY
CCS_sw_dict = {}
CCS_sw_iso_dict = {}
# EXTRACT DATA FOR SPECIFIC STATIONS AND MONTHS 
for row in rd_fil:
    if (row[2] == stat_ids[wt] and \
       #(70<=int(row[3][3:5])<80) and \
       (2<int(row[3][5:7])<6) and \
#       (5<int(row[3][5:7])<10) and \
       (len(row[51])>0) and \
       (len(row[52])>0) and \
       (len(row[53])>0) and \
       (len(row[56])>0) and \
       (len(row[59])>0) and \
       (len(row[58])>0)):

       # ALL DATA DICTIONARY
       for nvar in range(len(sw_vars)):
           insertIntoDataStruct(sw_vars[nvar], np.float64(row[nrow[nvar]]),CCS_sw_dict)

       # ISOPYCNAL RESTRICTION
       if ((np.float(row[53])>iso_pyc[wt]-mkt) and \
           (np.float(row[53])<iso_pyc[wt]+pkt)):

           # PUT SPECIFIC STATION/SEASONAL DATA IN DICTIONARY
           for nvar in range(len(sw_vars)): 
               insertIntoDataStruct(sw_vars[nvar], np.float64(row[nrow[nvar]]),CCS_sw_iso_dict)

# ANALYSIS AND PLOTTING

# PTEMP INCREASING MONOTONICLY
# Ind = np.argsort(CCS_sw_dict['PTEMP'])
# PTEMP = np.array(CCS_sw_dict['PTEMP'])[Ind]

Ind = np.argsort(CCS_sw_iso_dict['PTEMP'])
PTEMP = np.array(CCS_sw_iso_dict['PTEMP'])[Ind]

npt = len(PTEMP)

print 'PTEMP MIN', np.min(PTEMP)
print 'PTEMP MAX', np.max(PTEMP)

f, ax = plt.subplots(2, 2, sharey='row')
f2, ax2 = plt.subplots(2, 2, sharey='row')
# ITERATE OVER ALL VARS

n = 1 # FOR LINEAR FIT

# DEGREES OF FREEDOM
df = len(PTEMP)-(n+1)
# Vandermonde matrix of x
V = np.vander(PTEMP,n+1)
# QR decomposition: 
Q,R = np.linalg.qr(V)
# ERROR for x
err_fit = np.sqrt(1+np.sum(Q**2,axis=1))

err = np.ones(4)
for nv in range(4):    
    VAR = np.array(CCS_sw_iso_dict[sw_vars[nv+1]])[Ind]
    # VAR = np.array(CCS_sw_dict[sw_vars[nv+1]])[Ind]

    # FITTING POLYNOMIAL
    coeff = np.polyfit(PTEMP,VAR,n)

    # RESIDUALS CALCULATION
    r = VAR - np.dot(V,coeff)
    delta   = np.linalg.norm(r)/df*err_fit
    err[nv] = np.mean(delta)
     
    # SUBPLOTS
    a = nv/2; b = -2*(nv/2) + nv
    ax[a,b].plot(VAR,PTEMP,'b+')
    ax[a,b].plot(np.polyval(coeff,PTEMP),PTEMP,'r')
    ax2[a,b].plot(np.array(CCS_sw_dict[sw_vars[nv+1]]),np.array(CCS_sw_dict['PTEMP']),'o')
    ax2[a,b].plot(np.array(CCS_sw_iso_dict[sw_vars[nv+1]])[Ind],PTEMP,'bo',markeredgecolor='k')


    print '-------------------'
    print sw_vars[nv+1], 'ERROR = ', np.sum((np.polyval(np.polyfit(PTEMP, VAR, 1), PTEMP) - VAR)**2)
    print sw_vars[nv+1], 'MIN = ', np.min(np.polyval(coeff,PTEMP)) 
    print sw_vars[nv+1], 'MAX = ', np.max(np.polyval(coeff,PTEMP)) 

#plt.show()
# FIGURE DETAILS
print err
title_txt = 'CalCOFI LineStation ' + stat_ids[wt]
plt.suptitle(title_txt)

f.subplots_adjust(hspace=0.4)

ax[0,0].set_ylabel('Temperature (C)')
ax[1,0].set_ylabel('Temperature (C)')

ax[0,0].set_xlabel('Salinity')
ax[0,1].set_xlabel('Oxygen (ml/L)')
ax[1,0].set_xlabel('Silicate (umol/L)')
ax[1,1].set_xlabel('Phosphate (umol/L)')

f2.subplots_adjust(hspace=0.4)

ax2[0,0].set_ylabel('Temperature (C)')
ax2[1,0].set_ylabel('Temperature (C)')

ax2[0,0].set_xlabel('Salinity')
ax2[0,1].set_xlabel('Oxygen (ml/L)')
ax2[1,0].set_xlabel('Silicate (umol/L)')
ax2[1,1].set_xlabel('Phosphate (umol/L)')

f.suptitle(title_txt)
plt.show()

