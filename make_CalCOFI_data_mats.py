import csv
import scipy
import numpy as np
import seawater as sw
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

stat_ids = ['080.0 080.0', '093.3 030.0','093.3 110.0']
iso_pyc = [25.8,26.5,26.5]

# OPEN/READ IN THE CalCOFI DATA FILE
op_fil = open('/Users/elizabethdrenkard/Desktop/OMP/sample_bot_data.csv','rU') 
op_fil = open('/Users/elizabethdrenkard/Desktop/OMP/194903-201402_Bottle.csv','rU')
rd_fil = csv.reader(op_fil, delimiter=',')

# VARIABLES USED FOR WATER TYPE DEF
sw_vars = ['PTEMP','SALINITY','OXYGEN','NITRATE','PHOSPHATE']
nrow    = [50,53,56,60,59] 

# DEFINE DICTIONARY
CCS_sw_dict = {}
CCS_sw_iso_dict = {}
# EXTRACT DATA FOR SPECIFIC STATIONS AND MONTHS 
for row in rd_fil:
    if (row[2] == stat_ids[wt] and \
       #(70<=int(row[3][3:5])<80) and \
       (2<int(row[3][5:7])<6) and \
       (len(row[53])>0) and \
       (len(row[56])>0) and \
       (len(row[59])>0) and \
       (len(row[60])>0)):

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
# f, ax = plt.subplots(2, 2, sharey='row')

# PTEMP INCREASING MONOTONICLY
Ind = np.argsort(CCS_sw_iso_dict['PTEMP'])
PTEMP = np.array(CCS_sw_iso_dict['PTEMP'])[Ind]
npt = len(PTEMP)
print npt
print 'PTEMP MIN', np.min(PTEMP)
print 'PTEMP MAX', np.max(PTEMP)

f, ax = plt.subplots(2, 2, sharey='row')
f2, ax2 = plt.subplots(2, 2, sharey='row')
# ITERATE OVER ALL VARS
for nv in range(4):    
    VAR = np.array(CCS_sw_iso_dict[sw_vars[nv+1]])[Ind]

    # ERROR MINIMUM/REDUCTION FIGURES 
    # stor_err = np.zeros((npt-2,npt-2))
    # for strt_val in range(npt-2):
    #    for len_val in range(npt-(2+strt_val)): 
    #        c = strt_val + len_val + 3
    #        stor_err[strt_val,len_val] = np.sum((np.polyval(np.polyfit(PTEMP[strt_val:c], \
    #        VAR[strt_val:c], 1), PTEMP[strt_val:c]) - VAR[strt_val:c])**2)
    # plt.figure()
    # plt.pcolor(stor_err)
    # plt.figure()
    # plt.pcolor(np.diff(stor_err,axis=0))
    # plt.figure()
    # plt.pcolor(np.diff(stor_err,axis=1))
    # plt.colorbar()
    
    # POLYNOMIAL FITTING    
    # coeff, residuals, _, _, _   = np.polyfit(PTEMP, VAR,1,full=True)

    coeff = np.polyfit(PTEMP,VAR,1)
    # ALT ERROR CALC: np.sum((np.polyval(np.polyfit(PTEMP, VAR, 1), PTEMP) - VAR)**2)
    
    # VAR PLOT
    # plt.figure()
    # plt.plot(np.array(CCS_sw_dict[sw_vars[nv+1]]),np.array(CCS_sw_dict['PTEMP']),'o')
    # plt.plot(np.array(CCS_sw_iso_dict[sw_vars[nv+1]])[Ind],PTEMP,'bo',markeredgecolor='k')
     
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
title_txt = 'CalCOFI LineStation ' + stat_ids[wt]
plt.suptitle(title_txt)

f.subplots_adjust(hspace=0.4)

ax[0,0].set_ylabel('Temperature (C)')
ax[1,0].set_ylabel('Temperature (C)')

ax[0,0].set_xlabel('Salinity')
ax[0,1].set_xlabel('Oxygen (ml/L)')
ax[1,0].set_xlabel('Nitrate (umol/L)')
ax[1,1].set_xlabel('Phosphate (umol/L)')

f2.subplots_adjust(hspace=0.4)

ax2[0,0].set_ylabel('Temperature (C)')
ax2[1,0].set_ylabel('Temperature (C)')

ax2[0,0].set_xlabel('Salinity')
ax2[0,1].set_xlabel('Oxygen (ml/L)')
ax2[1,0].set_xlabel('Nitrate (umol/L)')
ax2[1,1].set_xlabel('Phosphate (umol/L)')

f.suptitle(title_txt)
plt.show()

