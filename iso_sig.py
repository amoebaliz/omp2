import csv
import scipy
import numpy as np
import seawater as sw
import matplotlib.pyplot as plt
import datetime as dt
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

CCS_sw_iso_dict = {}

# EXTRACT DATA FOR SPECIFIC STATIONS AND MONTHS 
for row in rd_fil:
    if (row[2] == stat_ids[wt] and \
       (int(row[3][:2] + row[3][3:5])>=1980) and \
       (len(row[53])>0) and \
       (len(row[56])>0) and \
       (len(row[59])>0) and \
       (len(row[60])>0)):


       # ISOPYCNAL RESTRICTION
       if ((np.float(row[53])>iso_pyc[wt]-mkt) and \
           (np.float(row[53])<iso_pyc[wt]+pkt)):

          yr_val = int(row[3][:2] + row[3][3:5])
          day_val = int(row[3][13:16])
          date_val = dt.datetime(yr_val,1,1)+dt.timedelta(days=day_val-1)
          insertIntoDataStruct('TIME', date_val,CCS_sw_iso_dict)

          # PUT SPECIFIC STATION/SEASONAL DATA IN DICTIONARY
          for nvar in range(len(sw_vars)):
              insertIntoDataStruct(sw_vars[nvar], np.float64(row[nrow[nvar]]),CCS_sw_iso_dict)
         
           #plt.plot(date_val,np.float64(row[56]),'bo') 
plt.figure()

date_val = np.array(CCS_sw_iso_dict['TIME'])
var_vals = np.array(CCS_sw_iso_dict[sw_vars[3]])
plt.plot(date_val,var_vals,'bo') 
#coeff = np.polyfit(date_val,var_vals,1)
#plt.plot(date_val,np.polyval(coeff,date_val),'r')


title_txt = 'Phosphate Concentration at ' + stat_ids[wt] + ', isopycnal ' + str(iso_pyc[wt]) 
plt.title(title_txt)
plt.ylabel('Phosphate (umol/L)')
plt.show()


