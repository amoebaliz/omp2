import csv
import scipy
import pickle
import numpy as np
import matplotlib.pyplot as plt

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

# SELECT LINE 
stat_ids = ['080.0 080.0', '093.3 030.0','093.3 110.0']
line_num = 93.3

# OPEN/READ IN THE CalCOFI DATA FILE
# SWFSC
dat_fil = open('/Users/liz.drenkard/external_data/CalCOFI/194903-201402_Bottle.csv','rU')
#pos_fil = open('/Users/liz.drenkard/external_data/CalCOFI/CalCOFIStaPosNDepth113.csv')
# MACBOOK
#dat_fil = open('/Users/elizabethdrenkard/external_data/CalCOFI/194903-201402_Bottle.csv','rU')
#pos_fil = open('/Users/elizabethdrenkard/external_data/CalCOFI/CalCOFIStaPosNDepth113.csv','rU')

rd_dat_fil = csv.reader(dat_fil, delimiter=',')
#rd_pos_fil = csv.reader(pos_fil, delimiter=',')

# VARIABLES USED FOR WATER TYPE DEF
sw_vars = ['DEPTH','PTEMP','SALINITY','PDENS','OXYGEN','SILICATE','PHOSPHATE','PRESS']
nrow    = [49,51,52,53,56,58,59,65] 

# DEFINE DICTIONARY
CCS_sw_dict = {}
s=''
# EXTRACT DATA FOR SPECIFIC STATIONS AND MONTHS 
sta=0
for row in rd_dat_fil:
    if (row[2][:5] == str(line_num).zfill(5) and \
       (len(row[51])>0) and \
       (len(row[52])>0) and \
       (len(row[53])>0) and \
       (len(row[56])>0) and \
       (len(row[58])>0) and \
       (len(row[59])>0)):
       # ADD LAT/LON BASED ON LINE AND STATION
       pos_fil = open('/Users/liz.drenkard/external_data/CalCOFI/CalCOFIStaPosNDepth113.csv','rU')
       rd_pos_fil = csv.reader(pos_fil, delimiter=',')
       for p_row in rd_pos_fil:
           if ((p_row[1] == str(line_num)) and (float(p_row[2]) == float(row[2][5:]))):
              insertIntoDataStruct('STATION', float(row[2][5:]),CCS_sw_dict)
              insertIntoDataStruct('LAT',float(p_row[3]),CCS_sw_dict)
              insertIntoDataStruct('LONG',float(p_row[4]),CCS_sw_dict)
              insertIntoDataStruct('YEAR', int(float(s.join((row[3][:2],row[3][3:5])))),CCS_sw_dict)
              insertIntoDataStruct('MONTH', int(float(row[3][5:7])),CCS_sw_dict)
              insertIntoDataStruct('DAY', int(float(row[3][13:15])),CCS_sw_dict)

              # ALL DATA DICTIONARY
              for nvar in range(len(sw_vars)):
                  insertIntoDataStruct(sw_vars[nvar], np.float64(row[nrow[nvar]]),CCS_sw_dict)

# SAVE LINE DATA AS DICTIONARY
print 'MEEP'

filname = 'CalCOFI_LINE_' + str(line_num).zfill(5) + '.npy'

np.save(filname, CCS_sw_dict)


