import csv
import scipy
import numpy as np
import netCDF4 as nc
import pickle
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


# PUT DATA VALUES IN DICTIONARY 
def insertIntoDataStruct(name,val,aDict):
    if not name in aDict:
        aDict[name] = [(val)]
    else:
        aDict[name].append((val))

# ITERRATE OVER ROWS IN A CSV FILE
def dic_from_csv(rd_fil):
    # EXTRACT DATA FOR SPECIFIC STATIONS AND MONTHS 
    for row in rd_fil:
        # IF CORRECT STATION
        if (row[2] == stat_ids[wt] and \
           # IF CORRECT YEAR 
           # (97<=int(row[3][3:5])<=99) and \
           # IF CORRECT MONTH
           (2<int(row[3][5:7])<6) and \
           # (8<int(row[3][5:7])) and \
           (len(row[51])>0) and \
           (len(row[52])>0) and \
           (len(row[53])>0) and \
           (len(row[56])>0) and \
           (len(row[59])>0) and \
           (len(row[58])>0)):

           # ALL DATA DICTIONARY
           for nvar in range(len(sw_vars)):
               insertIntoDataStruct(sw_vars[nvar], np.float64(row[nrow[nvar]]),CCS_sw_dict)

           #if ((int(row[3][3:5]) == 97) and (8<int(row[3][5:7]))):# or \
              #((int(row[3][3:5]) == 99) and (int(row[3][5:7])<=3)) :

               # ISOPYCNAL RESTRICTION
               if ((np.float(row[53])>iso_pyc[wt]-mkt) and \
                   (np.float(row[53])<iso_pyc[wt]+pkt)):
                  for nvar in range(len(sw_vars)): 
                      insertIntoDataStruct(sw_vars[nvar], np.float64(row[nrow[nvar]]),CCS_sw_iso_dict) 

           #elif (int(row[3][3:5]) == 98):
           #   for nvar in range(len(sw_vars)): 
           #       insertIntoDataStruct(sw_vars[nvar], np.float64(row[nrow[nvar]]),CCS_sw_iso_dict)

def dic_from_nc(fid):
    st_id = fid.variables['sta_id'][:]
    dp_id = fid.variables['depth_id'][:]
    s=''
    for row in range(len(st_id)):
        if (s.join(st_id[row][:]) == stat_ids[wt] and \
           (97<=int(s.join(dp_id[row][3:5]))<99) and \
           (2<int(s.join((dp_id[row][5:7])))<7)):
           print s.join(dp_id[row][3:5]), s.join(dp_id[row][5:7])
           for nvar in range(len(sw_vars)):
               insertIntoDataStruct(sw_vars[nvar], fid.variables[nc_vars[nvar]][row],CCS_sw_dict)
                
               # ISOPYCNAL RESTRICTION
               if ((fid.variables[nc_vars[-1]][row]>iso_pyc[wt]-mkt) and \
                   (fid.variables[nc_vars[-1]][row]<iso_pyc[wt]+pkt)):
                  for nvar in range(len(sw_vars)):
                      insertIntoDataStruct(sw_vars[nvar], np.float64(fid.variables[nc_vars[nvar]][row]),CCS_sw_iso_dict)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SELECT STATION FOR SPECIFIC WATER TYPE
wt = 2

# SELECT DENSITY RANGE
pkt = .3  # Buffer for potential density (isopycnal) values 
mkt = .3

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# WATER MASS SAMPLES:
wt_lab = ['Pacific Sub Arctic Water', \
          'Pacific Equatorial Water', \
          'North Pacific Central Water']

wt_abs = ['PSA_80.80','PEW_93.30','NPCW_93.110']

stat_ids = ['080.0 080.0', '093.3 030.0','093.3 110.0']

iso_pyc = [25.8,26.5,26.5]

# VARIABLES USED FOR WATER TYPE DEF
sw_vars = ['PTEMP','SALINITY','OXYGEN','SILICATE','PHOSPHATE']

# DATA FILES
# SWFSC
# csv_fils = ['/Users/liz.drenkard/external_data/CalCOFI/194903-201402_Bottle.csv',\
#             '/Users/liz.drenkard/external_data/CalCOFI/201402-201701_Bottle.csv']
# nc_fil = ''
# MACBOOK
csv_fils = ['/Users/elizabethdrenkard/external_data/CalCOFI/194903-201402_Bottle.csv',\
            '/Users/elizabethdrenkard/external_data/CalCOFI/201402-201701_Bottle.csv']

nc_fil = '/Users/elizabethdrenkard/external_data/CalCOFI/download_CalCOFI_CTD_surveys_feb96_to_oct99.nc'

# CSV COLUMN CONTAINING VARIABLE IN SPREADSHEET
nrow    = [51,52,56,58,59]

# NCFIL VARIABLES
nc_vars = ['r_temp', 'r_salinity', 'r_o2', 'r_sio3', 'r_po4','r_sigma']

# DEFINE DICTIONARY
CCS_sw_dict = {}
CCS_sw_iso_dict = {}

# OPEN/READ CalCOFI DATA FILES
if wt ==1 :
   fid = nc.Dataset(nc_fil) 
   dic_from_nc(fid)   
else:
   # USING SPREADSHEET DATA
   for nf in range(len(csv_fils)):
       op_fil = open(csv_fils[nf],'rU')
       rd_fil = csv.reader(op_fil, delimiter=',')
       dic_from_csv(rd_fil)

# ANALYSIS AND PLOTTING
# PTEMP INCREASING MONOTONICLY
Ind = np.argsort(CCS_sw_iso_dict['PTEMP'])
PTEMP = np.array(CCS_sw_iso_dict['PTEMP'])[Ind]

npt = len(PTEMP)

print 'PTEMP MIN', np.min(PTEMP)
print 'PTEMP MAX', np.max(PTEMP)
print 'PTEMP MEAN', np.mean(PTEMP)
f3, ax3 = plt.subplots(2, 4,figsize=(12,5),sharey='row')

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
    PTEMP_orig = np.array(CCS_sw_iso_dict['PTEMP'])[Ind]
    # DEALING WITH NAN's
    PTEMP = PTEMP_orig[~np.isnan(VAR)]
    VAR = VAR[~np.isnan(VAR)]
    # DEGREES OF FREEDOM
    df = len(PTEMP)-(n+1)
    # Vandermonde matrix of x
    V = np.vander(PTEMP,n+1)
    # QR decomposition: 
    Q,R = np.linalg.qr(V)
    # ERROR for x
    err_fit = np.sqrt(1+np.sum(Q**2,axis=1))

    # FITTING POLYNOMIAL
    coeff = np.polyfit(PTEMP,VAR,n)

    # RESIDUALS CALCULATION
    r = VAR - np.dot(V,coeff)
    delta   = np.linalg.norm(r)/df*err_fit
    err[nv] = np.mean(delta)
     
    # SUBPLOTS
#    a = nv/2; b = -2*(nv/2) + nv

    ax3[0,nv].plot(VAR,PTEMP,'b+')
    ax3[0,nv].plot(np.polyval(coeff,PTEMP),PTEMP,'r')
    ax3[1,nv].plot(np.array(CCS_sw_dict[sw_vars[nv+1]]),np.array(CCS_sw_dict['PTEMP']),'o')
    ax3[1,nv].plot(np.array(CCS_sw_iso_dict[sw_vars[nv+1]])[Ind],np.array(CCS_sw_iso_dict['PTEMP'])[Ind],'bo',markeredgecolor='k')

    print '-------------------'
    print sw_vars[nv+1], 'ERROR = ', np.sum((np.polyval(np.polyfit(PTEMP, VAR, 1), PTEMP) - VAR)**2)
    print sw_vars[nv+1], 'MIN = ', np.min(np.polyval(coeff,PTEMP)) 
    print sw_vars[nv+1], 'MAX = ', np.max(np.polyval(coeff,PTEMP)) 
    print sw_vars[nv+1], 'MEAN = ', np.mean(VAR)
print err

# FIGURE DETAILS
# YTICKS
ytick_r1 = [[5,10,15],[5,11,17],[6,12,18]]
ytick_r2 = [[9,10,11,12],[7.5,9,10.5],[7,8,9]]

title_txt = 'CalCOFI LineStation ' + stat_ids[wt]

f3.subplots_adjust(hspace=0.2)
f3.suptitle(title_txt,fontsize=15)
f3.text(0.08, 0.5, 'Temperature ($^\circ$C)', ha='center', va='center', rotation='vertical',fontsize=13)

ax3[1,0].set_xlabel('Salinity (g/kg)',fontsize=12)
ax3[1,1].set_xlabel('Oxygen (ml/L)',fontsize=12)
ax3[1,2].set_xlabel('Silicate (umol/L)',fontsize=12)
ax3[1,3].set_xlabel('Phosphate (umol/L)',fontsize=12)

ax3[1,0].set_yticks(ytick_r1[wt])
ax3[0,0].set_yticks(ytick_r2[wt])

fig_tit = wt_abs[wt] + '_regression.png'
plt.savefig(fig_tit)

plt.show()

