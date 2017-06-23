# ######### OMP analysis main program version 2  ###################
# 
# omp2auto.m
#     
# This is the control file version of an easy-to-handle package for the use 
# of OMP analysis to resolve fractions of water masses involved in the
# mixing of water masses at a given point in the ocean. The original
# version was prepared by Johannes Karstensen. This version incorporates
# improvements by Matthias Tomczak.
#
# This program will run without any changes, using the default settings
# supplied for all necessary input, and produce output based on
# the data file testdata.mat supplied with this package. For details
# see the README.ps or README.html files.
#
# Some preparation work is necessary if you want to use the program with
# your own data and water type definitions. Again, details can be found
# in the README.ps or README.html files.
#
#
# Function calls used: qwt2.m qwt_tst.m nansum.m (Philip Morgan, CSIRO)
#  sw_ptmp sw_dens0.m (Philip Morgan, CSIRO) may be called for some data files
#  sw_dist.m (Philip Morgan, CSIRO) is called through the contour2 call
# ---------------------------------------------
# This program is part of the OMP package from:
# GEOMAR
# Helmholtz Centre for Ocean Res. Kiel  FIAMS, Flinders University
# J. Karstensen                         Matthias Tomczak
# Duesternbrooker Weg 20				GPO Box 2100
# 24106 Kiel                            Adelaide, SA
# Germany                               Australia
#
# BUGS: jkarstensen@geomar.de
#   or  matthias.tomczak@flinders.edu.au
# --------------------------------------------
import scipy.io
import numpy as np
import matplotlib.pyplot as plt
from omp2 import omp2 #set up omp2.m as a py script
from qwt2 import qwt2
from sw_dist import sw_dist

print '  '
print 'OMP Analysis version 2 (March 1999)'
print '===================================  '
print '  '


# Loading the run parameters from the control file (call to incontr2)

### FROM incontr2: #### 
OMP = 'cla'
dataset = '/Users/elizabethdrenkard/TOOLS/omp2mats/testdata.mat'
selection=  '(pdens>26.3) & (pdens <27) & (oxy>=20) & (press>300) & (press<500)' 
#'(press>300) & (press<500)'
#'(pdens>26.3) & (pdens <27) & (oxy>=20) & (press>300) & (press<500)'
switchpot = 'n'
iox = 'y' # oxygen switch
iph = 'y' # phosphate switch
ini = 'y' # nitrate switch
isi = 'n' # silicate switch
weightset='/Users/elizabethdrenkard/TOOLS/omp2mats/testwght.mat'
swtypes = 'qwt2'
wm = 2
qwt_pos = [0,1,4,5] # changed from [1,2,3,4]
#####################

mat_dat=scipy.io.loadmat(dataset)
globals().update(mat_dat)

weight_dat=scipy.io.loadmat(weightset)
globals().update(weight_dat)

fig, ax = plt.subplots()
dist,phaseangle = sw_dist(lat,long,'km')
cumdist=np.append(0, np.cumsum(dist))
ax.plot(cumdist,pdens.squeeze(),'ko')
ax.invert_yaxis()
#plt.show()
#####################

#eex[:11] = [1,1,1,1,1,0,0,0,0,0,1]   # index of available variables
#esx[:11] = [1,1,1,1,1,0,0,0,0,0,1]   # index of selected variables

eex = np.array((1,1,1,1,1,0,0,0,0,0,1))   # index of available variables
esx = np.array((1,1,1,1,1,0,0,0,0,0,1))   # index of selected variables
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


if mat_dat.has_key('lat'):
   eex[0] = 1

if mat_dat.has_key('long'):
   eex[1] = 1

if mat_dat.has_key('press'):
   eex[2] = 1	

if mat_dat.has_key('sal'):
   eex[3] = 1

if mat_dat.has_key('ptemp'):
   eex[4] = 1

if mat_dat.has_key('oxy'):
   eex[5] = 1

if mat_dat.has_key('ph'):
   eex[6] = 1

if mat_dat.has_key('ni'):
   eex[7] = 1

if mat_dat.has_key('si'):
   eex[8] = 1

if mat_dat.has_key('pvort'):
   eex[9] = 1

if mat_dat.has_key('temp'):
   temp = mat_dat['temp']
else:
   temp = sw_temp(sal,ptemp,press,0)



#Check and if necessary calculate potential vorticity
if switchpot == 'y':
#Find top and bottom pressure for each station, calculate potential vorticity
	statind=[0, np.tanspose(np.where(np.diff(press)<0)), len(press)] 
	print 'gone through all right'
	vvort =[]
	pp = []
        bfrq,vort,p_ave = sw_bfrq(sal,temp,press,lat);
        for i in range(len(vort[:])):
            vvort = [vvort, vort[i]]
            pp    = [pp, p_ave[i]]
	vvort = 10E08*[vvort, 0]
	pp    = [pp, 10000]
	pvort = -999999*np.ones(press.shape)
        for i in range(1,len(statind[:])):
            pvort[statind[i-2]+2:statind[i]-1] = \
            interp1(pp[statind[i-1]:statind[i]-1],vvort[statind[i-1]:statind[i]-1],\
                        press[statind[i-1]+2:statind(i)-1])

	eex[9] = 1;
	del bfrq
	del vort
	del vvort
	del p_ave
	del pp
	pvort = abs(pvort)


# Determine the number of variables used in this run:
nvar = 3
if iox == 'y':
   nvar += 1 
   esx[5] = 1

if iph == 'y':
   nvar += 1
   esx[6] = 1

if ini == 'y':
   nvar += 1
   esx[7] = 1

if isi == 'y':
   nvar += 1 
   esx[8] = 1

if switchpot == 'y': 
   nvar += 1 
   esx[9] = 1

# Read the weight and Redfield ratio file
    

# Check which weights are needed and reset the diagonal:
A    = np.diag(Wx)
A.setflags(write=1)
A1   = A[7]  # change order of weights so that mass conservation is last
A[7] = A[6]
A[6] = A1
ratio = ratio.squeeze()

if esx[4] == 0:
   A[0] = 0
   ratio[0] = -99999                    # no pot. temperature weight if not needed
if esx[3] == 0:
   A[1] = 0
   ratio[1] = -99999                    # no salinity weight if not needed
if esx[5] == 0:
   A[2] = 0
   ratio[2] = -99999                    # no oxygen weight if no oxygen
if esx[6] == 0:
   A[3] = 0
   ratio[3] = -99999                    # no phosphate weight if no phosphate
if esx[7] == 0:
   A[4] = 0
   ratio[4] = -99999                    # no nitrate weight if no nitrate
if esx[8] == 0:
   A[5] = 0
   ratio[5] = -99999                    # no silicate weight if no silicate
if esx[9] == 0:
   A[6] = 0
   ratio[6] = -99999                    # no pot. vorticity weight if not needed

statind = np.where(A>0)[0]
Wx = np.diag(A[statind])
statind = np.where(ratio>-99999)[0]
redfrat = ratio[statind]         # Redfield ratio for selected variables only
print '  '
del A
# End of if statements for weights and Redfield ratio

# Read the water types
G0,wmnames,i = qwt2(qwt_pos,0)

wm_index = []
wm_ind0  = []
wm_ind1  = []
j = 0
print '  '
tit_index = []
print wmnames
for i in range(len(qwt_pos)):
    wm_ind1 = wmnames[qwt_pos[i]]
    k = (wm_ind0==wm_ind1)
    if not k:
       j = j+1
       tit_index.extend([wmnames[qwt_pos[i]]])
    wm_ind0 = wm_ind1
    wm_index.extend([j])

nr_of_wm = wm_index[len(wm_index)-1]

i = 2
#del G1
# G1 add rows
G1 = G0[:2,:]

if esx[5] == 1:
   G1 = np.concatenate((G1,np.array([G0[2,:]])),axis=0)
   i = i+1
if esx[6] == 1:
   G1 = np.concatenate((G1,np.array([G0[3,:]])),axis=0)   
   i = i+1
if esx[7] ==1:
   G1 = np.concatenate((G1,np.array([G0[4,:]])),axis=0)
   i = i+1
if esx[8] == 1:
   G1 = np.concatenate((G1,np.array([G0[5,:]])),axis=0)
   i = i+1
if esx[9] == 1:
   G1 = np.concatenate((G1,np.array([G0[7,:]])),axis=0)   
   i = i+1
G1 = np.concatenate((G1,np.array([G0[6,:]])),axis=0)
# This is the main part of it all: The call to omp2.m which does the analysis
omp2(OMP,nr_of_wm,tit_index,qwt_pos,wmnames,Wx,lat,switchpot,selection,long,esx,press,sal,oxy,ptemp,temp,pdens,ph,ni,G1,wm_index)
# It's all done. Documentation and display is all in omp2.m.
