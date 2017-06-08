# ######### OMP analysis main program version 2  ###################
#      
# omp2int.py 
#
# This is the interactive version of an easy-to-handle package for the use of 
# OMP analysis to resolve fractions of water masses involved in the
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
# sw_ptmp sw_dens0.m (Philip Morgan, CSIRO) may be called for some data files
# sw_dist.m (Philip Morgan, CSIRO) is called through the contour2 call
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
import numpy as np
from omp2 import omp2 #set up omp2.m as a py script

#clear all
#close all
print '  '
print 'OMP Analysis version 2 (March 1999)'
print '===================================  '
print '  '
print 'Note: Data sets for this program must contain the following information:'
print '  latitude: essential'
print '  longitude: essential'
print '  pressure: essential'
print '  salinity: essential'
print '  temperature: essential unless potential temperature is supplied'
print '  potential temperature: optional (will be calculated if not supplied)'
print '  density: optional (will be calculated if not supplied)'
print '  oxygen: optional'
print '  phosphate: optional'
print '  nitrate: optional'
print '  silicate: optional'
print '  potential vorticity: optional (will be calculated if necessary)'
print '===================================  '
print '  '
print 'Enter control values for this program run. Values in [] indicate default'
print 'values which will be used if no entry is supplied.'
print 'The run will issue a program run summary after successful completion.'
print 'Make sure that you retain a copy of the summary for later reference.'
print '  '
# choose basic or extended OMP (See the web manual for details)
OMP='cla'
incontrol = input('Do you want to apply basic or extended OMP analysis (b/e)?  [b]  ')

print '  '

if incontrol == 'e':
   OMP = 'ext'
   print 'YOU CHOSE TO USE EXTENDED OMP ANALYSIS.'
else:
   print 'YOU CHOSE TO USE BASIC OMP ANALYSIS.'

print '  '

#define your data set (this must be a *.mat file)
incontrol = input('Which data set do you want to use?  [testdata]  ')
if len(incontrol) > 0:
   dataset = incontrol
else:
   dataset = 'testdata'

print '  '

print 'YOU CHOSE THE DATASET:  ', dataset, '.'
eval('load dataset') # LD: Tricky if mat file - need to sort out file format

if not (('temp' in locals()) & ('ptemp' in locals())): 
   print 'WARNING: This dataset does not contain a variable recognised as temperature!'
if not 'sal' in locals():
   print 'WARNING: This dataset does not contain a variable recognised as salinity!'
if not 'long' in locals():	
   print 'WARNING: This dataset does not contain a variable recognised as longitude!'
if not 'lat' in locals(): 
   print 'WARNING: This dataset does not contain a variable recognised as latitude!'
if not 'press' in locals(): 
   print 'WARNING: This dataset does not contain a variable recognised as pressure!'

eex[:11] = [0,0,0,0,0,0,0,0,0,0,0];   # index of available variables
esx[:11] = [0,0,0,0,0,0,0,0,0,0,0];   # index of selected variables
                                      # 1: latitude
	     		              # 2: longitude
				      # 3: pressure
                                      # 4: salinity
				      # 5: potential temperature
				      # 6: oxygen
                                      # 7: phosphate
				      # 8: nitrate
				      # 9: silicate
                                      # 10: potential vorticity
				      # 11: temperature

# NOTE: For historical reasons the two columns mass conservation and potential vorticity are
# swapped in the program so that mass conservation is always the last column, after potential vorticity.
# The arrangement of the water type matrix and the weight vector thus differs from the description
# in the user manual. This should not be of concern but has to be watched when changing the code.
		
print 'This dataset contains the following variables:'
if 'lat' in locals():
   print '  latitude'
   eex[0] = 1
if 'long' in locals():
   print '  longitude' 
   eex[1] = 1
if 'press' in locals():
   print '  pressure' 
   eex[2] = 1
if 'temp' in locals():
   print '  temperature' 
   eex[3] = 1
else:
   temp = sw_temp(sal,ptemp,press,0)

eex[10] = 1
if 'sal' in locals():
   print '  salinity' 
   eex[3] = 1
if 'ptemp' in locals():
   print '  potential temperature'
   eex[4] = 1

if 'pdens' in locals():
   print '  density'
if 'oxy' in locals():
   print '  oxygen'
   eex[5] = 1
if 'ph' in locals():
   print '  phosphate'
   eex[6] = 1
if 'ni' in locals():
   print '  nitrate'
   eex[7] = 1
if 'si' in locals():
   print '  silicate'
   eex[8] = 1
if 'pvort' in locals():
   print '  potential vorticity' 
   eex[9] = 1

print '  '
if not 'ptemp' in locals():
   '  potential temperature is calculated'
if not 'pdens' in locals():
   print '  density is calculated'

#if exist('pvort') == 0

switchpot = 'n'
switchpot = input('Do you want to use potential vorticity in the analysis (y/n)? [n]  ')
#### RESUME PYTHON
if ((switchpot == 'y') & (eex[9]!=1)):
	print 'Potential vorticity will be calculated and included'
else:
	print 'Potential vorticity will not be included'
#end

# Sort out data through specific criteria; set the depth range
# (This assumes that negative oxygen and nutrient data indicate missing data.)

print '  '
print 'Specify a range for the analysis. For example '
print 'using only data in the density range 23 and 28 '
print 'with oxygen larger then 20 write:'
print 'pdens>=23&pdens<=28&oxy>=20'
print '  '

selection = 'press>=0'  # (just in case one ignores the above field)

incontrol = input('typpe your selection here: ')

if isempty(incontrol):
   incontrol=selection
else:
   selection=incontrol


#Check and if necessary calculate potential vorticity
if ((switchpot == 'y') & (eex[9] != 1)):

#Find top and bottom pressure for each station, calculate potential vorticity

	statind=[0, np.tanspose(np.where(np.diff(press)<0)), len(press)] 
	vvort =[]
	pp = []
	bfrq,vort,p_ave = sw_bfrq(sal,temp,press,lat)
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
	del bfrq
	del vort
	del vvort
	del p_ave
	del pp
	eex[9] = 1; esx[9] = 1

if esx[9] == 1:
   pvort = abs(pvort)

nvar = 3 
esx = [1,1,1,1,1,0,0,0,0,0,0]
print '  '
print 'Specify the data you want to use [default is yes = included in the analysis]:'
print 'longitude:   yes'
print 'latitude:    yes'
print 'pressure:    yes'
print 'salinity:    yes'
print 'potential temperature: yes'
iox = 'y'
iph = 'y'
ini = 'y'
isi = 'y'
if eex[5] == 1:
	incontrol = input('oxygen (y/n):  [y]  ','s')
	if length(incontrol) > 0:
           iox = incontrol
	if iox == 'y': 
           nvar+= 1
           esx[5] = 1

if eex[6] == 1:
	incontrol = input('phosphate (y/n):  [y]  ')
	if length(incontrol) > 0:
	   iph = incontrol
	if iph == 'y': 
           nvar+= 1
           esx[6] = 1

if eex[7] == 1:
	incontrol = input('nitrate (y/n):  [y]  ')
	if length(incontrol) > 0:
		ini = incontrol
	if ini == 'y': 
           nvar += 1 
           esx[7] = 1

if eex[8] == 1:
	incontrol = input('silicate (y/n):  [y]  ')
	if length(incontrol) > 0:
		isi = incontrol
	if ~isempty(isi)&isi == 'y':
           nvar += 1 
           esx[8] = 1

nvar+= 1 
esx[9] = 1

# ****************************************
#  Specify the Weigthing Matrix (a .mat file; see manual for details on how to calculate weights.)
print '  '
incontrol = 'f'
incontrol = input('Do you want to enter weights manually or from a file (m/f)?  [file]  ','s')

if ((len(incontrol) == 0) | (incontrol == 'f')):
	incontrol = input('Which file do you want to use to read the weights?  [testwght]  ')
	if len(incontrol) > 0:
		weightset = incontrol
	else:
		weightset = 'testwght'
	eval('load weightset')
	
	# Check which weights are needed and reset the diagonal:
	A = np.diag(Wx)
	A1 = A[7]  # change order of weights so that mass conservation is last
	A[7] = A[6]
	A[6] = A1;
	if esx[4] == 0:
           A[0] = 0
	   ratio[0] = -99999			# no pot. temperature weight if not needed
	if esx[3] == 0:
           A[1] = 0
	   ratio[1] = -99999      		# no salinity weight if not needed
	if esx[5] == 0: 
           A[2] = 0
	   ratio[2] = -99999    		# no oxygen weight if no oxygen
	if esx[6] == 0: 
           A[3] = 0
	   ratio[3] = -99999			# no phosphate weight if no phosphate
	if esx[7] == 0: 
           A[4] = 0
	   ratio[4] = -99999			# no nitrate weight if no nitrate
	if esx[8] == 0:
           A[5] = 0
	   ratio[5] = -99999			# no silicate weight if no silicate
	if esx[9] == 0:
           A[6] = 0
	   ratio[6] = -99999			# no pot. vorticity weight if not needed
else:
	A = [0,0,0,0,0,0,0,0]
	ratio = [0,0,-99999,-99999,-99999,-99999,0,0]
	A[0] = input('Enter weight for potential temperature:  ')
	A[1] = input('Enter weight for salinity:  ')
	if ((eex[5] == 1) & (iox == 'y')):
           A[2] = input('Enter weight for oxygen:  ')
	if ((eex[6] == 1) & (iph == 'y')): 
           A[3] = input('Enter weight for phosphate:  ')
	if ((eex[7] == 1) & (ini == 'y')): 
           A[4] = input('Enter weight for nitrate:  ')
	if ((eex[8] == 1) & (isi == 'y')): 
           A[5] = input('Enter weight for silicate:  ')
	if eex[9] == 1: 
           A[6] = input('Enter weight for potential vorticity:  ')
	A[7] = input('Enter weight for mass conservation:  ')
	if OMP == 'ext':
           if ((eex[5] == 1) & (iox == 'y')): 
	      ratio[2] = input('Enter Redfield ratio for oxygen (recommended -170):  ') 
	   if ((eex[6] == 1) & (iph == 'y')): 
	      ratio[3] = input('Enter  Redfield ratio for phosphate (should be 1):  ') 
	   if ((eex[7] == 1) & (ini == 'y')): 
	      ratio[4] = input('Enter  Redfield ratio for nitrate (recommended 16):  ') 
	   if ((eex[8] == 1) & (isi == 'y')): 
	      ratio[5] = input('Enter  Redfield ratio for silicate (recommended 40):  ') 

statind = np.where(A>0)[0]
Wx = np.diag(A[statind])
statind = np.where(ratio>-99999)[0]
redfrat = ratio[statind]	 # Redfield ratio for selected variables only
print '  '
print 'Your weight matrix is:'
print '  '
print Wx
del A


# *************************************************
# Select source water types from file
incontrol = input('Which routine do you want to use to define source water types?  [qwt2]  ');
if length(incontrol) > 0:
   source = incontrol
else:
   source = 'qwt2'

#First, display all available water types

qwt_pos = [1,2]
[G0,wmnames,k] = eval('source(qwt_pos,0)')
qwt_pos = []
for i in range(k):
    qwt_pos = [qwt_pos, i]
del G1
[G0,wmnames,i] = eval('source (qwt_pos,1)')
print '  '
print 'Here is a list of the available water type definitions.'
print '  '
print 'Water mass names (one for each row):'
print '  '
print wmnamesq
print '  '
print 'Water type definitions for the selected variables and mass conservation'
print '  '
i = 2
G1[0,:] = G0[0,:]
G1[1,:] = G0[1,:]
if esx[5] == 1:
   G1[2,:] = G0[2,:]
   i = i+1
if esx[6] == 1:
   G1[i,:] = G0[3,:]
   i = i+1
if esx[7] == 1:
   G1[i,:] = G0[4,:]
   i = i+1
if esx[8] == 1:
   G1[i,:] = G0[5,:]
   i = i+1
if esx[9] == 1:
   G1[i,:] = abs(G0[7,:])
   i = i+1
G1[i,:] = G0[6,:]
print G1
print '  '

# Now select appropriate source water types

wm = 4
incontrol = input('How many water types do you want for your analysis?  [4]  ')
if length(incontrol) > 0:
   wm = incontrol 
print '(The default for the next entries is 1, 2, 3 etc.'
print 'up to the number of water types selected.)'
qwt_pos = []
for i in range (wm):
	k = i
	incontrol = input('Select water type (row) number: ')
	if length(incontrol) > 0:
           k = incontrol
	qwt_pos = [qwt_pos, k]

del G1
G0,wmnames,i = eval('source[qwt_pos,0]')
print '  '
print 'You selected the following water type definitions.'
print '  '
print 'Water mass names (one for each row):'
wm_index = []
wm_ind0  = []
wm_ind1  = []
j = 0
print '  '
tit_index = []
for i in range(len(qwt_pos)):
    wm_ind1 = wmnames[5*(qwt_pos[i]-1):5*(qwt_pos[i]-1)+5]
    print wmnames[5*(qwt_pos[i]-1):5*(qwt_pos[i]-1)+5]
    k = (wm_ind0==wm_ind1)
    if not k:
       j = j+1
       tit_index = [tit_index, wmnames[5*(qwt_pos[i]-1):5*(qwt_pos[i]-1)+5]]
    wm_ind0 = wm_ind1
    wm_index = [wm_index, j]

nr_of_wm = wm_index[len(wm_index)-1]

print '  '
print 'Selected water type definitions:'
print '  '

i = 2
#LD: ?????? # del G1 # BEST WAY ?!?
G1[0,:] = G0[0,:]
G1[1,:] = G0[1,:]
if esx[5] == 1:
   G1[2,:] = G0[2,:]
   i = i+1
if esx[6] == 1:
   G1[i,:] = G0[3,:]
   i = i+1
if esx[7] ==1:
   G1[i,:] = G0[4,:]
   i = i+1
if esx[8] == 1:
   G1[i,:] = G0[5,:]
   i = i+1
if esx[9] == 1:
   G1[i,:] = G0[7,:]
   i = i+1
G1[i,:] = G0[6,:]
print G1


# This is the main part of it all: The call to omp2.m which does the analysis
omp2

# It's all done. Documentation and display is all in omp2.m.
