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

clear all
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
OMP='cla';
#### PICK UP WITH PYTHON
incontrol = input('Do you want to apply basic or extended OMP analysis (b/e)?  [b]  ','s');
print '  ')

switch(incontrol)
case 'e'
OMP = 'ext';
	print 'YOU CHOSE TO USE EXTENDED OMP ANALYSIS.'
otherwise
	print 'YOU CHOSE TO USE BASIC OMP ANALYSIS.'
end
print '  '

#define your data set (this must be a *.mat file)
incontrol = input('Which data set do you want to use?  [testdata]  ','s');
if length(incontrol) > 0
	dataset = incontrol
else
	dataset = 'testdata'
end

print '  '
print ['YOU CHOSE THE DATASET:  ' dataset '.']
eval(['load ' dataset])
if exist('temp') == 0 & exist('ptemp') == 0
	print 'WARNING: This dataset does not contain a variable recognised as temperature!'
end
if exist('sal') == 0
	print 'WARNING: This dataset does not contain a variable recognised as salinity!'
end
if exist('long') == 0
	print 'WARNING: This dataset does not contain a variable recognised as longitude!'
end
if exist('lat') == 0
	print 'WARNING: This dataset does not contain a variable recognised as latitude!'
end
if exist('press') == 0
	print 'WARNING: This dataset does not contain a variable recognised as pressure!'
end
eex(1:11) = [0 0 0 0 0 0 0 0 0 0 0];   # index of available variables
esx(1:11) = [0 0 0 0 0 0 0 0 0 0 0];   # index of selected variables
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
		
print 'This dataset contains the following variables:')
if exist('lat')   == 1 print '  latitude'); eex(1) = 1; end
if exist('long')  == 1 print '  longitude'); eex(2) = 1;  end
if exist('press') == 1 print '  pressure'); eex(3) = 1;  end
if exist('temp')  == 1
	print '  temperature');
else
	temp = sw_temp(sal,ptemp,press,0);
end
eex(11) = 1;
if exist('sal')   == 1 print '  salinity'; eex(4) = 1;  end
if exist('ptemp') == 1 print '  potential temperature'; eex(5) = 1;  end
if exist('pdens') == 1 print '  density'; end
if exist('oxy')   == 1 print '  oxygen'; eex(6) = 1;  end
if exist('ph')    == 1 print '  phosphate'; eex(7) = 1;  end
if exist('ni')    == 1 print '  nitrate'; eex(8) = 1;  end
if exist('si')    == 1 print '  silicate'; eex(9) = 1;  end
if exist('pvort') == 1 print '  potential vorticity'; eex(10) = 1;  end
print '  '
if exist('ptemp') == 0 print '  potential temperature is calculated'; end
if exist('pdens') == 0 print '  density is calculated'; end

#if exist('pvort') == 0

switchpot = 'n';
switchpot = input('Do you want to use potential vorticity in the analysis (y/n)? [n]  ','s');
if ~isempty(switchpot) & switchpot == 'y' & eex(10)~=1
	print 'Potential vorticity will be calculated and included';
else
	print 'Potential vorticity will not be included';
end
#end

# Sort out data through specific criteria; set the depth range
# (This assumes that negative oxygen and nutrient data indicate missing data.)

print '  '
print 'Specify a range for the analysis. For example '
print 'using only data in the density range 23 and 28 '
print 'with oxygen larger then 20 write:'
print 'pdens>=23&pdens<=28&oxy>=20'
print '  '

selection='press>=0';  # (just in case one ignores the above field)

incontrol= input('type your selection here: ','s');

if isempty(incontrol)
 incontrol=selection;
else
  selection=incontrol;
end


#Check and if necessary calculate potential vorticity
if ~isempty(switchpot)&switchpot == 'y' &eex(10)~=1

#Find top and bottom pressure for each station, calculate potential vorticity

	statind=[0 find(diff(press)<0).T length(press)]; 
	vvort =[];
	pp = [];
		[bfrq,vort,p_ave] = sw_bfrq(sal,temp,press,lat);
		for i = 1:size(vort(:))
			vvort = [vvort vort(i)];
			pp    = [pp p_ave(i)];
		end
	vvort = 10E08*[vvort 0];
	pp    = [pp 10000];
	pvort = -999999*ones(size(press));
	for i = 2:size(statind(:)):
		pvort(statind(i-1)+2:statind(i)-1) = ...
			interp1(pp(statind(i-1)+1:statind(i)-1),vvort(statind(i-1)+1:statind(i)-1),...
			press(statind(i-1)+2:statind(i)-1));
	end
	clear bfrq
	clear vort
	clear vvort
	clear p_ave
	clear pp
	eex(10) = 1; esx(10) = 1;
end 
if esx(10) == 1:
   pvort = abs(pvort)

nvar = 3; esx = [1 1 1 1 1 0 0 0 0 0 0];
print '  '
print 'Specify the data you want to use [default is yes = included in the analysis]:'
print 'longitude:   yes'
print 'latitude:    yes'
print 'pressure:    yes'
print 'salinity:    yes'
print 'potential temperature: yes';
iox = 'y';
iph = 'y';
ini = 'y';
isi = 'y';
if eex(6) == 1:
	incontrol = input('oxygen (y/n):  [y]  ','s');
	if length(incontrol) > 0:
           iox = incontrol
	if iox == 'y': 
           nvar+= 1
           esx(6) = 1

if eex(7) == 1
	incontrol = input('phosphate (y/n):  [y]  ','s')
	if length(incontrol) > 0:
	   iph = incontrol
	if iph == 'y' 
           nvar+= 1
           esx(7) = 1

if eex(8) == 1
	incontrol = input('nitrate (y/n):  [y]  ','s');
	if length(incontrol) > 0
		ini = incontrol;
	end
	if ini == 'y' nvar = nvar +1; esx(8) = 1; end
end
if eex(9) == 1
	incontrol = input('silicate (y/n):  [y]  ','s');
	if length(incontrol) > 0
		isi = incontrol;
	end
	if ~isempty(isi)&isi == 'y' nvar = nvar +1; esx(9) = 1; end
end

switch( switchpot)
case 'y' 
nvar+= 1 
esx(10) = 1;

# ****************************************
#  Specify the Weigthing Matrix (a .mat file; see manual for details on how to calculate weights.)
print '  '
incontrol = 'f'
incontrol = input('Do you want to enter weights manually or from a file (m/f)?  [file]  ','s');

if length(incontrol) == 0 | incontrol == 'f'
	incontrol = input('Which file do you want to use to read the weights?  [testwght]  ','s');
	if length(incontrol) > 0
		weightset = incontrol;
	else
		weightset = 'testwght';
	end
	eval(['load ' weightset]);
	
	# Check which weights are needed and reset the diagonal:
	A = diag(Wx);
	A1 = A(8);  # change order of weights so that mass conservation is last
	A(8) = A(7);
	A(7) = A1;
	if esx(5) == 0 A(1) = 0;
	ratio(1) = -99999; end			# no pot. temperature weight if not needed
	if esx(4) == 0 A(2) = 0;
	ratio(2) = -99999; end			# no salinity weight if not needed
	if esx(6) == 0 A(3) = 0;
	ratio(3) = -99999; end			# no oxygen weight if no oxygen
	if esx(7) == 0 A(4) = 0;
	ratio(4) = -99999;  end			# no phosphate weight if no phosphate
	if esx(8) == 0 A(5) = 0;
	ratio(5) = -99999;  end			# no nitrate weight if no nitrate
	if esx(9) == 0 A(6) = 0;
	ratio(6) = -99999;  end			# no silicate weight if no silicate
	if esx(10) == 0 A(7) = 0;
	ratio(7) = -99999;  end			# no pot. vorticity weight if not needed
else
	A = [0 0 0 0 0 0 0 0];
	ratio = [0  0  -99999  -99999  -99999  -99999  0  0];
	A(1) = input('Enter weight for potential temperature:  ');
	A(2) = input('Enter weight for salinity:  ');
	if (eex(6) == 1 & iox == 'y') A(3) = input('Enter weight for oxygen:  '); end
	if (eex(7) == 1 & iph == 'y') A(4) = input('Enter weight for phosphate:  ');  end
	if (eex(8) == 1 & ini == 'y') A(5) = input('Enter weight for nitrate:  ');  end
	if (eex(9) == 1 & isi == 'y') A(6) = input('Enter weight for silicate:  ');  end
	if eex(10) == 1 A(7) = input('Enter weight for potential vorticity:  ');  end
	A(8) = input('Enter weight for mass conservation:  ');
	if OMP == 'ext'
		if (eex(6) == 1 & iox == 'y') 
		ratio(3) = input('Enter Redfield ratio for oxygen (recommended -170):  '); 
		end
		if (eex(7) == 1 & iph == 'y') 
		ratio(4) = input('Enter  Redfield ratio for phosphate (should be 1):  ');  
		end
		if (eex(8) == 1 & ini == 'y') 
		ratio(5) = input('Enter  Redfield ratio for nitrate (recommended 16):  ');  
		end
		if (eex(9) == 1 & isi == 'y') 
		ratio(6) = input('Enter  Redfield ratio for silicate (recommended 40):  ');  
		end
	end
end

statind = find(A>0);
Wx = diag(A(statind))
statind = find(ratio>-99999);
redfrat = ratio(statind);	 # Redfield ratio for selected variables only
print '  '
print 'Your weight matrix is:'
print '  '
print Wx
clear A


# *************************************************
# Select source water types from file
incontrol = input('Which routine do you want to use to define source water types?  [qwt2]  ','s');
if length(incontrol) > 0:
	source = incontrol
else
	source = 'qwt2'



#First, display all available water types

qwt_pos = [1 2];
[G0,wmnames,k] = eval([source '(qwt_pos,0)']);
qwt_pos = [];
for i in range(k):
	qwt_pos = [qwt_pos i];
clear G1;
[G0,wmnames,i] = eval([source '(qwt_pos,1)']);
print '  '
print 'Here is a list of the available water type definitions.'
print '  '
print 'Water mass names (one for each row):'
print '  '
print wmnamesq
print '  '
print 'Water type definitions for the selected variables and mass conservation')
print '  '
i = 3;
G1(1,:) = G0(1,:);
G1(2,:) = G0(2,:);
if esx(6) == 1
	G1(3,:) = G0(3,:);
	i = i+1;
if esx(7) == 1
	G1(i,:) = G0(4,:);
	i = i+1;
if esx(8) == 1
	G1(i,:) = G0(5,:);
	i = i+1;
if esx(9) == 1
	G1(i,:) = G0(6,:);
	i = i+1;
if esx(10) == 1
	G1(i,:) = abs(G0(8,:));
	i = i+1;
G1(i,:) = G0(7,:);
print G1)
print '  ')

% Now select appropriate source water types

wm = 4;
incontrol = input('How many water types do you want for your analysis?  [4]  ');
if length(incontrol) > 0 wm = incontrol; 
print '(The default for the next entries is 1, 2, 3 etc.');
print 'up to the number of water types selected.)')
qwt_pos = [];
for i=1:wm
	k = i;
	incontrol = input('Select water type (row) number: ');
	if length(incontrol) > 0 k = incontrol; end
	qwt_pos = [qwt_pos k];
end

clear G1;
[G0,wmnames,i] = eval([source '(qwt_pos,1)']);
print '  ')
print 'You selected the following water type definitions.')
print '  ')
print 'Water mass names (one for each row):')
wm_index = [];
wm_ind0  = [     ];
wm_ind1  = [     ];
j = 0;
print '  ')
tit_index = [];
for i = 1:length(qwt_pos)
	wm_ind1 = wmnames(5*(qwt_pos(i)-1)+1:5*(qwt_pos(i)-1)+5);
	print wmnames(5*(qwt_pos(i)-1)+1:5*(qwt_pos(i)-1)+5))
	k = strcmp(wm_ind0,wm_ind1);
	if k == 0
		j = j+1;
		tit_index = [tit_index wmnames(5*(qwt_pos(i)-1)+1:5*(qwt_pos(i)-1)+5)];
		end
	wm_ind0 = wm_ind1;
	wm_index = [wm_index j];
end
nr_of_wm = wm_index(length(wm_index));

print '  ')
print 'Selected water type definitions:')
print '  ')
i = 3;
clear G1;
G1(1,:) = G0(1,:);
G1(2,:) = G0(2,:);
if esx(6) == 1
	G1(3,:) = G0(3,:);
	i = i+1;
end
if esx(7) == 1
	G1(i,:) = G0(4,:);
	i = i+1;
end
if esx(8) ==1
	G1(i,:) = G0(5,:);
	i = i+1;
end
if esx(9) == 1
	G1(i,:) = G0(6,:);
	i = i+1;
end
if esx(10) == 1
	G1(i,:) = G0(8,:);
	i = i+1;
end
G1(i,:) = G0(7,:);
print G1


# This is the main part of it all: The call to omp2.m which does the analysis
omp2

# It's all done. Documentation and display is all in omp2.m.
