% ######### OMP analysis main program version 2  ###################
%      
% omp2int.m 
%
% This is the interactive version of an easy-to-handle package for the use of 
% OMP analysis to resolve fractions of water masses involved in the
% mixing of water masses at a given point in the ocean. The original
% version was prepared by Johannes Karstensen. This version incorporates
% improvements by Matthias Tomczak.
%
% This program will run without any changes, using the default settings
% supplied for all necessary input, and produce output based on
% the data file testdata.mat supplied with this package. For details
% see the README.ps or README.html files.
%
% Some preparation work is necessary if you want to use the program with
% your own data and water type definitions. Again, details can be found
% in the README.ps or README.html files.
%
%
% Function calls used: qwt2.m qwt_tst.m nansum.m (Philip Morgan, CSIRO)
%  sw_ptmp sw_dens0.m (Philip Morgan, CSIRO) may be called for some data files
%  sw_dist.m (Philip Morgan, CSIRO) is called through the contour2 call
%---------------------------------------------
% This program is part of the OMP package from:
% GEOMAR
% Helmholtz Centre for Ocean Res. Kiel  FIAMS, Flinders University
% J. Karstensen                         Matthias Tomczak
% Duesternbrooker Weg 20				GPO Box 2100
% 24106 Kiel                            Adelaide, SA
% Germany                               Australia
%
% BUGS: jkarstensen@geomar.de
%   or  matthias.tomczak@flinders.edu.au
%--------------------------------------------

clear all
close all
disp('  ')
disp('OMP Analysis version 2 (March 1999)')
disp('===================================  ')
disp('  ')
disp('Note: Data sets for this program must contain the following information:')
disp('  latitude: essential')
disp('  longitude: essential')
disp('  pressure: essential')
disp('  salinity: essential')
disp('  temperature: essential unless potential temperature is supplied')
disp('  potential temperature: optional (will be calculated if not supplied)')
disp('  density: optional (will be calculated if not supplied)')
disp('  oxygen: optional')
disp('  phosphate: optional')
disp('  nitrate: optional')
disp('  silicate: optional')
disp('  potential vorticity: optional (will be calculated if necessary)')
disp('===================================  ')
disp('  ')
disp('Enter control values for this program run. Values in [] indicate default')
disp('values which will be used if no entry is supplied.')
disp('The run will issue a program run summary after successful completion.')
disp('Make sure that you retain a copy of the summary for later reference.')
disp('  ')
% choose basic or extended OMP (See the web manual for details)
OMP='cla';
incontrol = input('Do you want to apply basic or extended OMP analysis (b/e)?  [b]  ','s');
disp('  ')

switch(incontrol)
case 'e'
OMP = 'ext';
	disp('YOU CHOSE TO USE EXTENDED OMP ANALYSIS.')
otherwise
	disp('YOU CHOSE TO USE BASIC OMP ANALYSIS.')
end
disp('  ')

% define your data set (this must be a *.mat file)
incontrol = input('Which data set do you want to use?  [testdata]  ','s');
if length(incontrol) > 0
	dataset = incontrol;
else
	dataset = 'testdata';
end

disp('  ')
disp(['YOU CHOSE THE DATASET:  ' dataset '.'])
eval(['load ' dataset])
if exist('temp') == 0 & exist('ptemp') == 0
	disp('WARNING: This dataset does not contain a variable recognised as temperature!')
end
if exist('sal') == 0
	disp('WARNING: This dataset does not contain a variable recognised as salinity!')
end
if exist('long') == 0
	disp('WARNING: This dataset does not contain a variable recognised as longitude!')
end
if exist('lat') == 0
	disp('WARNING: This dataset does not contain a variable recognised as latitude!')
end
if exist('press') == 0
	disp('WARNING: This dataset does not contain a variable recognised as pressure!')
end
eex(1:11) = [0 0 0 0 0 0 0 0 0 0 0];   % index of available variables
esx(1:11) = [0 0 0 0 0 0 0 0 0 0 0];   % index of selected variables
                                        % 1: latitude
										% 2: longitude
										% 3: pressure
                                        % 4: salinity
										% 5: potential temperature
										% 6: oxygen
                                        % 7: phosphate
										% 8: nitrate
										% 9: silicate
                                        %10: potential vorticity
										%11: temperature

% NOTE: For historical reasons the two columns mass conservation and potential vorticity are
% swapped in the program so that mass conservation is always the last column, after potential vorticity.
% The arrangement of the water type matrix and the weight vector thus differs from the description
% in the user manual. This should not be of concern but has to be watched when changing the code.
		
disp('This dataset contains the following variables:')
if exist('lat')   == 1 disp('  latitude'); eex(1) = 1; end
if exist('long')  == 1 disp('  longitude'); eex(2) = 1;  end
if exist('press') == 1 disp('  pressure'); eex(3) = 1;  end
if exist('temp')  == 1
	disp('  temperature');
else
	temp = sw_temp(sal,ptemp,press,0);
end
eex(11) = 1;
if exist('sal')   == 1 disp('  salinity'); eex(4) = 1;  end
if exist('ptemp') == 1 disp('  potential temperature'); eex(5) = 1;  end
if exist('pdens') == 1 disp('  density'); end
if exist('oxy')   == 1 disp('  oxygen'); eex(6) = 1;  end
if exist('ph')    == 1 disp('  phosphate'); eex(7) = 1;  end
if exist('ni')    == 1 disp('  nitrate'); eex(8) = 1;  end
if exist('si')    == 1 disp('  silicate'); eex(9) = 1;  end
if exist('pvort') == 1 disp('  potential vorticity'); eex(10) = 1;  end
disp('  ')
if exist('ptemp') == 0 disp('  potential temperature is calculated'); end
if exist('pdens') == 0 disp('  density is calculated'); end

%if exist('pvort') == 0
switchpot = 'n';
switchpot = input('Do you want to use potential vorticity in the analysis (y/n)? [n]  ','s');
if ~isempty(switchpot) & switchpot == 'y' & eex(10)~=1
	disp('Potential vorticity will be calculated and included');
else
	disp('Potential vorticity will not be included');
end
%end

% Sort out data through specific criteria; set the depth range
% (This assumes that negative oxygen and nutrient data indicate missing data.)

disp('  ')
disp('Specify a range for the analysis. For example ');
disp('using only data in the density range 23 and 28 ')
disp('with oxygen larger then 20 write:')
disp('pdens>=23&pdens<=28&oxy>=20')
disp('  ')

selection='press>=0';  % (just in case one ignores the above field)

incontrol= input('type your selection here: ','s');

if isempty(incontrol)
 incontrol=selection;
else
  selection=incontrol;
end


%Check and if necessary calculate potential vorticity
if ~isempty(switchpot)&switchpot == 'y' &eex(10)~=1

%Find top and bottom pressure for each station, calculate potential vorticity

	statind=[0 find(diff(press)<0)' length(press)]; 
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
	for i = 2:size(statind(:))
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
if esx(10) == 1 pvort = abs(pvort); end

nvar = 3; esx = [1 1 1 1 1 0 0 0 0 0 0];
disp('  ')
disp('Specify the data you want to use [default is yes = included in the analysis]:')
disp('longitude:   yes')
disp('latitude:    yes')
disp('pressure:    yes')
disp('salinity:    yes')
disp('potential temperature: yes');
iox = 'y';
iph = 'y';
ini = 'y';
isi = 'y';
if eex(6) == 1
	incontrol = input('oxygen (y/n):  [y]  ','s');
	if length(incontrol) > 0
		iox = incontrol;
	end
	if iox == 'y' nvar = nvar +1; esx(6) = 1; end
end
if eex(7) == 1
	incontrol = input('phosphate (y/n):  [y]  ','s');
	if length(incontrol) > 0
		iph = incontrol;
	end
	if iph == 'y' nvar = nvar +1; esx(7) = 1; end
end
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
nvar = nvar +1; 
esx(10) = 1; 
end


%****************************************
% Specify the Weigthing Matrix (a .mat file; see manual for details on how to calculate weights.)
disp('  ')
incontrol = 'f';
incontrol = input('Do you want to enter weights manually or from a file (m/f)?  [file]  ','s');

if length(incontrol) == 0 | incontrol == 'f'
	incontrol = input('Which file do you want to use to read the weights?  [testwght]  ','s');
	if length(incontrol) > 0
		weightset = incontrol;
	else
		weightset = 'testwght';
	end
	eval(['load ' weightset]);
	
	%Check which weights are needed and reset the diagonal:
	A = diag(Wx);
	A1 = A(8);  % change order of weights so that mass conservation is last
	A(8) = A(7);
	A(7) = A1;
	if esx(5) == 0 A(1) = 0;
	ratio(1) = -99999; end			% no pot. temperature weight if not needed
	if esx(4) == 0 A(2) = 0;
	ratio(2) = -99999; end			% no salinity weight if not needed
	if esx(6) == 0 A(3) = 0;
	ratio(3) = -99999; end			% no oxygen weight if no oxygen
	if esx(7) == 0 A(4) = 0;
	ratio(4) = -99999;  end			% no phosphate weight if no phosphate
	if esx(8) == 0 A(5) = 0;
	ratio(5) = -99999;  end			% no nitrate weight if no nitrate
	if esx(9) == 0 A(6) = 0;
	ratio(6) = -99999;  end			% no silicate weight if no silicate
	if esx(10) == 0 A(7) = 0;
	ratio(7) = -99999;  end			% no pot. vorticity weight if not needed
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
redfrat = ratio(statind);	 % Redfield ratio for selected variables only
disp('  ')
disp('Your weight matrix is:')
disp('  ')
disp(Wx)
clear A


%*************************************************
% Select source water types from file
incontrol = input('Which routine do you want to use to define source water types?  [qwt2]  ','s');
if length(incontrol) > 0
	source = incontrol;
else
	source = 'qwt2';
end



% First, display all available water types

qwt_pos = [1 2];
[G0,wmnames,k] = eval([source '(qwt_pos,0)']);
qwt_pos = [];
for i=1:k
	qwt_pos = [qwt_pos i];
end
clear G1;
[G0,wmnames,i] = eval([source '(qwt_pos,1)']);
disp('  ')
disp('Here is a list of the available water type definitions.')
disp('  ')
disp('Water mass names (one for each row):')
disp('  ')
disp(wmnames)
disp('  ')
disp('Water type definitions for the selected variables and mass conservation')
disp('  ')
i = 3;
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
if esx(8) == 1
	G1(i,:) = G0(5,:);
	i = i+1;
end
if esx(9) == 1
	G1(i,:) = G0(6,:);
	i = i+1;
end
if esx(10) == 1
	G1(i,:) = abs(G0(8,:));
	i = i+1;
end
G1(i,:) = G0(7,:);
disp(G1)
disp('  ')

% Now select appropriate source water types

wm = 4;
incontrol = input('How many water types do you want for your analysis?  [4]  ');
if length(incontrol) > 0 wm = incontrol; end
disp('(The default for the next entries is 1, 2, 3 etc.');
disp('up to the number of water types selected.)')
qwt_pos = [];
for i=1:wm
	k = i;
	incontrol = input('Select water type (row) number: ');
	if length(incontrol) > 0 k = incontrol; end
	qwt_pos = [qwt_pos k];
end

clear G1;
[G0,wmnames,i] = eval([source '(qwt_pos,1)']);
disp('  ')
disp('You selected the following water type definitions.')
disp('  ')
disp('Water mass names (one for each row):')
wm_index = [];
wm_ind0  = [     ];
wm_ind1  = [     ];
j = 0;
disp('  ')
tit_index = [];
for i = 1:length(qwt_pos)
	wm_ind1 = wmnames(5*(qwt_pos(i)-1)+1:5*(qwt_pos(i)-1)+5);
	disp(wmnames(5*(qwt_pos(i)-1)+1:5*(qwt_pos(i)-1)+5))
	k = strcmp(wm_ind0,wm_ind1);
	if k == 0
		j = j+1;
		tit_index = [tit_index wmnames(5*(qwt_pos(i)-1)+1:5*(qwt_pos(i)-1)+5)];
		end
	wm_ind0 = wm_ind1;
	wm_index = [wm_index j];
end
nr_of_wm = wm_index(length(wm_index));

disp('  ')
disp('Selected water type definitions:')
disp('  ')
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
disp(G1)


% This is the main part of it all: The call to omp2.m which does the analysis
omp2

% It's all done. Documentation and display is all in omp2.m.
