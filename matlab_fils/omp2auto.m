% ######### OMP analysis main program version 2  ###################
% 
% omp2auto.m
%     
% This is the control file version of an easy-to-handle package for the use 
% of OMP analysis to resolve fractions of water masses involved in the
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


disp('  ')
disp('OMP Analysis version 2 (March 1999)')
disp('===================================  ')

% Loading the run parameters from the control file (call to incontr2)
eval(newfile);
eval(['load ' dataset])

eex(1:11) = [1 1 1 1 1 0 0 0 0 0 1];   % index of available variables
esx(1:11) = [1 1 1 1 1 0 0 0 0 0 1];   % index of selected variables
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
		
if exist('temp')  == 0
	temp = sw_temp(sal,ptemp,press,0);
end

if exist('sal')   == 1 eex(4) = 1;  end
if exist('ptemp') == 1 eex(5) = 1;  end
if exist('oxy')   == 1 eex(6) = 1;  end
if exist('ph')    == 1 eex(7) = 1;  end
if exist('ni')    == 1 eex(8) = 1;  end
if exist('si')    == 1 eex(9) = 1;  end
if exist('pvort') == 1 eex(10) = 1;  end

%Check and if necessary calculate potential vorticity
if switchpot == 'y'
%Find top and bottom pressure for each station, calculate potential vorticity
	statind=[0 find(diff(press)<0)' length(press)]; 
	disp('gone through all right')
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
	eex(10) = 1;
	clear bfrq
	clear vort
	clear vvort
	clear p_ave
	clear pp
	pvort = abs(pvort);
end


% Determine the number of variables used in this run:
nvar = 3;
if iox == 'y' nvar = nvar +1; esx(6) = 1; end
if iph == 'y' nvar = nvar +1; esx(7) = 1; end
if ini == 'y' nvar = nvar +1; esx(8) = 1; end
if isi == 'y' nvar = nvar +1; esx(9) = 1; end
if switchpot == 'y' nvar = nvar +1; esx(10) = 1; end

% Read the weight and Redfield ratio file
eval(['load ' weightset]);

%Check which weights are needed and reset the diagonal:
A = diag(Wx);
A1 = A(8);  % change order of weights so that mass conservation is last
A(8) = A(7);
A(7) = A1;
if esx(5) == 0 A(1) = 0;
	ratio(1) = -999; end			% no pot. temperature weight if not needed
if esx(4) == 0 A(2) = 0;
	ratio(2) = -999; end			% no salinity weight if not needed
if esx(6) == 0 A(3) = 0;
	ratio(3) = -999; end			% no oxygen weight if no oxygen
if esx(7) == 0 A(4) = 0;
	ratio(4) = -999;  end			% no phosphate weight if no phosphate
if esx(8) == 0 A(5) = 0;
	ratio(5) = -999;  end			% no nitrate weight if no nitrate
if esx(9) == 0 A(6) = 0;
	ratio(6) = -999;  end			% no silicate weight if no silicate
if esx(10) == 0 A(7) = 0;
	ratio(7) = -999;  end			% no pot. vorticity weight if not needed
statind = find(A>0);
Wx = diag(A(statind));
statind = find(ratio>-999);
redfrat = ratio(statind);	 % Redfield ratio for selected variables only
disp('  ')
clear A
% End of if statements for weights and Redfield ratio

% Read the water types
clear G1;
[G0,wmnames,i]=eval([swtypes '(qwt_pos,1)']);
wm_index = [];
wm_ind0  = [     ];
wm_ind1  = [     ];
j = 0;
tit_index = [];
for i = 1:length(qwt_pos)
	wm_ind1 = wmnames(5*(qwt_pos(i)-1)+1:5*(qwt_pos(i)-1)+5);
	k = strcmp(wm_ind0,wm_ind1);
	if k == 0
		j = j+1;
		tit_index = [tit_index wmnames(5*(qwt_pos(i)-1)+1:5*(qwt_pos(i)-1)+5)];
		end
	wm_ind0 = wm_ind1;
	wm_index = [wm_index j];
end
nr_of_wm = wm_index(length(wm_index));

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


close all

% This is the main part of it all: The call to omp2.m which does the analysis
omp2
% It's all done. Documentation and display is all in omp2.m.
