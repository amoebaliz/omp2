% OMP2.M: OMP analysis main program version 2
%     
% This is an updated version of an easy-to-handle package for the use of 
% OMP analysis to resolve fractions of water masses involved in the
% mixing of water masses at a given point in the ocean. The original
% version was prepared by Johannes Karstensen. This version incorporates
% improvements by Matthias Tomczak.
%
% This program is called by omp2int.m, omp2gui.m and omp2auto.m and will not
% run before one of these programs is called and placed all necessary
% variables into the workspace.
%
%
% Function calls used: qwt2.m qwt_tst.m nansum.m (Philip Morgan, CSIRO)
%  sw_ptmp sw_dens0.m (Philip Morgan, CSIRO) may be called for some data files
%  sw_dist.m (Philip Morgan, CSIRO) is called through the contour2 call
%
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

disp('  ')
disp(['OMP analysis now running. ' num2str(length(lat)) ' data points found.'])
disp('  ')
starttime = clock;
gap=0;

%% which MatLab version is used?
vers=version;vers=vers(1);

% Cut the data to the selcted range 
% 
% set potential vorticity to positive values (independent of hemisphere) if required
disp('Screening the data and reducing them to the selected range.')
disp('  ')
switch(switchpot)
case 'y'
	eval(['index=find(imag(pvort)==0&pvort<100&' selection ');'])
	pvort = abs(pvort);
otherwise
	eval(['index=find(' selection  ');'])
end
      lat=   lat(index)';
    press= press(index)';
     long=  long(index)';
      sal=   sal(index)';
      
	if exist('temp') == 1
	      temp= temp(index)';
	end
	if ~isempty(switchpot) & switchpot == 'y'
	      pvort=pvort(index)';
	end
	 if exist('ptemp') == 1
	      ptemp=ptemp(index)';
	else
             ptemp = sw_ptmp(sal,temp,press,0);
	end
	if exist('pdens') == 1
              pdens=pdens(index)';
	else
		pdens = sw_dens0(sal,temp) - 1000;
	end
	if esx(6)  == 1 oxy=  oxy(index)';  end
	if esx(7)  == 1 ph =   ph(index)';  end
	if esx(8)  == 1 ni =   ni(index)';  end
	if esx(9)  == 1 si =   si(index)';  end
clear index
disp(['OMP analysis now running. ' num2str(length(lat)) ' data points to be analysed.'])
disp('  ')

[m,n]=size(G1); % n = number of water types, m = number of equations

% normalise the source water matrix (get meanG, get stdG for weighting):
[G,mG,stdG]=norm_qwt(G1);

% EXTENDED OMP switch:
switch(OMP)
case 'ext'
   % Adding Redfield ratio to the system, ratio comes from weight file
 	G1(1:m,n+1)=[redfrat(1:m)]';
	% normalisation of the ratios:
	%---------------------------------
 	for rr=1:(m-1)
  		G(rr,n+1)=redfrat(rr)*(max(G(rr,1:n))-min(G(rr,1:n)))...
           /(max(G1(rr,1:n))-min(G1(rr,1:n)));
	end
	G(m,n+1)=0;
end

% adding weights
G2=Wx*G;


gap=0;

%***********************************************************
% This is the main loop for each data point; k = point index
% First some initial settings
   err=zeros(m,length(lat))-nan; 

switch(OMP)
case 'ext';
  biogeo=zeros(1,length(lat))-nan; 
end



A(1:wm_index(length(wm_index)),1:length(lat)) = ...
                 zeros(wm_index(length(wm_index)),length(lat));

oxy_dat = [];
ph_dat  = [];
ni_dat  = [];
si_dat  = [];
pv_dat  = [];

%Vector of each datapoint (btst) is build here
for k=1:length(lat);
	% selecting the correct parameters
     p_dat	= press(k);
     t_dat	= ptemp(k); 
     s_dat	=   sal(k); 
     if esx(6) == 1 oxy_dat	=   oxy(k); end
     if esx(7) == 1  ph_dat	=    ph(k); end
     if esx(8) == 1  ni_dat	=    ni(k); end
	 if esx(9) == 1  si_dat    =    si(k); end
     if exist('pdens') == 1  pden_dat	= pdens(k); end
	 if esx(10) == 1 pv_dat = pvort(k); end
     kon=1;
	
	btst= [t_dat,s_dat,oxy_dat,ph_dat,ni_dat,si_dat,pv_dat,kon];

	if etime(clock,starttime) > 5
		disp([num2str(k) ' data points analysed so far.'])
		starttime = clock;
	end
%looking for GAP (indicated through NaN):
   index1=find(~isnan(btst));
   index0=find(isnan(btst));

	cutit=n;
% using extended OMP we need one parameter more 
% (because we have one unknown more!)
	if OMP(1:3)=='ext' cutit=n+1; end

	if length(index1) < cutit+1 %if1
% not enough parameters to find a NNLS fit
% DATA point not successful analysed
		disp(['ANALYSIS of the datapoint failed, not enough parameters available !!'])
   			A(1:nr_of_wm,k) = nan;
 			Dual(1:nr_of_wm,k) = nan;
 		gap=gap+1;

	else
     %new data without GAP:
		b1 = btst(index1);
		mG1=   mG(index1);
		stdG1= stdG(index1);

% standardize the data:
   for i=1:length(b1)-1
     b(i,1)=(b1(i)-mG1(i))/stdG1(i);
   end
   b(length(b1))=b1(length(b1));

% add weights:  
	b2=Wx(index1,index1)*b;

%% use either nnls.m or lsqnonneg.m depending on MatLab version
if str2num(vers)<6
    [x,dual] = nnls(G2(index1,:),b2);
else
	[x,resnorm] = lsqnonneg(G2(index1,:),b2);    
end

% calculate residuals for individual parameters
	err(index1,k)= G1(index1,:)*x-btst(index1)'; 

%add contributions from identically named water masses
	for i=1:n
		A(wm_index(i),k)    = A(wm_index(i),k) + x(i);
	end  

% in case of extended OMP analysis the biogeochemical part is 
% stored:
% NOTE: this has to be referenced with the appropriate ratio to
% convert into "mixage"
% default is changes in oxygen UNIT= ?mol/kg!!! and NOT years!!!
     switch(OMP)
      case 'ext'
       biogeo(k)=x(length(x))*(-ratio(3));
      end
	
	clear b
	end                   %end of loop with enough data

end  % %end of data point loop

%summary of run:
disp('  ')
disp('  ')
disp('  ')
disp('P R O G R A M   R U N   S U M M A R Y :')
disp('---------------------------------------')

switch(OMP)
case 'ext';
	disp('Method used:   EXTENDED OMP ANALYSIS.')
otherwise
	disp('Method used:   BASIC OMP ANALYSIS.')
end

disp(['Dataset used:  ' dataset '.'])
disp(['Selected data range:  ' selection])

disp('Parameters used:')
disp('  potential temperature')
disp('  salinity')
if esx(6)   == 1 disp('  oxygen'); end
if esx(7)   == 1 disp('  phosphate'); end
if esx(8)   == 1 disp('  nitrate'); end
if esx(9)   == 1 disp('  silicate'); end
if esx(10)   == 1 disp('  potential vorticity'); end
disp('  mass conservation');
disp('Weights used (variables as listed):')
disp(diag(Wx))
disp('  ')
disp('Water types used:')
disp('  ')
for i = 1:length(qwt_pos)
	disp(wmnames(5*(qwt_pos(i)-1)+1:5*(qwt_pos(i)-1)+5))
end
disp('  ')
disp('Water type definitions for the selected variables and mass conservation')
if OMP == 'ext' disp('Last column gives Redfield ratios'); end
disp('  ')
disp(G1)

disp(['successfully analysed datapoints:' num2str(100-100*gap/k) ' %' ])

disp('  ')
disp('Print this summary for reference and check that the results make sense.')
disp('Press any key to see a graph of the total residual')
disp('(mass conservation residual) against density.')
pause


% plotting residuals
figure
plot(100*err(m,:),pdens,'.','markers',10),axis('ij')
xlabel('mass conservation residual of fit (%)')

disp('  ')
j = 'n';
incontrol = input('Do you want to see more graphic output (y/n)?  [n]  ','s');
if length(incontrol) > 0 j = incontrol; end

if j == 'y'
%plotting water mass fractions
	for i = 1:nr_of_wm
		ctpara = i;
		tit_str = [tit_index(5*(i-1)+1:5*(i-1)+5)];
		figure
		contour2
		pause(2)
	end
% add a biogeochemistry plot if extended OMP 
switch(OMP)
case 'ext'
figure
  contour_bio
end
end

disp('  ')

%storing data in directory/folder OUTPUT
j = 'y';
incontrol = input('Do you want to store your results (y/n)?  [y]  ','s');
if length(incontrol) > 0 j = incontrol; end

if j == 'y'
	drswitch('Output is stored in');
		disp('  ')
		vname = 'result';
		incontrol = input('Give a file name for output storage: [result]  ','s');
		if length(incontrol) > 0 vname = incontrol; end
		incontrol = vname;
		lv = length(vname)+1;
switch(OMP)
case 'ext'	
disp('extended results written')
   vname = [vname '  nr_of_wm tit_index A err esx lat long press pdens biogeo'];
otherwise		
   vname = [vname '  nr_of_wm tit_index A err esx lat long press pdens'];
end
		
		if esx(4)  == 1, vname = [vname ' sal'];    end
		if esx(5)  == 1, vname = [vname ' ptemp'];  end
		if esx(6)  == 1, vname = [vname ' oxy'];    end
		if esx(7)  == 1, vname = [vname ' ph'];     end
		if esx(8)  == 1, vname = [vname ' ni'];     end
		if esx(9)  == 1, vname = [vname ' si'];     end
		if esx(10) == 1, vname = [vname ' pvort'];   end
		
		sout = sprintf('save %s',vname);
		eval(sout);
		disp('  ')
		disp(['File ' incontrol ' created and saved as: '  vname(1:lv) '.mat'])
		disp([' in:  ' pwd])
end

disp('  ')
disp('E N D   O F   O M P   A N A L Y S I S')
