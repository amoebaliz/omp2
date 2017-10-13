% How to define the source water masses ? Just an example
%
% EXAMPLE: data from throughflow region (Timor Sea)
%
%
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
%
% OUTPUT to be used in OMP packet:
% *******************************************
% G:  	value to fill in qwt_step.m
% err:	value to caluculate weigths
%
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% EDITING SECTION %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% input dataset from source region of a water mass 
%%%%% (*.MAT file in standard format):
load source_data

%%%%%%%%%% flags: insert here what kind of flag is used (eg. -9, NaN, 99999, etc.) 
flag=nan;

%%%%%%%% name of "independent" parameter (as string variable ' '!!):
nam_indep='ptemp';

%%%%%%%% choose your own fit range:
fit_max=16.4; %
fit_min=10;  %

%%%%%%%% names of "dependent" parameter to be fitted (as string variable ' ';' ')
para=[' sal';' oxy';'  ph';'  ni'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%  FROM HERE ON:  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% EDITING ONLY IF YOU KNOW WHAT YOU ARE DOING !! %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eval(['in_depend=' nam_indep ';'])

%%%%%%%% sort by independent parameter
[in_depend,I]=sort(in_depend);

%%%%%%% number of independent parameter
no_para=size(para,1);

%%%%%%% find independent parameter range:
ind=find(in_depend<=fit_max&in_depend>=fit_min);

%%%%% check for nan  
check=find(isnan(in_depend));
if ~isempty(check)
 disp(' ')
 disp([' Flag check: ' nam_indep])
 disp(' Remove NaN from data for use in OMP analysis!')

end 


%%%%%%%%% MIN/MAX indices of indenpendent parameter
ind_low=min(ind);
ind_upp=max(ind);

G(1,1:2)=[fit_max fit_min];

%%%%% split subplots:
sub1=ceil((no_para)/2);

disp(' ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% loop to fit all parameters  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for main=1:no_para  %%% START OF MAIN LOOP

%%%%%%% sort all variables accordingly
eval([para(main,:) '=' para(main,:) '(I);'])


disp([' Now fitting: ' para(main,:)])
 eval(['parafit=' para(main,:) ';'])

     subplot(sub1,2,main) 
   plot(parafit,in_depend,'c.','markersize',15) 
   xlabel(para(main,:)) 
   ylabel(nam_indep)

%%%%% check for nan  
check=find(isnan(parafit));
if ~isempty(check)
 disp(' Remove NaN from data for use in OMP analysis!')
 disp(' ')
end 
 

%%%%%%%%%%%%%%%%%
ix=find(in_depend(ind)~=flag&~isnan(in_depend(ind))&parafit(ind)~=flag&~isnan(parafit(ind)));

%%%%%%%%  fit to data and error estimate  %%%%%%%%%%   
[coeff,er]=polyfit(in_depend(ind(ix)),parafit(ind(ix)),1);   
[da,err_fit]=polyval(coeff,in_depend(ind(ix)),er);
hold on,plot(da,in_depend(ind(ix)),'r-','linewidth',2)

% sample fit to SWT matrix:
G(main+1,1:2)=polyval(coeff,[fit_max fit_min]);
err(main)=mean(err_fit);
%%% nicer axis of plot
  set(gca,'ylim',[fit_min-2 fit_max+2])

end  %%%% OF MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp('Fitted Source water matrix (G) is:')
G

disp('  ')
disp('Mean fit error (err) of dependent variables: ')
err