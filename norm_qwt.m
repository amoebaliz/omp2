function [G,mG,stdG]=norm_qwt(G1)
% Normalises the source water type (SWT) matrix G1
% and calculate standarddeviation and mean 
%
%  INPUT:
%  	G1	: Input nonnormalized SWt matrix
%
% OUTPUT:
%	G	: normalized SWT matrix
%	mG	: mean original SWT matrix
%	stdG	: standrddeviation of original SWT matrix
% 
% called by OMP_MAIN.M
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

[m,n]=size(G1); % number of water types

% mean and standarddeviation of SWT
  mG = mean(G1');
stdG =  std(G1');

% standardize QWT for m eq. 7
for i=1:n
  for kkk=1:m-1
   G(kkk,i)=(G1(kkk,i)-mG(kkk))/stdG(kkk);
  end  
   G(m,i)=G1(m,i);  % mass untouched
end
