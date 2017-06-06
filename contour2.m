%  This routine contours water mass content calculated
%                    by omp2.m
%
%  NOTES: 1. This routine is called from omp2.m on request.
%         2. To use contourp2.m as a separate call from the
%            command window you have to run omp2.m first
%            and keep all variables in the workspace.
%
%---------------------------------------------
% CALL: sw_dist.m form CSIRO seawater package
% external variables required: ctpara tit_str A lat long press
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

% select parameter to be plotted:
para=A(ctpara,:)*100;

%calculate distance between stations
[dist,phaseangle] = sw_dist(lat,long,'km');

if size(dist,1)>1
 dist=dist';
end

cumdist=[0 cumsum(dist)];

%check for duplicates and separate them by 0.5 m.
for i=2:length(press)
	if cumdist(i-1) == cumdist(i) & press(i-1) == press(i)
		press(i-1) = press(i-1) - 0.5;
	end
end


% create regular grid:
XI=linspace(min(cumdist),max(cumdist),20)';
YI=linspace(min(press),max(press),20);

% remove NaN do apply griddata.m
press1=press(~isnan(para));
cumdist1=cumdist(~isnan(para));
para1=para(~isnan(para));


% interpolate to regular grid:
para2=griddata(cumdist1,press1,para1,XI,YI,'v4');

contour(linspace(min(cumdist1),max(cumdist1),20), ...
                 linspace(min(press1),max(press1),20),para2,[0:10:100])
shading flat
hold on

[C,h]=contour(linspace(min(cumdist1),max(cumdist1),20), ...
                 linspace(min(press),max(press),20),para2,[0:10:100]);
clabel(C,[0 20 40 60 80])
colormap gray
axis ij
plot(cumdist,press,'k.','markersize',15)
caxis([-50 100])
set(gca,'position',[.1 .1 .85 .5])

xlabel('distance (km)')
ylabel('pressure (dbar)')
title([tit_str ' water mass content (percent)'])
set(gca,'box','on')

