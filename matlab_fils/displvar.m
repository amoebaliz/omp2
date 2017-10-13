%  DISPLVAR.M: displays individual stations from omp2 produced by statns2.m.
%
%  displvar.m
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

disp('  ')
disp('available variables:')
who

%define input file
incontrol = input('Which variable for display?  [salinity]  ','s');
if length(incontrol) > 0
	dataset = incontrol;
else
	dataset = 'salinity';
end

incontrol = input('Overplotting (y/n)?  [n]  ','s');
if incontrol == 'y'
	oplt = incontrol;
else
	oplt = 'n';
end

figure
hold off
if oplt == 'y' hold on; end
dept = -levels;
for i = 1:length(salinity(1,:))
	s = ['plot(' dataset '(:,i),dept)'];
	eval(s);
	pause
end
