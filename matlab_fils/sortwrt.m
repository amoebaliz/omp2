%  This routine is used by sortomp2.m to store individual
%  data types. It is mainly required because OMP Analysis
%  version 1 defined the names of variables with different
%  lengths.
%
%  sortwrt.m
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


if j < k
	if nn == 1
		ptsort(i-1,j:k) = interp1(press(statind(i-1)+1:statind(i)), ...
			lat(statind(i-1)+1:statind(i)),levels(j:k))';
	end
	
	if nn == 2
		ptsort(i-1,j:k) = interp1(press(statind(i-1)+1:statind(i)), ...
			long(statind(i-1)+1:statind(i)),levels(j:k))';
	end
	
	if nn == 3
		ptsort(i-1,j:k) = interp1(press(statind(i-1)+1:statind(i)), ...
			press(statind(i-1)+1:statind(i)),levels(j:k))';
	end
	
	if nn == 4
		ptsort(i-1,j:k) = interp1(press(statind(i-1)+1:statind(i)), ...
			sal(statind(i-1)+1:statind(i)),levels(j:k))';
	end
	
	if nn == 5
		ptsort(i-1,j:k) = interp1(press(statind(i-1)+1:statind(i)), ...
			ptemp(statind(i-1)+1:statind(i)),levels(j:k))';
	end
	
	if nn == 6
		ptsort(i-1,j:k) = interp1(press(statind(i-1)+1:statind(i)), ...
			oxy(statind(i-1)+1:statind(i)),levels(j:k))';
	end
	
	if nn == 7
		ptsort(i-1,j:k) = interp1(press(statind(i-1)+1:statind(i)), ...
			ph(statind(i-1)+1:statind(i)),levels(j:k))';
	end
	
	if nn == 8
		ptsort(i-1,j:k) = interp1(press(statind(i-1)+1:statind(i)), ...
			ni(statind(i-1)+1:statind(i)),levels(j:k))';
	end
	
	if nn == 9
		ptsort(i-1,j:k) = interp1(press(statind(i-1)+1:statind(i)), ...
			si(statind(i-1)+1:statind(i)),levels(j:k))';
	end
	
	if nn == 10
		ptsort(i-1,j:k) = interp1(press(statind(i-1)+1:statind(i)), ...
			pvort(statind(i-1)+1:statind(i)),levels(j:k))';
	end
end
