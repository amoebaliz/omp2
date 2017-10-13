function [G1,wm,allsize]=qwt2(wm_row,ict)
%  QWT2.M: summary of water type definitions used by OMP2.M. 
%  Calling qwt2.m from the main OMP analysis program omp2.m will 
%  produce the source water matrix G1.
%
%  You will most likely have to edit this file for your own
%  application. We recommend that you save your edited file under
%  a different file name. See the web manual for details.
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

if nargin<1
 disp(' ')
 disp(['     Please give the rwo index of SWT definitions you want to use:'])
 disp([' e.g. qwt_step([1 2 3 6]) activates row 1, 2, 3, and 6 '])
 disp('  '),return
end

% The following lines identify the water masses defined through each water type.
% There is one water mass name for each water type definition row.
% Each name has to consist of exactly 5 letters including blanks; DO NOT CHANGE
% THE LENGTH unless you want to edit the main program omp2.m, too!
%
% NOTES:  1. You can store more than one water type definition for a water mass
%            and select which one you want to use at run time. This example file
%            demonstrates this by offering two definitions for ICW. Calling
%            qwt2([1 2 3 4]) activates rows 1, 2, 3, and 4 (AAMW and the first
%            ICW set); calling qwt2([1 2 5 6]) activates row 1, 2, 5, and 6
%            (AAMW and the second ICW set).
%         2. You can use this function to build up your own inventory of water
%            type definitions by including definitions for the Atlantic (NACW,
%            SACW, MedW etc.) and other oceans (AAIW, AABW etc.). This will
%            allow you to run OMP analysis for any ocean region by selecting the
%            appropriate rows in the function call.
%         3. Water types with identical names define water masses by property-
%            property relationships (e.g. in the thermocline); their contributions
%            will be added to all contributions of water masses with the same name.
%            For example, the call qwt2([1 2 3 4 7]) will produce three results,
%            AAMW (added from 1 and 2), ICW (added from 3 and 4) and AAIW.
%            For this to work all water types with identical names for which you
%            want contributions added in the result have to be called in an
%            uninterrupted sequence. (For example, calling qwt2([1 2 7]) will add
%            the AAMW contributions, calling qwt2([1 7 2]) will not.)
%

  wm(1:5)   = ' AAMW'; %                      first row
  wm(6:10)  = ' AAMW'; %                      second row
  wm(11:15) = '  ICW'; %upper, first set      third row
  wm(16:20) = '  ICW'; %lower, first set      forth row
  wm(21:25) = '  ICW'; %upper, second set     fifth row
  wm(26:30) = '  ICW'; %lower, second set     sixth row
  wm(31:35) = ' AAIW'; %                      seventh row
  wm(36:40) = '  IEW'; %                      eigth row

all=[...

%  The following lines define the water types. The order of parameters is
%  ptemp    sal      oxy    PO4     NO3    Si    mass   pvort
%  Note: potential vorticity is multiplied by 10*8.

  10       34.56      91    2.1    30     40     1.0    0.03;...  %1 lower AAMW
  16.4     34.55     100    1.4    19     25     1.0    1.12;...  %2 upper AAMW  
   9       34.65     260    1.1    15      5     1.0    0.03;...  %3 lower ICW, first set
  18       35.8      230    0       0      0.5   1.0    0.05;...  %4 upper ICW, first set
   9       34.72     209    1.47   20      5     1.0    0.03;...   %5 lower ICW, second set
  14.35    35.4      224    0.6     6.5    0.5   1.0    0.05;...  %6 upper ICW, second set
   4.5     34.35     210    2.2    32     35     1.0    0.30;...  %7 AAIW
   8.5     35         60    2.5    35     60     1.0    0.04;...  %8 IEW

];

G1=all(wm_row,:)';
allsize = size(all);
