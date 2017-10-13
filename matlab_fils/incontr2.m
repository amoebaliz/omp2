%  This is the control file for automatic mode of OMP analysis. 
%  Calling incontr2.m from the main OMP analysis program omp2auto.m 
%  will set all parameters needed to run omp2.m.
%
%  You will have to edit this file for your own application.
%  We recommend that you save your edited file under
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

% Set type of OMP analysis by placing/deleting comment (%) sign:
 OMP = 'cla'; % classical OMP analysis
% OMP = 'ext'; % extended OMP analysis

% Define the file which contains the data by replacing the word testdata:
dataset = 'testdata';

%%% define a specific range of data to be used:
%   as a string variable
%
%** for example: to use only data between potential 
%   density 23.12 and 28.03 and oxygen larger then 20 please write 
%   selection='pdens>=23&pdens<=28&oxy>=20';
selection='pdens>=23&pdens<=28&oxy>=20';


% Select/deselect potential vorticity by setting switchpot to 'y' or 'n':
switchpot = 'n';

% Select/deselect variables by setting corresponding switches to 'y' or 'n':
iox = 'y'; % oxygen switch
iph = 'y'; %phosphate switch
ini = 'y'; % nitrate switch
isi = 'n'; %silicate switch

% Define the file which contains the weights by replacing the word testwght: 
weightset = 'testwght';

% Define the routine which contains the source water type definitions:
swtypes = 'qwt2';

% Set the total number of water masses to be included in the analysis
wm = 2;

%  Select the water type numbers (row in the water type matrix)
% (The square brackets have to contain wm numbers):

qwt_pos = [1 2 3 4];




