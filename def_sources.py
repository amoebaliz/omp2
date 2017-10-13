# How to define the source water masses ? Just an example
#
# EXAMPLE: data from throughflow region (Timor Sea)
#
# -------------------------------------------
# This program is part of the OMP package from:
# GEOMAR
# Helmholtz Centre for Ocean Res. Kiel  FIAMS, Flinders University
# J. Karstensen                         Matthias Tomczak
# Duesternbrooker Weg 20				GPO Box 2100
# 24106 Kiel                            Adelaide, SA
# Germany                               Australia
#
# BUGS: jkarstensen@geomar.de
#   or  matthias.tomczak@flinders.edu.au
# -------------------------------------------
#
# OUTPUT to be used in OMP packet:

# -------------------------------------------
# G:  	value to fill in qwt_step.m
# err:	value to caluculate weigths

import scipy
import numpy as np
import matplotlib.pyplot as plt

# -------------------------------------------
# |             EDITING SECTION             | 
# -------------------------------------------

# load dataset from source region of a water mass 
# (*.MAT file in standard format):
# load source_data

# PYTHON READ IN
#CCS_sw_dict = np.load()

# flags: insert here what kind of flag is used (eg. -9, NaN, 99999, etc.) 
flag=np.nan

# name of "independent" parameter (as string variable ' '!!):
nam_indep='ptemp'

# choose your own fit range:
fit_max=16.4
fit_min=10

# names of "dependent" parameter to be fitted (as string variable ' ';' ')
# para=[' sal',' oxy','  ph','  ni']

# ------------------------------------------------------
# |                    FROM HERE ON:                   | 
# |   EDITING ONLY IF YOU KNOW WHAT YOU ARE DOING !!   |
# ------------------------------------------------------

eval(['in_depend=' nam_indep ';'])

# sort by independent parameter
[in_depend,I]=sort(in_depend)

# number of independent parameter
no_para=len(para)

# find independent parameter range:
ind = np.where((in_depend<=fit_max)&&(in_depend>=fit_min))[0]

# check for nan  
check=np.isnan(in_depend)

if ~isempty(check):

   print ' '
   print ' Flag check: ', nam_indep
   print ' Remove NaN from data for use in OMP analysis!'



# MIN/MAX indices of indenpendent parameter
ind_low=np.min(ind)
ind_upp=np.max(ind)

G[0,:2]=[fit_max, fit_min]

# split subplots:
sub1=ceil((no_para)/2.)

print ' '

# -------------------------------------------
# |       loop to fit all parameters        |
# -------------------------------------------

for main in range(no_para):  #  START OF MAIN LOOP

# sort all variables accordingly
eval([para(main,:) '=' para(main,:) '(I);'])


disp([' Now fitting: ' para(main,:)])
 eval(['parafit=' para(main,:) ';'])

     subplot(sub1,2,main) 
   plot(parafit,in_depend,'c.','markersize',15) 
   xlabel(para(main,:)) 
   ylabel(nam_indep)

# check for nan  
check=find(isnan(parafit))
if ~isempty(check):
 print ' Remove NaN from data for use in OMP analysis!')
 print ' '
 

#
ix=find(in_depend(ind)~=flag&~isnan(in_depend(ind))&parafit(ind)~=flag&~isnan(parafit(ind)))

# fit to data and error estimate    
[coeff,er]=polyfit(in_depend(ind(ix)),parafit(ind(ix)),1)  
[da,err_fit]=polyval(coeff,in_depend(ind(ix)),er)
hold on,plot(da,in_depend(ind(ix)),'r-','linewidth',2)

# sample fit to SWT matrix:
G(main+1,1:2)=polyval(coeff,[fit_max fit_min])
err(main)=np.mean(err_fit)
# nicer axis of plot
  set(gca,'ylim',[fit_min-2 fit_max+2])

# END OF MAIN LOOP

print ' '
print 'Fitted Source water matrix (G) is:'
print G

print '  '
print 'Mean fit error (err) of dependent variables: '
print err

