# OMP2.M: OMP analysis main program version 2
#     
# This is an updated version of an easy-to-handle package for the use of 
# OMP analysis to resolve fractions of water masses involved in the
# mixing of water masses at a given point in the ocean. The original
# version was prepared by Johannes Karstensen. This version incorporates
# improvements by Matthias Tomczak.
#
# This program is called by omp2int.m, omp2gui.m and omp2auto.m and will not
# run before one of these programs is called and placed all necessary
# variables into the workspace.
#
#
# Function calls used: qwt2.m qwt_tst.m nansum.m (Philip Morgan, CSIRO)
#  sw_ptmp sw_dens0.m (Philip Morgan, CSIRO) may be called for some data files
#  sw_dist.m (Philip Morgan, CSIRO) is called through the contour2 call
#
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

# 
def omp2():

    print '  '
    print 'OMP analysis now running. ', str(len(lat)) + ' data points found.'
    print '  '

    #starttime = clock;
    gap=0

    # Cut the data to the selcted range 
    # 
    # set potential vorticity to positive values (independent of hemisphere) if required

    print 'Screening the data and reducing them to the selected range.'
    print '  '

    switch(switchpot)
    if switchpot == 'y':
	 eval('index=np.where((imag[pvort]==0) & (pvort<100) & (selection ))')
	 pvort = abs(pvort)
    else:
	 eval('index=np.where(selection)[0]')

    lat   =  np.transpose(lat[index])
    press =  np.transpose(press[index])
    long  =  np.transpose(long[index])
    sal   =  np.transpose(sal[index])
      
    if 'temp' in locals():
       temp = np.transpose(temp[index])
    if ~isempty(switchpot) & switchpot == 'y':
       pvort = np.transpose(pvort(index))
    if 'ptemp' in locals():
       ptemp = np.transpose(ptemp[index])
    else:
       ptemp = sw_ptmp(sal,temp,press,0)
    if 'pdens' in locals():
       pdens = np.transpose(pdens[index])
    else:
       pdens = sw_dens0(sal,temp) - 1000
    if esx[6] == 1: 
       oxy =  np.transpose(oxy[index])
    if esx[7] == 1: 
       ph  =  np.transpose(ph[index])
    if esx[8] == 1: 
       ni  =  np.transpose(ni[index])
    if esx[9] == 1: 
       si  =  np.transpose(si[index])
    del index

    print 'OMP analysis now running. ', num2str(length(lat)), ' data points to be analysed.'
    print '  '

    m,n = G1.shape[:] # n = number of water types, m = number of equations

    # normalise the source water matrix (get meanG, get stdG for weighting):
    G, mG, stdG = norm_qwt(G1)

    # EXTENDED OMP switch:
    switch(OMP)
    if OMP == 'ext':
         # Adding Redfield ratio to the system, ratio comes from weight file
         G1[:m,n]=np.transpose(redfrat[:m])
	 # normalisation of the ratios:
	 # ---------------------------------
 	 for rr in range(m-1):
  	     G[rr,n] = redfrat[rr]*(max(G[rr,:n])-min(G[rr,:n])) \
           /(max(G1[rr,:n])-min(G1[rr,:n]))
	 G[m,n] = 0

    # adding weights
    G2=Wx*G

    gap=0
    # PYTH
    # ***********************************************************
    # This is the main loop for each data point; k = point index
    # First some initial settings
    err=zeros((m,length(lat)))-nan; 

    if OMP == 'ext':
       biogeo=np.zeros[1,len(lat)]-nan 

    A[:wm_index[len(wm_index)-1],:len(lat)] = \
                 np.zeros((wm_index[len(wm_index)-1],length(lat)));

    oxy_dat = []
    ph_dat  = []
    ni_dat  = []
    si_dat  = []
    pv_dat  = []

    # Vector of each datapoint (btst) is build here
    for k in len(lat):
	# selecting the correct parameters
        p_dat	= press[k]
        t_dat	= ptemp[k] 
        s_dat	=   sal[k] 
    if esx[5] == 1: 
       oxy_dat	= oxy[k]
    if esx[6] == 1:  
       ph_dat	= ph[k]
    if esx[7] == 1:  
       ni_dat	= ni[k]
    if esx[8] == 1:  
       si_dat  = si[k]
    if 'pdens'in locals():
       pden_dat = pdens[k]
    if esx[9] == 1:
       pv_dat = pvort[k]
    kon=1
	
    btst= [t_dat,s_dat,oxy_dat,ph_dat,ni_dat,si_dat,pv_dat,kon]

    if etime(clock,starttime) > 5:
       print num2str[k], ' data points analysed so far.'
       starttime = clock
    # looking for GAP (indicated through NaN):
    index1=find(~isnan(btst))
    index0=find(isnan(btst))

    cutit=n;
    # using extended OMP we need one parameter more 
    # (because we have one unknown more!)
    if OMP[:3]=='ext':
       cutit=n+1

    if length(index1) < cutit+1 %if1 :
       # not enough parameters to find a NNLS fit
       # DATA point not successful analysed
       print 'ANALYSIS of the datapoint failed, not enough parameters available !!'
       A[:nr_of_wm,k] = nan
       Dual[:nr_of_wm,k] = nan
       gap+=1

    else:
       #new data without GAP:
       b1    = btst[index1]
       mG1   =   mG[index1]
       stdG1 = stdG[index1]

    # standardize the data:
    for i in range(length(b1)-1):
        b[i,0]=(b1[i]-mG1[i])/stdG1[i]
    b[len(b1)-1]=b1[len(b1)-1]

    # add weights:  
    b2=Wx[index1,index1]*b

    ## use either nnls.m or lsqnonneg.m depending on MatLab version
    if str2num(vers)<6:
       x,dual = nnls[G2[index1,:],b2]
    else:
       x,resnorm = lsqnonneg[G2[index1,:],b2]    

    # calculate residuals for individual parameters
    err[index1,k] = np.transpose(G1[index1,:]*x-btst[index1]) 

    #add contributions from identically named water masses
    for i in range(n):
	A[wm_index[i],k]    = A[wm_index[i],k] + x[i]

    # in case of extended OMP analysis the biogeochemical part is 
    # stored:
    # NOTE: this has to be referenced with the appropriate ratio to
    # convert into "mixage"
    # default is changes in oxygen UNIT= ?mol/kg!!! and NOT years!!!
    if OMP == 'ext':
       biogeo[k]=x[len(x)-1]*(-ratio[3])
    del b
        #end of loop with enough data

       ## end of data point loop

    #summary of run:
    print '  '
    print '  '
    print '  '
    print 'P R O G R A M   R U N   S U M M A R Y :'
    print '---------------------------------------'

    if OMP == 'ext':
       print 'Method used:   EXTENDED OMP ANALYSIS.'
    else:
       print 'Method used:   BASIC OMP ANALYSIS.'

       print 'Dataset used:   dataset .'
       print 'Selected data range:   selection'

    print 'Parameters used:'
    print '  potential temperature'
    print '  salinity'
    if esx[5]   == 1: print '  oxygen'
    if esx[6]   == 1: print '  phosphate'
    if esx[7]   == 1: print '  nitrate'
    if esx[8]   == 1: print '  silicate'
    if esx[9]   == 1: print '  potential vorticity'
    print '  mass conservation'
    print 'Weights used (variables as listed):'
    print np.diag(Wx)
    print '  '
    print 'Water types used:'
    print '  '
    for i in range(len(qwt_pos)):
	print wmnames[5*(qwt_pos[i]-1):5*(qwt_pos[i]-1)+5]
    print '  '
    print 'Water type definitions for the selected variables and mass conservation'
    if OMP == 'ext': print 'Last column gives Redfield ratios'
    print '  '
    print G1

    print 'successfully analysed datapoints:', num2str(100-100*gap/k), ' %' 

    print '  '
    print 'Print this summary for reference and check that the results make sense.'
    print 'Press any key to see a graph of the total residual'
    print '(mass conservation residual) against density.'


    # plotting residuals
    plt.figure()
    plt.plot(100*err[m,:],pdens,'.','markers',10),axis('ij')
    xlabel('mass conservation residual of fit (%)')

    print '  '
    j = 'n'
    incontrol = input('Do you want to see more graphic output (y/n)?  [n]  ','s')
    if length(incontrol) > 0: j = incontrol

    if j == 'y':
       # plotting water mass fractions
       for i in range (nr_of_wm):
	   ctpara = i
	   tit_str = tit_index[5*(i-1):5*(i-1)+5]
	   np.figure()
	   contour2
    # add a biogeochemistry plot if extended OMP 
    if OMP == 'ext':
       plt.figure()
       contour_bio

    print '  '

    # storing data in directory/folder OUTPUT
    j = 'y'
    incontrol = input('Do you want to store your results (y/n)?  [y]  ','s')
    if length(incontrol) > 0:
       j = incontrol

    if j == 'y':
       drswitch('Output is stored in')
       print '  '
       vname = 'result'
       incontrol = input('Give a file name for output storage: [result]  ','s')
       if length(incontrol) > 0:
           vname = incontrol
       incontrol = vname
       lv = length(vname)+1

    if OMP == 'ext':
       print 'extended results written'
       vname = vname + '  nr_of_wm tit_index A err esx lat long press pdens biogeo'
    else: 
       vname = vname + '  nr_of_wm tit_index A err esx lat long press pdens'
		
    if esx[3]  == 1: vname = vname + ' sal'
    if esx[4]  == 1: vname = vname + ' ptemp'
    if esx[5]  == 1: vname = vname + ' oxy'
    if esx[6]  == 1: vname = vname + ' ph'
    if esx[7]  == 1: vname = vname + ' ni'
    if esx[8]  == 1: vname = vname + ' si'
    if esx[9] == 1: vname = vname + ' pvort'
		
    # SAVE SOMETHING 
    eval('save %s vname')
    print '  '
    print 'File ', incontrol, ' created and saved as: ',  vname[:lv], '.mat'
    print ' in:  ', pwd

    print '  '
    print 'E N D   O F   O M P   A N A L Y S I S'
