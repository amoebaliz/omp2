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
def omp2(OMP,nr_of_wm,tit_index,qwt_pos,wmnames,Wx,lat,switchpot,selection,lon,esx,press,sal,oxy,ptemp,pdens,ph,si,G1,wm_index):    
    from norm_qwt import norm_qwt
    import scipy
    import numpy as np
    import matplotlib.pyplot as plt
    from contour2 import contour2
    print '  '
    print 'OMP analysis now running. ', str(len(lat)) + ' data points found.'
    print '  '
    gap=0

    print 'Screening the data and reducing them to the selected range.'
    print '  '

    if switchpot == 'y':
	 eval('index=np.where((imag[pvort]==0) & (pvort<100) & (selection ))')
	 pvort = abs(pvort)
    else:
         print selection
         print np.min(pdens), np.max(pdens), np.min(press),np.max(press)
         exec('index = np.where(' + selection + ')[0]')

    lat   = lat[index]
    press = press[index]
    lon   = lon[index]
    sal   = sal[index]

    #lat   =  np.transpose(lat[index])
    #press =  np.transpose(press[index])
    #long  =  np.transpose(long[index])
    #sal   =  np.transpose(sal[index])
      
    if 'temp' in locals():
       temp = temp[index]
    if switchpot == 'y':
       pvort = pvort[index]
    if 'ptemp' in locals():
       ptemp = ptemp[index]
    else:
       ptemp = sw_ptmp(sal,temp,press,0)
    if 'pdens' in locals():
       pdens = pdens[index]
    else:
       pdens = sw_dens0(sal,temp) - 1000
    if esx[5] == 1: 
       oxy =  oxy[index]
    if esx[6] == 1: 
       ph  =  ph[index]
    if esx[7] == 1: 
       ni  =  ni[index]
    if esx[8] == 1: 
       si  =  si[index]
    del index

    print 'OMP analysis now running. ', str(len(lat)), ' data points to be analysed.'
    print '  '

    m,n = G1.shape[:] # n = number of water types, m = number of equations

    # normalise the source water matrix (get meanG, get stdG for weighting):
    G, mG, stdG = norm_qwt(G1)

    # EXTENDED OMP switch:
    # switch(OMP)
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
    G2=np.dot(Wx,G)

    gap=0
    # PYTH
    # ***********************************************************
    # This is the main loop for each data point; k = point index
    # First some initial settings
    err=np.zeros((m,len(lat)))-np.nan; 
    if OMP == 'ext':
       biogeo=np.zeros[1,len(lat)]-nan 
    A = np.zeros((wm_index[len(wm_index)-1],len(lat)))
    # Vector of each datapoint (btst) is build here
    for k in range(len(lat)):
	# selecting the correct parameters
        btst = np.append(ptemp[k],sal[k])
        if esx[5] == 1: 
           btst = np.append(btst,oxy[k])
        if esx[6] == 1:  
           btst = np.append(btst,ph[k])
        if esx[7] == 1:  
           btst = np.append(btst,ni[k])
        if esx[8] == 1:
           btst = np.append(btst,si[k])  
        if 'pdens'in locals():
           pden_dat = pdens[k]
        if esx[9] == 1:
           btst = np.append(btst,pvort[k])

        btst = np.append(btst,1)

        index1=np.where(~np.isnan(btst))[0]
        index0=np.where(np.isnan(btst))[0]
        cutit=n
        # (because we have one unknown more!)
        if OMP[:3]=='ext':
           cutit=n+1

        if len(index1) < cutit+1: # if 
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
        b = np.zeros(len(b1))
        for i in range(len(b1)-1):
            b[i]=(b1[i]-mG1[i])/stdG1[i]
        b[len(b1)-1]=b1[len(b1)-1]

        # add weights:  
        b2=Wx[index1,index1]*b

        ## use either nnls.m or lsqnonneg.m depending on MatLab version
        x, rnorm = scipy.optimize.nnls(G2[index1,:],b2)
        # calculate residuals for individual parameters
        err[index1,k] = np.dot(G1[index1,:],x) - np.transpose(btst[index1]) 
        #add contributions from identically named water masses
        for i in range(n):
	    A[wm_index[i]-1,k]    = A[wm_index[i]-1,k] + x[i]

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
    for i in qwt_pos:
	print wmnames[i]
    print '  '
    print 'Water type definitions for the selected variables and mass conservation'
    if OMP == 'ext': print 'Last column gives Redfield ratios'
    print '  '
    print G1

    print 'successfully analysed datapoints:', str(100-100*gap/(k+1)), ' %' 

    print '  '
    print 'Print this summary for reference and check that the results make sense.'
    print 'Press any key to see a graph of the total residual'
    print '(mass conservation residual) against density.'


    # plotting residuals
    fig, ax = plt.subplots()
    ax.plot(100*err[m-1,:],pdens,'o')
    #ax.set_ylim([22,28])
    #ax.set_xlim([-50,150])
    ax.invert_yaxis()
    ax.set_xlabel('mass conservation residual of fit (%)')
    ax.set_ylabel('density')
    print '  '
    j = 'n'
    incontrol = input('Do you want to see more graphic output (y/n)?  [n]  ')
    if len(incontrol) > 0: 
       j = incontrol

    if j == 'y':
       # plotting water mass fractions
       for i in range (nr_of_wm):
	   ctpara = i
	   tit_str = tit_index[i]
	   contour2(ctpara, tit_str, A, lat, lon,press)
    # add a biogeochemistry plot if extended OMP 
    if OMP == 'ext':
       plt.figure()
       contour_bio
    print tit_index
    print '  '

    # storing data in directory/folder OUTPUT
    j = 'y'
    #incontrol = input('Do you want to store your results (y/n)?  [y]  ')
    #if len(incontrol) > 0:
    #   j = incontrol

    #if j == 'y':
    #   drswitch('Output is stored in')
    #   print '  '
    #   vname = 'result'
    #   incontrol = input('Give a file name for output storage: [result]  ','s')
    #   if len(incontrol) > 0:
    #       vname = incontrol
    #   incontrol = vname
    #   lv = len(vname)+1

    #if OMP == 'ext':
    #   print 'extended results written'
    #   vname = vname + '  nr_of_wm tit_index A err esx lat long press pdens biogeo'
    #else: 
    #   vname = vname + '  nr_of_wm tit_index A err esx lat long press pdens'
		
    #if esx[3]  == 1: vname = vname + ' sal'
    #if esx[4]  == 1: vname = vname + ' ptemp'
    #if esx[5]  == 1: vname = vname + ' oxy'
    #if esx[6]  == 1: vname = vname + ' ph'
    #if esx[7]  == 1: vname = vname + ' ni'
    #if esx[8]  == 1: vname = vname + ' si'
    #if esx[9] == 1: vname = vname + ' pvort'
		
    # SAVE SOMETHING 
    #eval('save %s vname')
    #print '  '
    #print 'File ', incontrol, ' created and saved as: ',  vname[:lv], '.mat'
    #print ' in:  ', pwd

    print '  '
    print 'E N D   O F   O M P   A N A L Y S I S'
    plt.show()
