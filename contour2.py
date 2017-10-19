#  This routine contours water mass content calculated
#                    by omp2.m
#
#  NOTES: 1. This routine is called from omp2.m on request.
#         2. To use contourp2.m as a separate call from the
#            command window you have to run omp2.m first
#            and keep all variables in the workspace.
#
# ---------------------------------------------
# CALL: sw_dist.m form CSIRO seawater package
# external variables required: ctpara tit_str A lat long press
#
# ---------------------------------------------
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
# --------------------------------------------

def contour2(ctpara, tit_str, A, lat, lon, press):
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.interpolate import griddata
    from sw_dist import sw_dist 
    # select parameter to be plotted:
    para=A[ctpara,:]*100

    # calculate distance between stations
    dist = sw_dist(lat,lon,'km')

    #if dist.size[0]>1:
    #   dist=dist'
    #end
    cumdist=np.append(0, np.cumsum(dist))

    #check for duplicates and separate them by 0.5 m.
    for i in range(1,len(press)):
	if ((cumdist[i-1] == cumdist[i]) & (press[i-1] == press[i])):
           press[i-1] = press[i-1] - 0.5
	# end
    #end

    # create regular grid:
    XI=np.linspace(min(cumdist),max(cumdist),20)
    YI=np.linspace(min(press),max(press),20)
    X,Y = np.meshgrid(XI, YI)  
    # remove NaN do apply griddata.m
    press1=press[~np.isnan(para)].squeeze()
    cumdist1=cumdist[~np.isnan(para)].squeeze()
    para1=para[~np.isnan(para)]
    # interpolate to regular grid:
    points = np.column_stack((cumdist1, press1))
    para2=griddata(points,para1,(X,Y),method='cubic')
    levs = [0,10,20,30,40,50,60,70,80,90,100]   
 
    fig, ax = plt.subplots()
    ax.contour(np.linspace(min(cumdist1),max(cumdist1),20), \
                 np.linspace(min(press1),max(press1),20),para2,levs,shading='flat') #,[0:10:100])

    C = ax.contour(np.linspace(min(cumdist1),max(cumdist1),20), \
                 np.linspace(min(press),max(press),20),para2,levs,shading='flat') # ,[0:10:100]);

    plt.clabel(C, inline=1, fontsize=10)
    #clabel(C,[0 20 40 60 80])
    #colormap gray
    #axis ij
    plt.plot(cumdist,press,'ko')
    #caxis([-50 100])
    #set(gca,'position',[.1 .1 .85 .5])
    #ax.invert_yaxis()
    ax.set_ylim([300,500])
    ax.invert_yaxis()
    ax.set_xlabel('distance (km)') 
    ax.set_ylabel('pressure (dbar)')
    tit_text = tit_str + ' water mass content (percent)'
    plt.title(tit_text)
    #set(gca,'box','on')
