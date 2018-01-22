import numpy as np
import datetime as dt
import matplotlib.dates as pltd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from math import ceil
import numpy as np
from scipy import linalg


def lowess(x, y, f=2. / 3., iter=3):
    # Lowess smoother: Robust locally weighted regression.
    # The lowess function fits a nonparametric regression curve to a scatterplot.
    # The arrays x and y contain an equal number of elements; each pair
    # (x[i], y[i]) defines a data point in the scatterplot. The function returns
    # the estimated (smooth) values of y.

    # The smoothing span is given by f. A larger value for f will result in a
    # smoother curve. The number of robustifying iterations is given by iter. The
    # function will run faster with a smaller number of iterations.
    n = len(x)
    print 'x',x
    r = int(ceil(f * n))
    h = [np.sort(np.abs(x - x[i]))[r] for i in range(n)]
    w = np.clip(np.abs((x[:, None] - x[None, :]) / h), 0.0, 1.0)  
    w = (1 - w ** 3) ** 3
    yest = np.zeros(n)
    delta = np.ones(n)
    for iteration in range(iter):
        for i in range(n):
            weights = delta * w[:, i]
            b = np.array([np.sum(weights * y), np.sum(weights * y * x)])
            A = np.array([[np.sum(weights), np.sum(weights * x)],\
                [np.sum(weights * x), np.sum(weights * x * x)]])
            print weights
            print y
            print 'b',b
            print 'A',A
            beta = linalg.solve(A, b)
            yest[i] = beta[0] + beta[1] * x[i]
        residuals = y - yest
        s = np.median(np.abs(residuals))
        delta = np.clip(residuals / (6.0 * s), -1, 1)
        delta = (1 - delta ** 2) ** 2
    
    return yest

#stations = [26.7,28,30,35,40,45,50,55,60,70,80,90,100,110,120] #LINE 93.3
#stations = [51,55,60,70,80,90,100]                             #LINE 80.0

years = mdates.YearLocator(2)

ln = 1
lns = ['93.0','80.0']
wms = ['PSA','PEW','NPCW']
stations  = [[26.7,28,30,35,40,45,50,55,60,70,80,90,100,110,120],\
             [51,55,60,70,80,90,100]]
sub_stats = [[30,60,90,120],[60,80,100]] #LINE 93.3; short

# LOAD DATA FOR PLOTTING
data_fil = 'water_mass_fractions_200m.npz'
dic_vals = np.load(data_fil)

time_vals = dic_vals['arr_1']
surf_wm_fracs = dic_vals['arr_0']

# FIGURE

f,ax = plt.subplots(len(sub_stats[ln]), 1,figsize=(4,10),sharex='col')
n = 0
for sta in range(surf_wm_fracs.shape[2]):
    if stations[ln][sta] == sub_stats[ln][n]:
       print 'MEEEP:',np.nanmax(np.sum(surf_wm_fracs[:,:,sta],axis=1)),np.nanmin(np.sum(surf_wm_fracs[:,:,sta],axis=1))

       for nwm in range(surf_wm_fracs.shape[1]):
           tms, = ax[n].plot(time_vals,surf_wm_fracs[:,nwm,sta],'o',label=wms[nwm])
            
       tit_str = 'Station ' + str(stations[ln][sta])
       ax[n].set_title(tit_str)
       ax[n].set_ylim(0,80)
       ax[n].set_yticks([0,20,40,60,80])
       
       if sub_stats[ln][n] == sub_stats[ln][-1]:
          break
       n+=1
    
ax[0].xaxis_date()
ax[0].set_xlim([pltd.date2num(dt.datetime(yr,1,1)) for yr in [1981,2017]])
ax[0].set_xticks([pltd.date2num(dt.datetime(yr,1,1)) for yr in range(1984,2017,10)])
ax[0].xaxis.set_minor_locator(years)

f.subplots_adjust(hspace=0.35)

# LEGEND
if ln == 0:
   plt.legend(loc='center',bbox_to_anchor=(0.5, 5.45),ncol=3)
elif ln == 1:
   plt.legend(loc='center',bbox_to_anchor=(0.5, 4.05),ncol=3)

# SAVE PNG FIGURE
fig_tit = lns[ln] + '_wt_frac.png'
plt.savefig(fig_tit)

# PLOTTING LINES
#col = ['cornflowerblue','orange','forestgreen']
#f=0.25
#print surf_wm_fracs.shape
#for sta in range(surf_wm_fracs.shape[2]):
#    print 'MEEP', sta
#    fig = plt.figure(figsize=(8, 4))
#    ax = fig.add_subplot(111)
#    for nwm in range(surf_wm_fracs.shape[1]):
#        wm_frac = surf_wm_fracs[:,nwm,sta].squeeze()
#        avg_frac = np.nanmean(wm_frac)
        #print wm_frac[11:]
        #yest = lowess(time_vals[11:], wm_frac[11:], f=f, iter=2)
#        P=ax.plot(time_vals,wm_frac,'o')
#        ax.plot([time_vals[0],time_vals[-1]],avg_frac*np.ones(2),\
 #               linestyle='--',color=col[nwm],linewidth=1,zorder=-1)
    
        #ax.plot(time_vals[11:],yest,color=col[nwm])
#    ax.xaxis_date()
#    ax.set_ylim(0,1)
#    ax.set_xlim(time_vals[4],time_vals[-1])
    #tit_str = str(stations[ln][sta])
#    plt.title(tit_str)
plt.show()
