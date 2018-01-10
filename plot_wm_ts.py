import numpy as np
import datetime as dt
import matplotlib.dates as pltd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

#stations = [26.7,28,30,35,40,45,50,55,60,70,80,90,100,110,120] #LINE 93.3
#stations = [51,55,60,70,80,90,100]                             #LINE 80.0

minorLocator = AutoMinorLocator()

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
       for nwm in range(surf_wm_fracs.shape[1]):
           tms, = ax[n].plot(time_vals,surf_wm_fracs[:,nwm,sta],'o',label=wms[nwm])
            
       tit_str = 'Station ' + str(stations[ln][sta])
       ax[n].set_title(tit_str)
       ax[n].set_ylim(0,1)
       ax[n].set_yticks([0,.5,1])
       ax[n].yaxis.set_minor_locator(minorLocator)
       if sub_stats[ln][n] == sub_stats[ln][-1]:
          break
       n+=1
ax[0].xaxis_date()
ax[0].set_xlim([pltd.date2num(dt.datetime(yr,1,1)) for yr in [1981,2017]])
ax[0].set_xticks([pltd.date2num(dt.datetime(yr,1,1)) for yr in range(1985,2017,10)])
ax[0].xaxis.set_minor_locator(minorLocator)

#handles, labels = ax[0].get_legend_handles_labels()
#labels = ['PSA','PEW','NPCW']
#plt.legend([handles,labels])

f.subplots_adjust(hspace=0.35)

#ax[-1].legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),
#          fancybox=True, shadow=True, ncol=3)


plt.legend(loc='center',bbox_to_anchor=(0.5, 5.45),ncol=3)
plt.legend(loc='center',bbox_to_anchor=(0.5, 4.05),ncol=3)

fig_tit = lns[ln] + '_wt_frac.png'
plt.savefig(fig_tit)

plt.show()
