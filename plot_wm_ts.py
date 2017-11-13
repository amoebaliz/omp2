import numpy as np
import matplotlib.pyplot as plt


stations = [26.7,28,30,35,40,45,50,55,60,70,80,90,100,110,120] #LINE 93.3
#stations = [51,55,60,70,80,90,100]                             #LINE 80.0
#stations = [51,55,60,70,80,90]                                 #LINE 66.7

data_fil = 'water_mass_fractions_75m.npz'
data_fil = 'water_mass_fractions_200m.npz'
dic_vals = np.load(data_fil)

time_vals = dic_vals['arr_1']
surf_wm_fracs = dic_vals['arr_0']

for sta in range(surf_wm_fracs.shape[2]):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for nwm in range(surf_wm_fracs.shape[1]):
        plt.plot(time_vals,surf_wm_fracs[:,nwm,sta],'o')

    ax.xaxis_date()
    ax.set_ylim(0,1)
    tit_str = str(stations[sta])
    plt.title(tit_str)
plt.show()
