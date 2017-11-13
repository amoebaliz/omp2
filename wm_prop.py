# This program is designed to calculate the fraction
# of water mass over specified depth range
# ---------------------------------------------
# CALL: sw_dist.m form CSIRO seawater package
# external variables required: ctpara tit_str A lat long press
# ---------------------------------------------

def wm_prop(para, stations, stats, lat, lon, press):
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.interpolate import griddata
    from sw_dist import sw_dist 


    # remove NaN do apply griddata.m
    # press1=press[~np.isnan(para)].squeeze()
    # cumdist1=cumdist[~np.isnan(para)].squeeze()
    para1=para[~np.isnan(para)]

    # DEPTH RANGES: above 75 for surface; 100m-200m for sub-surface
    slim = 200.
    mld1 = 100.
    mld2 = 200.

    # INITIALIZE WATER FRACTION ARRAY: allocate space for each station
    # and assign to nan in case particular station is not taken during 
    # a given survey
    surf_frac = np.empty(len(stations))
    surf_frac[:] = np.NAN

    # IDENTIFY THE BEGINNING OF EACH CAST BASED ON DEPTH
    # (press==0 @ surface)
    Istart = np.where(press==0)[0]

    # ITTERATE OVER STATIONS
    for ns in range(len(stations)):
        # CHECK FOR MULTIPLE CASTS AT THE SAME STATION
        ncast = np.sum(stats[Istart]==stations[ns])
        # INDEX THE SURFACE MEASUREMENTS FOR THESE CASTS
        IC = np.where(stats[Istart]==stations[ns])[0]
        # ALLOCATE SPACE FOR AVERAGING ACROSS STATION CASTS 
        tmp_str = np.zeros(ncast)

        # ITTERATE OVER CASTS AT A GIVEN STATION:
        # Most times, ncast = 1
        for nc in range(ncast): 

            if Istart[IC[nc]] == Istart[-1]:
               # IF ON LAST STATION ON THE LINE
               sta_press = press[Istart[IC[nc]]:]
               sta_param = para1[Istart[IC[nc]]:]
            else: 
               sta_press = press[Istart[IC[nc]]:Istart[IC[nc]+1]]
               sta_param = para1[Istart[IC[nc]]:Istart[IC[nc]+1]]

            # SURFACE RANGE INDICES   
            Iys = np.where(sta_press<=slim)[0] 
            # SUBSURFACE RANGE
            Iym = np.where((sta_press>=mld1) & (sta_press<=mld2))[0]

            # GET THE DEPTH DIST. BETWEEN MEASUREMENTS FOR DIFFERENT STATION SCENARIOS
            
            # The last indexed value will be less than the total number of depth measurements
            # if the entire cast does not fall in the desired range: Most likely at near-shore stations 
            # SURFACE
            if Iys[-1]<(len(sta_press)-1): 
               # surface differencing, including 1 idepth measurement below the bottom of range  
               sdep_dif = np.diff(sta_press[np.append(Iys,Iys[-1]+1)])
                   
               # If the halfway point between last depth measurement in range and first
               # measurement out of range is below the depth limit, only need influence 
               # of last measurement in range. Change value last diff value to distance 
               # between bottom depth and bottom limit
               if (sta_press[Iys[-1]+1]+sta_press[Iys[-1]])/2>slim:
                  sdep_dif[-1] = slim-sta_press[Iys[-1]]
               # ELSE - include some influence below last in-range depth measurement
               else:
                  sdep_dif = np.append(sdep_dif,slim-(sta_press[Iys[-1]]+sdep_dif[-1]/2)) 
            
            else: #LAST STATION MEASUREMENT IS IN (SURFACE) RANGE
               sdep_dif = np.diff(sta_press[Iys])
               sdep_dif = np.append(sdep_dif,np.min((sdep_dif[-1],slim-sta_press[-1])))
            
            # calculate total depth ranges for the surface later
            sdep_wts = (sdep_dif + np.append(0,sdep_dif[:-1]))/2.
            Iwm = np.where(sta_param[Iys]>50)[0]
            tmp_str[nc] = np.sum(sdep_wts[Iwm])/np.sum(sdep_wts)
        # AVERAGE over casts for a given station
        surf_frac[ns] = np.mean(tmp_str)
    
    return surf_frac
