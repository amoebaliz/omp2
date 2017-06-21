# function [G,mG,stdG]=norm_qwt(G1)
# Normalises the source water type (SWT) matrix G1
# and calculate standarddeviation and mean 
#
#  INPUT:
#  	G1	: Input nonnormalized SWt matrix
#
# OUTPUT:
#	G	: normalized SWT matrix
#	mG	: mean original SWT matrix
#	stdG	: standrddeviation of original SWT matrix
# 
# called by OMP_MAIN.M
#
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
def norm_qwt(G1):
    import numpy as np
    m,n = G1.shape # number of water types

    # mean and standarddeviation of SWT
    mG = np.mean(np.transpose(G1),axis=0)
    stdG =  np.std(np.transpose(G1),axis=0)
    G = np.zeros(G1.shape)
    # standardize QWT for m eq. 7
    for i in range(n):
        for kkk in range(m):
            G[kkk,i]=(G1[kkk,i]-mG[kkk])/stdG[kkk]
        G[m-1,i]=G1[m-1,i]  # mass untouched
    return G,mG,stdG
