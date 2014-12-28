import numpy as np

def WhiteCovar(x,bl,bl_long):

    #BASELINE LENGTH MUST BE EVEN TO SUBTRACT PAIRS
    nbaselines = np.int(len(x)/bl)
    nbl_long = np.int(len(x)/bl_long)

    if np.mod(bl,2) != 0:
        bl_safe = bl_long-1
    else:
        bl_safe = bl_long*1

    pairs = np.zeros(bl_safe/2)
    C_N = np.zeros(nbaselines,dtype='Float64')

    #FOR EACH BASELINE SUBTRACT INDEPENDENT PAIRS TO ESTIMATE WHITE NOISE
    for i in np.arange(nbl_long):

        for j in np.arange(bl_safe/2):
            pairs[j] = x[i*bl_safe + 2.*j] - x[i*bl_safe + 2.*j+1]


        psa = pairs.argsort()
        interquartile = np.std(pairs[psa[int(len(pairs)*0.25):int(len(pairs)*0.75)]])

        hivals = np.squeeze(np.where((np.abs(pairs) > interquartile*5.)))
        if hivals.size > 1:
            tp = pairs[np.squeeze(np.where((np.abs(pairs) < interquartile*5.)))]
        else:
            tp = pairs

        if i < nbl_long-1:
            C_N[i*bl_long/bl:(i+1)*bl_long/bl] = np.std(tp)/np.sqrt(2.)
        else:
            C_N[i*bl_long/bl:] = np.std(tp)/np.sqrt(2.)

    nans = np.squeeze(np.where(np.isnan(C_N) == 1))

    if nans.size > 0:
        notnans   = np.squeeze(np.where(np.isnan(C_N) == 0))
        C_N[nans] = np.mean(C_N[notnans])

    zeros = np.squeeze(np.where(C_N**2 == 0))

    if zeros.size > 0:
        notzeros   = np.squeeze(np.where(C_N**2 > 0))
        C_N[zeros] = np.mean(C_N[notzeros])

    return C_N**2
