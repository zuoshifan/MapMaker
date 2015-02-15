import numpy as np
from mpi4py import MPI
import MPI_tools



def WhiteCovar(x,bl,bl_long,comm=MPI.COMM_WORLD):

    #BASELINE LENGTH MUST BE EVEN TO SUBTRACT PAIRS
    nbaselines = np.int(len(x)/bl)
    nbl_long = np.int(len(x)/bl_long)

    if np.mod(bl_long,2) != 0:
        bl_safe = bl_long-1
    else:
        bl_safe = bl_long*1

    pairs = np.zeros(bl_safe/2)
    C_N = np.zeros(nbaselines,dtype='Float64')


    #FOR EACH BASELINE SUBTRACT INDEPENDENT PAIRS TO ESTIMATE WHITE NOISE
    for i in np.arange(nbl_long):

        #Pair up values within 1 baseline
        pairs = x[i*bl_safe:i*bl_safe + bl_safe-1] - x[i*bl_safe+1:i*bl_safe + bl_safe]

        #Sort pairs and find interquartile range
        psa = pairs.argsort()
        interquartile = np.std(pairs[psa[int(len(pairs)*0.25):int(len(pairs)*0.75)]])

        
        #Remove spikes
        tp = pairs[(np.abs(pairs) < interquartile*5.)]

        if (i+1)*bl_long/bl <= C_N.size: 
            C_N[i*bl_long/bl:(i+1)*bl_long/bl] = np.std(tp)/np.sqrt(2.)
        else:
            C_N[i*bl_long/bl:] = np.std(tp)/np.sqrt(2.)
    
    nans = (np.isnan(C_N) == 1)

    if nans.size > 0:
        notnans   = (np.isnan(C_N) == 0)
        C_N[nans] = MPI_tools.MPI_sum(comm,C_N[notnans])/MPI_tools.MPI_len(comm,C_N[notnans])

    zeros = (C_N**2 == 0)

    if zeros.size > 0:
        notzeros   = (C_N**2 > 0)
        C_N[zeros] = MPI_tools.MPI_sum(comm,C_N[notzeros])/MPI_tools.MPI_len(comm,C_N[notzeros])

    return C_N**2
