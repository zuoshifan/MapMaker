import numpy as np
from mpi4py import MPI
import MPI_tools

from nCGM import CGM
from Mapping import MapsClass
from DesFuncs import bFunc,AXFunc

import nBinning as Binning
import WhiteCovar

#------------------------------------------------------------------------#
#----------------------DESTRIPER FUNCTIONS-------------------------------#
#------------------------------------------------------------------------#

def InitGuess(tod,baselength):
    '''Estimate baseline values using medians as the Destripers initial guess.

    Arguments
    tod -- Input data from instrument
    baselenth -- Length of baselines in samples
    '''
    bl = int(baselength)

    n  = tod.size/bl #Number of baselines
    a0 = np.zeros(n)
    
    for i in range(n): #Loop through median values of tod and estimate initial baseline values
        a0[i] = np.median(tod[i*bl:(i+1)*bl])

    return a0

def Destriper(tod,bl,pix,npix,comm=MPI.COMM_WORLD,bl_long=None):
    '''Execute the destriper functions. Returns the map and hit map.
    
    '''
    # Switch on MPI 
    size = comm.Get_size()
    rank = comm.Get_rank()

    #If user provides no 'long' baselength, estimate it to be 10*baselength or tod.size
    if bl_long == None: 
        bl_long = np.min([tod.size,int(bl) * 10])


    tod -= np.median(tod)

    #Define inital guess:
    a0 = InitGuess(tod,bl)
    
    #Estimate white-noise level of the data:
    cn = WhiteCovar.WhiteCovar(tod,bl,bl_long)

    #Generate Maps:
    Maps = MapsClass(npix,rank=rank) #If rank == 0, generate extra maps

    #Return the Destriper derived values for baselines:
    ain = a0*1.
    a0[:] = CGM(a0,bFunc,AXFunc,args=(tod,bl,pix,cn,Maps),comm=comm)

    from matplotlib import pyplot
    pyplot.plot(np.squeeze(ain),np.squeeze(a0),',')
    pyplot.show()
    
    Binning.BinMap_with_ext(tod,np.squeeze(a0),bl,pix,cn,Maps.m,
                            sw=Maps.sw,
                            hw=Maps.hw,
                            swroot=Maps.swroot,
                            hwroot=Maps.hwroot,
                            comm=comm)

    if rank == 0:
        return Maps.m
    else:
        return None
