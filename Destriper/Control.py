#Standard modules:
import numpy as np
from mpi4py import MPI

#Map-making modules:
from MapMaker.Tools import MPI_tools
from MapMaker.Tools.Mapping import MapsClass
from MapMaker.Tools import nBinning as Binning
from MapMaker.Tools import WhiteCovar

from MapMaker.CGM.nCGM import CGM

#Destriper modules:
from DesFuncs import bFunc,AXFunc


#------------------------------------------------------------------------#
#----------------------DESTRIPER FUNCTIONS-------------------------------#
#------------------------------------------------------------------------#

def InitGuess(tod,baselength):
    '''
    Return estimate of offsets values with medians.

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

def Destriper(tod,bl,pix,npix,comm=MPI.COMM_WORLD,bl_long=None,Verbose=False,maxiter=300,Medians=False):
    '''
    Return Destriped maps for a given set of TOD.

    Arguments
    tod -- Time-ordered single-dish radio telescope scanning data. (Array)
    pix -- Time-ordered pixel coordinate data. (Array)
    bl  -- Length of baselines used to Destripe the tod. (Integer)
    npix -- Total number of pixels for output map.

    Keyword Arguments
    bl_long -- Length of baseline used to estimate weights of samples.
    comm    -- MPI.COMM_WORLD
    Verbose -- Output debugging information (Default: False).
    
    '''
    # Switch on MPI 
    size = comm.Get_size()
    rank = comm.Get_rank()

    #If user provides no 'long' baselength, estimate it to be 10*baselength or tod.size
    if bl_long == None: 
        bl_long = np.min([tod.size,int(bl) * 10])


    a0   = InitGuess(tod,bl)
    tod -= MPI_tools.MPI_sum(comm,a0)/MPI_tools.MPI_len(comm,a0)

    #Define inital guess:
    a0 = InitGuess(tod,bl)
    
    #Estimate white-noise level of the data:
    cn = WhiteCovar.WhiteCovar(tod,bl,bl_long,comm=comm)

    #Generate Maps:
    Maps = MapsClass(npix,rank=rank) #If rank == 0, generate root maps (swroot, hwroot)

    if not Medians:
        #Return the Destriper derived values for baselines:
        a0[:] = CGM(a0,bFunc,AXFunc,args=(tod,bl,pix,cn,Maps),comm=comm,Verbose=Verbose,maxiter=maxiter)
    
    #Binning.BinMap_with_ext(tod,np.squeeze(a0),bl,pix,cn,Maps.m,
    #                        sw=Maps.sw,
    #                        hw=Maps.hw,
    #                        swroot=Maps.swroot,
    #                        hwroot=Maps.hwroot,
    #                        comm=comm)

    tod = tod-np.repeat(np.squeeze(a0),bl)

    #Estimate white-noise level of the data:
    cn = WhiteCovar.WhiteCovar(tod,bl,bl_long,comm=comm)
    
    Binning.BinMap(tod,bl,pix,cn,Maps.m,
                   sw=Maps.sw,
                   hw=Maps.hw,
                   swroot=Maps.swroot,
                   hwroot=Maps.hwroot,
                   comm=comm)

    if rank == 0:
        return Maps
    else:
        return None
