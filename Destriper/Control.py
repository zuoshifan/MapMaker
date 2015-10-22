#Standard modules:
import numpy as np
try:
    from mpi4py import MPI
    f_found=True
    from ..Tools import MPI_tools
except ImportError:
    f_found=False

#Map-making modules:
from ..Tools.Mapping import MapsClass
from ..Tools import nBinning as Binning
from ..Tools import WhiteCovar

from ..CGM.nCGM import CGM 

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

def Destriper(tod,bl,pix,npix,comm=None,bl_long=None,Verbose=False,maxiter=300,Medians=False,mask=None,cn=None,BinOnly=False,ReturnOffsets=False):
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
    if f_found:
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        rank = comm.Get_rank()
    else:
        rank = 0
        pass

    #If user provides no 'long' baselength, estimate it to be 10*baselength or tod.size
    if bl_long == None: 
        bl_long = np.min([tod.size,int(bl) * 10])


    a0   = InitGuess(tod,bl)

    if not BinOnly:
        if f_found:
            tod -= MPI_tools.MPI_sum(comm,a0)/MPI_tools.MPI_len(comm,a0)
        else:
            tod -= np.median(tod)

    #Define inital guess:
    a0 = InitGuess(tod,bl)
    
    #Estimate white-noise level of the data:
    if isinstance(cn,type(None)):
        cn = WhiteCovar.WhiteCovar(tod,bl,bl_long,comm=comm)
        cn = np.repeat(cn,bl)
        
    if not isinstance(mask,type(None)):
        cn[(mask == False)] = 1e24

    #Generate Maps:
    Maps = MapsClass(npix,rank=rank) #If rank == 0, generate root maps (swroot, hwroot)

    #from matplotlib import pyplot
        #tod[(mask == True)] = 0.
    #pyplot.plot(tod,',')        
    #pyplot.plot(mask,'-r')
    #pyplot.show()

    if not Medians:
        #Return the Destriper derived values for baselines:
        a0[:] = CGM(a0,bFunc,AXFunc,args=(tod,bl,pix,cn,Maps),comm=comm,Verbose=Verbose,maxiter=maxiter)

    if not BinOnly:
        
        for i in range(len(a0)):
            tod[i*bl:(i+1)*bl] -= a0[i]

            
    #if not isinstance(mask,type(None)):

    
    Binning.BinMap(tod,bl,pix,cn,Maps.m,
                   sw=Maps.sw,
                   hw=Maps.hw,
                   hits=Maps.hits,                   
                   swroot=Maps.swroot,
                   hwroot=Maps.hwroot,
                   hitmap=Maps.hitmap,                                      
                   comm=comm)

    if rank == 0:
        if ReturnOffsets:
            return Maps, a0
        else:
            return Maps
    else:
        return None
