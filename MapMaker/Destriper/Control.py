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

    n  = int(np.ceil(tod.size/float(bl))) #Number of baselines
    a0 = np.zeros(n)
    
    for i in range(n): #Loop through median values of tod and estimate initial baseline values
        a0[i] = np.median(tod[i*bl:(i+1)*bl])

    return a0

def Destriper(tod,bl,pix,npix,comm=None,bl_long=None,Verbose=False,maxiter=300,Medians=False,mask=None,cn=None,BinOnly=False,ReturnOffsets=False,poorMask=None):
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
    #a0 = np.zeros(tod.size//bl)#InitGuess(tod,bl)
    

    #Estimate white-noise level of the data:
    if isinstance(cn,type(None)):
        cn = WhiteCovar.WhiteCovar(tod,bl,bl_long,comm=comm)
        print 'Created estimated CN'
        cn = np.concatenate( (np.repeat(cn,bl),np.ones(tod.size-bl*cn.size)*cn[-1]) )

    #if not isinstance(mask,type(None)):
    #    cn[(mask == False)] = 1e24


    #Want to calculate an estimate of ratio of 1/f noise to white noise

    nBaselines     = int(np.ceil(tod.size / float(bl)))
    nLongBaselines = int(np.ceil(tod.size / float(bl_long)))
    bl2long = bl_long // bl
    for i in range(nLongBaselines):
        a0[i*bl2long:(i+1)*bl2long] -= np.median(a0[i*bl2long:(i+1)*bl2long])


    MAD_a0 = np.median(np.abs(a0-np.median(a0))) * 1.4826 # Factor for MAD of normal distribution
    noiseRatio = 6.#np.max([MAD_a0**2/np.median(cn)/15.,1.])

    if rank == 0:
        print 'Subtracted medians'
        print 'NOISE RATIO:', rank, noiseRatio,MAD_a0,np.std(a0),np.median(cn),a0.size

    #Generate Maps:
    Maps = MapsClass(npix,rank=rank) #If rank == 0, generate root maps (swroot, hwroot)

    #Make an array of prior 'long' baselines from medians
    b = np.zeros(a0.size) #Variations around zero
    for i in range(nLongBaselines):
        if i == nLongBaselines-1:
            hi = b.size 
        else:
            hi = (i+1) * bl2long

        b[i*bl2long:hi] = np.median(tod[i*bl_long:(i+1)*bl_long])
   #     tod[i*bl_long:(i+1)*bl_long] -= np.median(tod[i*bl_long:(i+1)*bl_long])

    a0 = np.copy(b)


    if not BinOnly:

        if not Medians:
        #Return the Destriper derived values for baselines:
            a0[:] = CGM(a0,bFunc,AXFunc,args=(tod,bl,pix,cn,Maps,b,noiseRatio,poorMask),comm=comm,Verbose=Verbose,maxiter=maxiter)
        
            if isinstance(poorMask, type(None)):
                for i in range(len(a0)):
                    if i < a0.size-1:
                        hi = (i+1)*bl
                    else:
                        hi = tod.size
                    tod[i*bl:hi] -= a0[i]
            else:
                a0 = (np.repeat(a0,bl))[:tod.size]#np.concatenate( (np.repeat(a0,bl),np.ones(tod.size-bl*a0.size)*a0[-1]) )
                from scipy.interpolate import interp1d
                aind = np.arange(a0.size)
                aterp = interp1d(aind[poorMask],a0[poorMask],kind='nearest',bounds_error=False)
                a0 = aterp(aind)
                tod -= a0


                
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
