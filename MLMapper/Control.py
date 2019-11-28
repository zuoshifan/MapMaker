#Standard modules:
try:
    from mpi4py import MPI
    f_found=True
    from ..Tools import MPI_tools
except ImportError:
    f_found=False

import numpy as np


import scipy.fftpack as sfft
import scipy.interpolate as interp

#Map-making modules:
#from MapMaker.Tools import MPI_tools
from MapMaker.Tools.Mapping import MapsClass
from MapMaker.Tools import nBinning as Binning
from MapMaker.Tools import WhiteCovar

from MapMaker.CGM.nCGM import CGM

#Destriper modules:
from MLFuncs import bFunc,AXFunc


#------------------------------------------------------------------------#
#-----------------------------ML FUNCTIONS-------------------------------#
#------------------------------------------------------------------------#

def EstimateModel(resid):
    '''
    Return PSD of residual noise vector.

    Arguments
    resid -- Residual noise vector 
    '''

    #FFT's are faster if length of array is a power of 2
    # bit2 = 2**np.ceil(np.log10(resid.size)/np.log10(2))
    bit2 = 2**np.int(np.ceil(np.log10(resid.size)/np.log10(2)))

    #FFT the residual noise vector:
    fa = sfft.fft(resid,n=bit2)
    fabs = np.abs(fa)

    #Remove any very small spikes in the PSD:
    bd = np.where(fabs < fabs.mean()/1e10)[0]
    fabs[bd] = fabs[bd+1]

    #Bin PSD:
    bmdl = Binning.DownSample(fabs[0:fabs.size/2],(fabs.size/2)/20)
    bmid = Binning.DownSample(np.log10(np.arange(fabs.size/2)+1),(fabs.size/2)/20)

    #Remove the first few bins:
    end = 10.**np.mean(np.log10(bmdl[0:3]))
    bmdl = bmdl[3:]
    bmid = bmid[3:]
    bmdl[0] = end
    bmid[0] = 0.
    bmid[-1] = np.log10(fabs.size/2+1)

    #Interpolate between bins:
    mdlfunc = interp.interp1d(bmid,np.log10(bmdl),bounds_error=False,kind='linear')
    mdl     =(10.**mdlfunc(np.log10(np.arange(fabs.size/2)+1)))

    
    if mdl.size*2 == fabs.size:
        return np.append(mdl,mdl[::-1])
    else:
        return np.append(np.append(mdl,mdl[-1]),mdl[::-1])



def MLMapper(tod,pix,npix,comm=None,bl_long=None,Verbose=False,maxiter=3,cn=None):
    '''
    Return Destriped maps for a given set of TOD.

    Arguments
    tod -- Time-ordered single-dish radio telescope scanning data. (Array)
    pix -- Time-ordered pixel coordinate data. (Array)
    npix -- Total number of pixels for output map.

    Keyword Arguments
    bl_long -- Number of samples used to estimate weights of samples.
    comm    -- MPI.COMM_WORLD
    Verbose -- Output debugging information (Default: False).
    maxiter -- Number of iterations used to estimate noise covariance (Default: 3).
    
    '''

    pix = pix.astype('i')
    npix = int(npix)

    # Switch on MPI 
    if f_found:
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        rank = comm.Get_rank()
    else:
        comm = None
        rank = 0
        pass

    #If user provides no 'long' baselength, estimate it to be a quarter of the tod size or 2 samples.
    if bl_long == None: 
        bl_long = np.max([tod.size/4,2])


    tod -= np.median(tod)
    
    #Estimate white-noise level of the data:
    bl = 1
    null = np.ones(tod.size/bl_long + 1)
    #cn   = WhiteCovar.WhiteCovar(tod,bl,bl_long,comm=comm)
    if isinstance(cn,type(None)):
        cn = WhiteCovar.WhiteCovar(tod,bl_long,bl_long,comm=comm)
        # cn = np.repeat(cn,bl_long)
        cn = np.repeat(cn,bl_long+1)[:len(tod)]

    #Generate Maps:
    Maps = MapsClass(npix,rank=rank) #If rank == 0, generate root maps (swroot, hwroot)
    Maps.GoodPixels(pix)

    #Initial noise residual == input tod:
    resid = tod*1.

    #Initial guess at map vector:
    Binning.BinMap(tod,bl,pix,cn,Maps.m,
                   sw=Maps.sw,
                   hw=Maps.hw,
                   hits=Maps.hits,
                   swroot=Maps.swroot,
                   hwroot=Maps.hwroot,
                   hitmap=Maps.hitmap,                   
                   comm=comm)

    m0 = Maps.m[Maps.gd] #Save initial guess to separate variable
    

    for iteration in range(maxiter):
        if Verbose:
            print 'NOISE ITERATION: ', iteration

        #bit2 = 2**np.ceil(np.log10(resid.size)/np.log10(2))
        noiseModel = EstimateModel(resid)


        #Return the Residual fit to gain drifts:
        Maps.m[Maps.gd] = CGM(m0,bFunc,AXFunc,args=(tod,pix,noiseModel,null,cn,Maps,bl_long),comm=comm,Verbose=Verbose)

            
        mapMean = np.nanmean(Maps.m[Maps.gd])

        #Calculate new noise residual:
        resid = tod - (Maps.m[pix] - mapMean)

    return Maps
