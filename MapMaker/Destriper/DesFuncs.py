#Holds all the function for solving Ax and b.
#Standard modules:
import numpy as np
from scipy.interpolate import interp1d

try:
    from mpi4py import MPI
    f_found=True
    from ..Tools import MPI_tools
except ImportError:
    f_found=False

#Map-making modules:
from ..Tools.Mapping import MapsClass
from ..Tools import nBinning as Binning
from ..Tools import fBinning

import time

def bFunc(a0,tod,bl,pix,cn,Maps,b,noiseRatio,poorMask,comm=None):
    '''
    Returns solution for Ft Z d. 

    Arguments
    a0 -- Offsets, not used.
    tod -- input data
    bl  -- baseline length
    pix -- pixel coordinate vector (tod.size)
    cn -- estimated white-noise variance vector
    Maps -- Object holding map data
    
    '''
    if comm:
        size = comm.Get_size()
        rank = comm.Get_rank()
    else:
        rank = 0
        size = 1


    Binning.BinMap(tod,bl,pix,cn,Maps.m,
                   sw=Maps.sw,
                   hw=Maps.hw,
                   hits=Maps.hits,
                   swroot=Maps.swroot,
                   hwroot=Maps.hwroot,
                   hitmap=Maps.hitmap, 
                   poorMask=poorMask,
                   comm=comm)


    FtZd  = Ft(tod,bl,cn,poorMask=poorMask) 
    print 'step1',np.sum(FtZd)

    FtZd -= FtP(Maps.m,pix,bl,cn,Maps.hits,a0.size,poorMask=poorMask)
    print 'step2'
    FtZd += Ft_ext(b,bl,cn*noiseRatio**2,poorMask=poorMask)#b/(cn[::bl]*noiseRatio**2) * float(bl)
    #FtZd +=b/(cn[::bl]*noiseRatio**2) * float(bl)

    print 'DONE!'
    return np.reshape(FtZd,(FtZd.size,1))


def AXFunc(a,FtZFa,tod,bl,pix,cn,Maps,b,noiseRatio,poorMask,comm=None):
    '''
    Returns solution for Ft Z F a

    Arguments
    a0 -- Offsets for this CGM iteration.
    FtZFa -- This iterations guess at the conjugate vector to a0. Modified in place.

    tod -- input data
    bl  -- baseline length
    pix -- pixel coordinate vector (tod.size)
    cn -- estimated white-noise variance vector
    Maps -- Object holding map data
    

    '''
    #Make a map of the baselines:    
    Binning.BinMap_ext(a[:,0],bl,pix,cn,Maps.m,
                       sw=Maps.sw,
                       hw=Maps.hw,
                       swroot=Maps.swroot,
                       hwroot=Maps.hwroot,
                       poorMask=poorMask,
                       comm=comm)

    
    #Get some of the baselines, e.g. 1^T a :
    if not isinstance(comm,type(None)):
        asum = MPI_tools.MPI_sum(comm,a)
    else:
        asum = np.nansum(a)


    #Calculate weighted baselength values and subtract the sum of each baseline of pixels:
    
    #First two steps are calculating th
    FtZFa[:,0]  = Ft_ext(a[:,0],bl,cn,poorMask=poorMask) 
    FtZFa[:,0] -= FtP(Maps.m,pix,bl,cn,Maps.hits,a.size,poorMask=poorMask)


    
    #Now subtract the prior information (this is a gaussian prior)
    FtZFa[:,0] += Ft_ext(a[:,0],bl,cn*noiseRatio**2,poorMask=poorMask)#a[:,0]/(cn[::bl]*noiseRatio**2) * float(bl)
    #FtZFa[:,0] += a[:,0]/(cn[::bl]*noiseRatio**2) * float(bl)

    #FtZFa[:,0] -= asum


def FtP(m,p,bl,cn,hits,asize,poorMask=None):
    '''
    Returns stretched out map binned into baselines

    Arguments
    m -- map vector
    p -- pixel vector
    bl -- baseline length
    cn -- white-noise variances vector
    hits -- min hits (always set to 0)
    '''

    limit = 0

    if isinstance(poorMask,type(None)):
        x = fBinning.bin_to_baselines(m.astype('d')   ,
                                      p.astype('i')   ,
                                      int(bl)         ,
                                      cn.astype('d')  ,
                                      asize)
    else:
        x = fBinning.bin_to_baselines_pmask(m.astype('d')   ,
                                            p.astype('i')   ,
                                            int(bl)         ,
                                            poorMask.astype('i'),
                                            cn.astype('d')  ,
                                            asize           )
        

    return x

def Ft(x,bl,cn,poorMask=None):
    '''
    Return bin data into baselines

    x -- tod to be binned into baselines
    bl -- baseline length
    C_N -- white-noise variances vector
    '''

    #BIN TOD TO BASELINES
    #n = int(np.ceil(len(x)/float(bl)))
    #out = np.zeros(n)
    #for i in range(n):
    #    out[i] = np.sum(x[i*bl : (i+1)*bl]/cn[i*bl : (i+1)*bl])

    #BIN TOD TO BASELINES
    n = int(np.ceil(len(x)/float(bl)))
    if isinstance(poorMask,type(None)):
        out = fBinning.bin_ft(x,cn,bl,n)
    else:
        out = fBinning.bin_ft_pmask(x,cn,bl,n,poorMask.astype('i'))
    return out


    return out


def Ft_ext(x,bl,cn,poorMask=None):
    '''
    Return bin data into baselines

    x -- tod to be binned into baselines
    bl -- baseline length
    C_N -- white-noise variances vector
    '''

    #BIN TOD TO BASELINES
    if isinstance(poorMask,type(None)):
        out = fBinning.bin_ft_ext(x,cn,bl)
    else:
        print x.size,cn.size,bl,poorMask.size
        out = fBinning.bin_ft_ext_pmask(x,cn,bl,poorMask.astype('i'))
    return out
