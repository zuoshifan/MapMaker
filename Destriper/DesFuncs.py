#Holds all the function for solving Ax and b.
#Standard modules:
import numpy as np
from mpi4py import MPI

#Map-making modules:
from ..Tools.Mapping import MapsClass
from ..Tools import nBinning as Binning
from ..Tools import fBinning
from ..Tools import MPI_tools


def bFunc(a0,tod,bl,pix,cn,Maps,comm=None):
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
                   comm=comm)

    
    FtZd = Ft(tod,bl,cn) - FtP(Maps.m,pix,bl,cn,Maps.hits)
    return np.reshape(FtZd,(FtZd.size,1))


def AXFunc(a,FtZFa,tod,bl,pix,cn,Maps,comm=None):
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
    Binning.BinMap_ext(a,bl,pix,cn,Maps.m,
                       sw=Maps.sw,
                       hw=Maps.hw,
                       swroot=Maps.swroot,
                       hwroot=Maps.hwroot,
                       comm=comm)
    
    #Get some of the baselines, e.g. 1^T a :
    if comm:
        asum = MPI_tools.MPI_sum(comm,a)
    else:
        asum = np.nansum(a)
        

    #Calculate weighted baselength values and subtract the sum of each baseline of pixels:    
    FtZFa[:,0] =  Ft_ext(np.squeeze(a),bl,cn) - FtP(Maps.m,pix,bl,cn,Maps.hits) - asum
    

def FtP(m,p,bl,cn,hits):
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

    x = fBinning.bin_to_baselines(m.astype('d')   ,
                                  p.astype('i')   ,
                                  int(bl)         ,
                                  cn.astype('d')  ,
                                  len(p)/int(bl))


    return x

def Ft(x,bl,cn):
    '''
    Return bin data into baselines

    x -- tod to be binned into baselines
    bl -- baseline length
    C_N -- white-noise variances vector
    '''

    #BIN TOD TO BASELINES
    n = len(x)/int(bl)
    out = np.zeros(n)
    for i in range(n):
        out[i] = np.sum(x[i*bl : (i+1)*bl]/cn[i*bl : (i+1)*bl])


    return out


def Ft_ext(x,bl,cn):
    '''
    Return bin data into baselines

    x -- tod to be binned into baselines
    bl -- baseline length
    C_N -- white-noise variances vector
    '''

    #BIN TOD TO BASELINES
    out = fBinning.bin_ft_ext(x,cn,bl)
    return out
