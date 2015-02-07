#Holds all the function for solving Ax and b.
#Standard modules:
import numpy as np
from mpi4py import MPI

#Map-making modules:
from MapMaker.Tools.Mapping import MapsClass
from MapMaker.Tools import nBinning as Binning
from MapMaker.Tools import fBinning
from MapMaker.Tools import MPI_tools


def bFunc(a0,tod,bl,pix,cn,Maps,comm=None):
    '''Returns solution for Ft Z d. 
    
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
    '''Returns solution for Ft Z F a
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
    FtZFa[:,0] = np.squeeze(a)*np.float(bl)/cn - FtP(Maps.m,pix,bl,cn,Maps.hits) - asum


def FtP(m,p,bl,cn,hits):

    limit = 0
    x = fBinning.bin_to_baselines(m.astype('d')   ,
                                  p.astype('i')   ,
                                  int(bl)         ,
                                  cn.astype('d')  ,
                                  hits.astype('i'),limit)

    return x

def Ft(x,bl,C_N):

    #BIN TOD TO BASELINES
    n = len(x)/int(bl)
    out = np.zeros(n)
    for i in range(n):
        out[i] = np.sum(x[i*bl : (i+1)*bl])/ C_N[i]

    return out
