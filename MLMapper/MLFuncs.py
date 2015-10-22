#Holds all the function for solving Ax and b.
#Standard modules:
import numpy as np
try:
    from mpi4py import MPI
    f_found=True
    from ..Tools import MPI_tools
except ImportError:
    f_found=False

import scipy.fftpack as sfft

#Map-making modules:
from MapMaker.Tools.Mapping import MapsClass
from MapMaker.Tools import nBinning as Binning
from MapMaker.Tools import fBinning


def InvertCirculant(psd,tod):
    '''
    Return dot product of inverse circulant matrix N with vector d.

    psd -- First row of circulant matrix N
    tod -- Vector containing time-ordered data d.

    Notes: Instead of flipping the whole matrix N, only need to flip the first row.
           This means a circulant matrix can be stored as just a single row vector.
           Nf = np.fft.fft(np.fft.fft(N[::1])
    '''
    
    #Zero pad fft's to speed them up a bobbin:
    bit2 = 2**np.ceil(np.log10(tod.size)/np.log10(2))
    ftod  = sfft.fft(tod.astype(np.float64),n=bit2)/psd
    
    #Divide each fft(d) by fft(N), inverse transform and shift vector by -1 (for reasons unknown)
    temp = np.real(sfft.ifft(ftod))*float(tod.size)**2

    return temp[0:tod.size]



def bFunc(m0,tod,pix,model,null,cn,Maps,bl,comm=None):
    '''
    Return solution for Pt N d. 
    
    '''
    if comm:
        size = comm.Get_size()
        rank = comm.Get_rank()
    else:
        rank = 0
        size = 1


    Nd = InvertCirculant(model,tod)

    print cn.shape

    Binning.BinMap(Nd,bl,pix,cn,Maps.m,
                   sw=Maps.sw,
                   hw=Maps.hw,
                   hits=Maps.hits,
                   swroot=Maps.swroot,
                   hwroot=Maps.hwroot,
                   hitmap=Maps.hitmap,                   
                   comm=comm)
    
    return np.reshape(Maps.m[Maps.gd],(Maps.gd.size,1))


def AXFunc(m0,PtNPm,tod,pix,model,null,cn,Maps,bl,comm=None):
    '''
    Return solution for Pt N P m
    '''
    
    if comm:
        size = comm.Get_size()
        rank = comm.Get_rank()
    else:
        rank = 0
        size = 1


    #Set map to new map iteration:
    Maps.m[Maps.gd] = m0[:,0]

    #Sample this iteration
    NPm = InvertCirculant(model,Maps.m[pix])

    Binning.BinMap(NPm,bl,pix,cn,Maps.m,
                   sw=Maps.sw,
                   hw=Maps.hw,
                   hits=Maps.hits,
                   swroot=Maps.swroot,
                   hwroot=Maps.hwroot,
                   hitmap=Maps.hitmap,                   
                   comm=comm)
    
        
    #Needs to return as a column vector
    PtNPm[:,0] = Maps.m[Maps.gd]
