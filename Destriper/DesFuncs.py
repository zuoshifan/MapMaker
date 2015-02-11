#Holds all the function for solving Ax and b.
#Standard modules:
import numpy as np
from mpi4py import MPI

#Map-making modules:
from Polarisation.Tools.Mapping import MapsClass
from Polarisation.Tools import nBinning as Binning
from Polarisation.Tools import fBinning
from Polarisation.Tools import MPI_tools


def bFunc(a0,tod,bl,pix,cn,cn_mask,Maps,phi,comm=None):
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

   
    #Make the angle maps:
    Binning.BinMapPol_Angs(tod,bl,pix,phi,cn_mask,Maps)

    #Make noisey Q and U maps of tod:
    Binning.BinMapPol(tod,bl,pix,phi,cn_mask,Maps)
    #import healpy as hp
    #from matplotlib import pyplot
    #hp.mollview(Maps.q,rot=[83,22],max=1,min=-1)
    #hp.mollview(Maps.u,rot=[83,22],max=1,min=-1)
    #pyplot.show()

    FtZd = Ft(tod,bl,cn) - FtP(Maps,phi,pix,bl,cn_mask,Maps.hits)    
    return np.reshape(FtZd,(FtZd.size,1))


def AXFunc(a,FtZFa,tod,bl,pix,cn,cn_mask,Maps,phi,comm=None):
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
    Binning.BinMapPol(np.repeat(np.squeeze(a),bl),bl,pix,phi,cn_mask,Maps)

    #import healpy as hp
    #from matplotlib import pyplot
    #hp.mollview(Maps.q,rot=[83,22],max=1,min=-1)
    #hp.mollview(Maps.u,rot=[83,22],max=1,min=-1)
    ##pyplot.show()
    
    #Get some of the baselines, e.g. 1^T a :
    if comm:
        asum = MPI_tools.MPI_sum(comm,a)
    else:
        asum = np.nansum(a)

    #pyplot.plot(FtP(Maps,phi,pix,bl,cn_mask,Maps.hits),',')
    #pyplot.show()

        
    #Calculate weighted baselength values and subtract the sum of each baseline of pixels:
    cn1 = np.reshape(cn,(bl,cn.size/bl))
    cn1 = np.sum(cn1,axis=0)
    FtZFa[:,0] = np.squeeze(a)*np.float(bl)/cn1 - FtP(Maps,phi,pix,bl,cn_mask,Maps.hits) - asum


def FtP(m,phi,p,bl,cn,hits):
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
    x = fBinning.bin_to_baselines(m.q.astype('d')   ,
                                  m.u.astype('d')   ,
                                  phi.astype('d')   ,                                 
                                  p.astype('i')   ,                                  
                                  int(bl)         ,
                                  cn.astype('d')  ,
                                  hits.astype('i'),limit,phi.size/int(bl))

    return x

def Ft(x,bl,cn):
    '''
    Return bin data into baselines

    x -- tod to be binned into baselines
    bl -- baseline length
    cn -- white-noise variances vector
    '''

    #BIN TOD TO BASELINES
    n = len(x)/int(bl)
    out = np.zeros(n)
    cn1 = np.reshape(x/cn,(bl,cn.size/bl))
    out = np.sum(cn1,axis=0)
   
    #for i in range(n):
    #    out[i] = np.sum(x[i*bl : (i+1)*bl])/ cn1[i]

    return out
