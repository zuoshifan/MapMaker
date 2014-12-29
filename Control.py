import numpy as np
from mpi4py import MPI
import MPI_tools

import WhiteCovar
import MapBinning
import Binning

#------------------------------------------------------------------------#
#----------------------DESTRIPER FUNCTIONS-------------------------------#
#------------------------------------------------------------------------#

def Fold(comm,a,baselength,pix,C_N,MM):
    #FOLD: THE PROCEDURE FOR TURNING AMPLITUDES INTO
    # THE FTZFa VECTOR.

    size = comm.Get_size()
    rank = comm.Get_rank()

    #STRETCH THE AMPLITUDES OUT

    fa = F(a,baselength)
    #BIN THIS INTO A MAP
    MM.MakeMap(comm,fa,pix,C_N)
    

    #BIN EACH PIXEL INTO A BASELINE

    asum  = MPI_tools.MPI_sum(comm,a)
    FtZFa = FtP(MM.m,pix,baselength,C_N,MM.hits) + asum
      
    #SUBTRACT OFF THE BASELINE LENGTH
    FtZFa = a*np.float(baselength)/C_N - FtZFa
 

    return FtZFa




def FtP(m,p,bl,C_N,hitmap):


    #x = np.zeros(len(p)/bl)

    #BIN MAP PIXELS INTO BASELINES
    #for ipix, pix in enumerate(p):
    #    x[ipix/bl] += m[pix]

    #x /= C_N

    limit = 1
    x = Binning.bin_to_baselines(m,p,bl,C_N,hitmap,limit)

    return x

def F(x,bl):
    #STRETCH BASELINES TO TOD LENGTH
    return np.repeat(np.squeeze(x),bl)

def Ft(x,bl,C_N):

    #BIN TOD TO BASELINES
    nbaselines = len(x)/bl
    out = np.zeros(nbaselines)
    for i in range(nbaselines.astype('Int32')):
        out[i] = np.sum(x[i*bl : (i+1)*bl])/ C_N[i]

    return out


def CGM(comm,tod,pix,a0,baselength,C_N,MM):


    size = comm.Get_size()
    rank = comm.Get_rank()

    #MAKE FtZd
    MM.MakeMap(comm,tod,pix,C_N)
    
    FtZd = FtP(MM.m,pix,baselength,C_N,MM.hits)
    Ftd = Ft(tod,baselength,C_N)
    FtZd = Ftd - FtZd

    mpisum = MPI_tools.MPI_sum(comm,FtZd)

    #Guess First FtZFa
    FtZFa = Fold(comm,a0,baselength,pix,C_N,MM)

    mpisum = MPI_tools.MPI_sum(comm,FtZFa)

    
    #MAKE COLUMN VECTORS
    r  = np.array([FtZd - FtZFa]).T
    d  = np.array([FtZd - FtZFa]).T
    a  = np.array([a0]).T

    #INITIAL THRESHOLD
    del0 = r.T.dot(r)
    del0 = MPI_tools.MPI_sum(comm,del0)

    maxiter = 200
    for i in range(maxiter):

        Ad = np.array([Fold(comm,np.squeeze(d.T),baselength,pix,C_N,MM)]).T
        
        sumr2 = r.T.dot(r)
        sumdAd = d.T.dot(Ad)

        sumr2  = MPI_tools.MPI_sum(comm,sumr2)
        sumdAd = MPI_tools.MPI_sum(comm,sumdAd)

        alpha = sumr2/sumdAd

        a = a + alpha*d

        if np.mod(i+1,200) == 0:
            rnew = np.array([FtZd - Fold(comm,np.squeeze(a.T),baselength,pix,C_N,MM)]).T
        else:
            rnew = r - alpha*Ad

        sumnewr2 = rnew.T.dot(rnew) 
        sumnewr2 = MPI_tools.MPI_sum(comm,sumnewr2)
        
        beta = sumnewr2/sumr2
        
        d = rnew + beta*d
        r = rnew*1.

        lim = rnew.T.dot(rnew)
        lim = MPI_tools.MPI_sum(comm,lim)

        asum = MPI_tools.MPI_sum(comm,a)
        if rank == 0:
            print 1e-15*del0/lim,asum,'iteration: ', i

        if  lim < 1e-15*del0:
            if rank == 0:
                print 'Limit reached'

            break

    return a

def Destriper(comm,tod,baselength,pix,npix,medians=False,jackknife=False):

    """
    Destriper:

    Execute the destriper functions. Returns the map and hit map.

    Variables:
    ----------

    Inputs : tod,baselength,pix,npix ,(pixout, npixout, feed)

    tod = input data in order of time.
    baselength = baseleng
    pix = image pixel numbers in order of time.
    npix = total number of map pixels

    Optional inputs:
    ---------------
    pixout = image pixel numbers in order of time for output map.
    npixout = total number of map pixels for output map.
    feed   = data to use when calculating hitmap.
    
    """

    # Switch on MPI 
    size = comm.Get_size()
    rank = comm.Get_rank()

    #Map making class:
    MM = MapBinning.MapMaking(npix,int(baselength))
    MM.SortPixels(pix)
    #Define data:
    baselength = np.double(baselength)
    a0 = np.zeros(tod.size/baselength)
    
    for i in np.arange(tod.size/baselength):
        tod[i*baselength:(i+1)*baselength] -= np.median(tod[i*baselength:(i+1)*baselength])
        a0[i] = np.median(tod[i*baselength:(i+1)*baselength])

    C_N = WhiteCovar.WhiteCovar(tod,baselength)


    if medians:
        a = a0
    else:
        a = CGM(comm,tod,pix,a0,baselength,C_N,MM)

    MM.MakeMap(comm,tod-F(a,baselength),pix,C_N)

    if jackknife:
        MM.MakeJackKnives(comm,tod-F(a,baselength),pix,C_N)

        return MM.m, MM.hits, MM.jk1, MM.jk2
    else:

        return MM.m, MM.hits
