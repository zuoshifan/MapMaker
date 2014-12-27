import numpy as np
from mpi4py import MPI
import MPI_tools

import WhiteCovar
import MapBinning
import Binning

from CGM import CGM


class CGMSolver:

    def __init__(self,comm=MPI.COMM_WORLD):
    def AXFunc(x):


class Destriper:

    """ Create instance of destriping map-maker to act on a set of data.
    
    
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

    def __init__(self,tod,pix,baselength,bl_long,npix,medians=False,comm=MPI.COMM_WORLD):

        self.size = comm.Get_size()
        self.rank = comm.Get_rank()
        
        self.MM = MapBinning.MapMaking(npix,int(baselength))
        self.tod = tod
        self.pix = pix

        #Define data:
        self.baselength = np.double(baselength)


        #Make guess of 1/f using long baselines:
        self.a0 = np.zeros(self.tod.size/baselength)


        lastmed = 0.
        self.tod -= np.median(self.tod)
        for i in np.arange(self.tod.size/bl_long):
            if i < self.tod.size/bl_long - 1:
                lastmed = np.median(self.tod[i*bl_long:(i+1)*bl_long])
                self.a0[i*bl_long/baselength:(i+1)*bl_long/baselength] = np.median(self.tod[i*bl_long:(i+1)*bl_long])
            else:
                self.a0[i*bl_long/baselength:] = self.a0[i*bl_long/baselength-1]


        self.C_N = WhiteCovar.WhiteCovar(self.tod,self.baselength,bl_long)


        if medians:
            self.a = self.a0
        else:
            self.a = CGM(comm,tod,pix,a0,baselength,C_N,MM)
