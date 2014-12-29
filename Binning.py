import numpy as np
from mpi4py import MPI
import MPI_tools

import fBinning
import time

def BinMap(x,bl,p,cn,m,sw=None,hw=None,hits=None,swroot=None,hwroot=None,comm=comm):
    '''Return data binned into map pixels.
    '''

        size = comm.Get_size()
        rank = comm.Get_rank()
        
        m[:] = 0.

        if sw:
            sw[:]= 0.
        else:
            sw = np.zeros(len(m))

        if hw:
            hw[:]= 0.        
        else:
            hw = np.zeros(len(m))

        if hits:
            hits[:] = 0.


        #LOOP THROUGH EACH MAP PIXEL TO PRODUCE MAP (AND WEIGHT MAP)
        if hits:
            sw[:],hw[:],hits[:] = fBinning.bin_pix_hits(x,p,cn,int(bl),sw,hw,hits)
        else:
            sw[:],hw[:] = fBinning.bin_pix(x,p,cn,int(bl),sw,hw)
            
        #Sum up on root node:
        comm.Reduce(sw  ,swroot  ,op=MPI.SUM,root=0)
        comm.Reduce(hw  ,hwroot  ,op=MPI.SUM,root=0)

        if hits:
            comm.Reduce(hits,hitmap ,op=MPI.SUM,root=0)

        if rank==0:
            meh = (hwroot>0) & (np.isnan(hwroot) == 0) & (np.isinf(hwroot) ==0)
            m[meh] = swroot[meh]/hwroot[meh]
                

        # Then broadcast it back out to all nodes
        m[:] = comm.bcast(m,root=0)
        if hits:
            hits[:] = comm.bcast(hitmap,root=0)


def BinMap_ext(a,bl,p,cn,m,sw=None,hw=None,swroot=None,hwroot=None,comm=comm):
    '''Return baselines binned into a map.
    '''

        size = comm.Get_size()
        rank = comm.Get_rank()
        
        m[:] = 0.

        if sw:
            sw[:]= 0.
        else:
            sw = np.zeros(len(m))

        if hw:
            hw[:]= 0.        
        else:
            hw = np.zeros(len(m))


        #LOOP THROUGH EACH MAP PIXEL TO PRODUCE MAP (AND WEIGHT MAP)
        sw[:],hw[:] = fBinning.bin_pix_ext(a,p,cn,sw,hw,int(bl))
            
        #Sum up on root node:
        comm.Reduce(sw  ,swroot  ,op=MPI.SUM,root=0)
        comm.Reduce(hw  ,hwroot  ,op=MPI.SUM,root=0)

        if rank==0:
            meh = (hwroot>0) & (np.isnan(hwroot) == 0) & (np.isinf(hwroot) ==0)
            m[meh] = swroot[meh]/hwroot[meh]
                

        # Then broadcast it back out to all nodes
        m[:] = comm.bcast(m,root=0)
