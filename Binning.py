import numpy as np
from mpi4py import MPI
import MPI_tools

import fBinning
import time

def BinMap(x,bl,p,cn,m,sw=[None],hw=[None],hits=[None],swroot=None,hwroot=None,comm=None):
    '''Return data binned into map pixels.
    '''

    size = comm.Get_size()
    rank = comm.Get_rank()
    
    m[:] = 0.
    
    if sw[0]:
        sw[:]= 0.
    else:
        sw = np.zeros(len(m))

    if hw[0]:
        hw[:]= 0.        
    else:
        hw = np.zeros(len(m))

    if hits[0]:
        hits[:] = 0.


    #LOOP THROUGH EACH MAP PIXEL TO PRODUCE MAP (AND WEIGHT MAP)
    if hits[0]:
        sw[:],hw[:],hits[:] = fBinning.bin_pix_hits(x,p,cn,int(bl),sw,hw,hits)
    else:
        sw[:],hw[:] = fBinning.bin_pix(x,p,cn,int(bl),sw,hw)
            
    #Sum up on root node:
    comm.Reduce(sw  ,swroot  ,op=MPI.SUM,root=0)
    comm.Reduce(hw  ,hwroot  ,op=MPI.SUM,root=0)

    if hits[0]:
        comm.Reduce(hits,hitmap ,op=MPI.SUM,root=0)

    if rank==0:
        meh = (hwroot>0) & (np.isnan(hwroot) == 0) & (np.isinf(hwroot) ==0)
        m[meh] = swroot[meh]/hwroot[meh]
                

    # Then broadcast it back out to all nodes
    m[:] = comm.bcast(m,root=0)
    if hits[0]:
        hits[:] = comm.bcast(hitmap,root=0)
        

def BinMap_ext(a,bl,p,cn,m,sw=[None],hw=[None],swroot=None,hwroot=None,comm=None):
    '''Return baselines binned into a map.
    '''

    size = comm.Get_size()
    rank = comm.Get_rank()
        
    m[:] = 0.

    if sw[0]:
        sw[:]= 0.
    else:
        sw = np.zeros(len(m))

    if hw[0]:
        hw[:]= 0.        
    else:
        hw = np.zeros(len(m))


    #LOOP THROUGH EACH MAP PIXEL TO PRODUCE MAP (AND WEIGHT MAP)
    print p.shape,np.squeeze(a).shape
    sw[:],hw[:] = fBinning.bin_pix_ext(np.squeeze(a).astype('d'),p.astype('i'),cn.astype('d'),sw.astype('d'),hw.astype('d'),int(bl))
            
    #Sum up on root node:
    comm.Reduce(sw  ,swroot  ,op=MPI.SUM,root=0)
    comm.Reduce(hw  ,hwroot  ,op=MPI.SUM,root=0)
    
    if rank==0:
        meh = (hwroot>0) & (np.isnan(hwroot) == 0) & (np.isinf(hwroot) ==0)
        m[meh] = swroot[meh]/hwroot[meh]
                

    # Then broadcast it back out to all nodes
    m[:] = comm.bcast(m,root=0)
