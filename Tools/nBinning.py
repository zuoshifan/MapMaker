import numpy as np
try:
    from mpi4py import MPI
    f_found=True
    from ..Tools import MPI_tools
except ImportError:
    f_found=False

import fBinning
import time

def BinMap(x,bl,p,cn,m,sw=[None],hw=[None],hits=[None],swroot=None,hwroot=None,hitmap=None,comm=None,poorMask=None):
    '''Return data binned into map pixels.
    '''

    # Switch on MPI 
    if f_found:
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        rank = comm.Get_rank()
    else:
        rank = 0
        pass

    m[:] = 0.
    
    if sw[0] != None:
        sw[:]= 0.
    else:
        sw = np.zeros(len(m))

    if hw[0] != None:
        hw[:]= 0.        
    else:
        hw = np.zeros(len(m))

    if hits[0] != None:
        hits[:] = 0.


    #LOOP THROUGH EACH MAP PIXEL TO PRODUCE MAP (AND WEIGHT MAP)
    if hits[0] != None:
        if isinstance(poorMask, type(None)):
            sw[:],hw[:],hits[:] = fBinning.bin_pix_hits(x,p,cn,sw,hw,hits)
        else:
            sw[:],hw[:],hits[:] = fBinning.bin_pix_hits_pmask(x,p,cn,sw,hw,hits,poorMask.astype('i'))
    else: 
        sw[:],hw[:] = fBinning.bin_pix(x,p,cn,sw,hw)
            



    if f_found:
        comm.Reduce(sw  ,swroot  ,op=MPI.SUM,root=0)
        comm.Reduce(hw  ,hwroot  ,op=MPI.SUM,root=0)

        if hits[0] != None:
            comm.Reduce(hits,hitmap ,op=MPI.SUM,root=0)        
    else:
        swroot = sw*1.
        hwroot = hw*1.
        if hits[0] != None:
            hitmap = hits*1.


    if rank==0:
        meh = (hwroot > 0) & (np.isnan(hwroot) == 0) & (np.isinf(hwroot) == 0)
        m[meh] = swroot[meh]/hwroot[meh]
                

    # Then broadcast it back out to all nodes
    if f_found:
        m[:] = comm.bcast(m,root=0)
        if hits[0] != None:
            hits[:] = comm.bcast(hitmap,root=0)

def BinMapPol_Angs(x,bl,p,phi,cn,map,sw=[None],hw=[None],swroot=None,hwroot=None,hitmap=None,comm=None):
    '''Return data binned into map pixels.
    '''

#     size = comm.Get_size()
#     rank = comm.Get_rank()
    
#     map[:] = 0.
    
#     if sw[0] != None:
#         sw[:]= 0.
#     else:
#         sw = np.zeros(len(m))

#     if hw[0] != None:
#         hw[:]= 0.        
#     else:
#         hw = np.zeros(len(m))
    

    c2 = map*0. 
    s2 = map*0. 
    sc = map*0.
    sw = np.zeros(len(map))
    hw = np.zeros(len(map))

    #Loop through angle maps
    #noweights = cn*0. + 1.
    maps = [sc,s2,c2]
    funcs = [lambda p0:  np.sin(p0)*np.cos(p0),
             lambda p0:  np.sin(p0)*np.sin(p0),
             lambda p0:  np.cos(p0)*np.cos(p0)]

    for i,m in enumerate(maps):
        sw[:] = 0.
        hw[:] = 0.

        sw[:],hw[:] = fBinning.bin_pix(funcs[i](phi),p,cn,sw,hw)
        
        if i == 0:
            gd = np.where(hw != 0)[0]

        m[gd] = sw[gd]/hw[gd]
        

    map[:] = (c2 * s2 - sc**2 )        
        

def BinMap_ext(a,bl,p,cn,m,sw=[None],hw=[None],swroot=None,hwroot=None,comm=None,hits=[None],poorMask=None):
    '''Return baselines binned into a map.
    '''
    # Switch on MPI 
    if f_found:
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        rank = comm.Get_rank()
    else:
        rank = 0
        pass
        
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
    if isinstance(poorMask, type(None)):
        sw[:],hw[:] = fBinning.bin_pix_ext(a.astype('d'),p.astype('i'),cn.astype('d'),sw.astype('d'),hw.astype('d'),int(bl))
    else:
        print a.size,p.size,cn.size,sw.size,hw.size,poorMask.size
        sw[:],hw[:] = fBinning.bin_pix_ext_pmask(a.astype('d'),p.astype('i'),cn.astype('d'),sw.astype('d'),hw.astype('d'),poorMask.astype('i'),int(bl))

    
    #Sum up on root node:
    if f_found:
        comm.Reduce(sw  ,swroot  ,op=MPI.SUM,root=0)
        comm.Reduce(hw  ,hwroot  ,op=MPI.SUM,root=0)

        if hits[0] != None:
            comm.Reduce(hits,hitmap ,op=MPI.SUM,root=0)        
    else:
        swroot = sw*1.
        hwroot = hw*1.
        if hits[0] != None:
            hitmap = hits*1.

    if rank==0:
        meh = (hwroot > 0) & (np.isnan(hwroot) == 0) & (np.isinf(hwroot) == 0)
        m[meh] = swroot[meh]/hwroot[meh]
                
    # Then broadcast it back out to all nodes
    if f_found:
        m[:] = comm.bcast(m,root=0)
        if hits[0] != None:
            hits[:] = comm.bcast(hitmap,root=0)


def BinMap_with_ext(tod,a,bl,p,cn,m,sw=[None],hw=[None],swroot=None,hwroot=None,comm=None,hits=[None]):
    '''Return baselines binned into a map.
    '''
    # Switch on MPI 
    if f_found:
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        rank = comm.Get_rank()
    else:
        rank = 0
        pass
        
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
    sw[:],hw[:] = fBinning.bin_pix_with_ext(tod.astype('d'),np.squeeze(a).astype('d'),p.astype('i'),cn.astype('d'),sw.astype('d'),hw.astype('d'),int(bl))
            
    #Sum up on root node:
    #Sum up on root node:
    if f_found:
        comm.Reduce(sw  ,swroot  ,op=MPI.SUM,root=0)
        comm.Reduce(hw  ,hwroot  ,op=MPI.SUM,root=0)

        if hits[0] != None:
            comm.Reduce(hits,hitmap ,op=MPI.SUM,root=0)        
    else:
        swroot[:] = sw[:]
        hwroot[:] = hw[:]
        if hits[0] != None:
            hitmap[:] = hits[:]

    if rank==0:
        meh = (hwroot > 0) & (np.isnan(hwroot) == 0) & (np.isinf(hwroot) == 0)
        m[meh] = swroot[meh]/hwroot[meh]
                

    # Then broadcast it back out to all nodes
    if f_found:
        m[:] = comm.bcast(m,root=0)
        if hits[0] != None:
            hits[:] = comm.bcast(hitmap,root=0)

def DownSample(a,newlen,Errors=False):
    '''
    Return binned version of input array

    Arguments
    a -- Input array
    newlen -- Length of binned output array

    Keyword Arguments
    Errors -- Return errors on bin values (Default: False)


    Notes: The last datum of the output array may contain less samples
    than the rest of the data in the output array.

    '''

    if newlen > a.size:
        print 'WARNING: BIN SIZE GREATER THAN ARRAY SIZE'
        return None

    bins,errs = fBinning.downsample(a,newlen)

    if Errors:
        return bins, errs
    else:
        return bins
