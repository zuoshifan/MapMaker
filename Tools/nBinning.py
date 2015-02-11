import numpy as np
from mpi4py import MPI
import MPI_tools

from Polarisation.Tools import fBinning
import time

def BinMap(x,bl,p,cn,m,sw=[None],hw=[None],hits=[None],swroot=None,hwroot=None,hitmap=None,comm=None):
    '''Return data binned into map pixels.
    '''

    if type(comm) != type(None):
        size = comm.Get_size()
        rank = comm.Get_rank()
    else:
        rank = 0
    
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
        sw[:],hw[:],hits[:] = fBinning.bin_pix_hits(x,p,cn,int(bl),sw,hw,hits)
    else: 
        sw[:],hw[:] = fBinning.bin_pix(x,p,cn,int(bl),sw,hw)
            
    #Sum up on root node:
    if type(comm) != type(None):
        comm.Reduce(sw  ,swroot  ,op=MPI.SUM,root=0)
        comm.Reduce(hw  ,hwroot  ,op=MPI.SUM,root=0)

    if hits[0] != None:
        comm.Reduce(hits,hitmap ,op=MPI.SUM,root=0)

    if (type(comm) != type(None)) & (rank==0):
        meh = (hwroot>0) & (np.isnan(hwroot) == 0) & (np.isinf(hwroot) ==0)
        m[meh] = swroot[meh]/hwroot[meh]
    else:
        meh = (hw>0) & (np.isnan(hw) == 0) & (np.isinf(hw) ==0)
        m[meh] = sw[meh]/hw[meh]
        
                

    # Then broadcast it back out to all nodes
    if (type(comm) != type(None)):
        m[:] = comm.bcast(m,root=0)
        
    if hits[0] != None:
        hits[:] = comm.bcast(hitmap,root=0)

def BinMapPol_Angs(x,bl,p,phi,cn,Maps):
    '''Return data binned into map pixels.
    '''

    Maps.c2[:] = 0. 
    Maps.s2[:] = 0. 
    Maps.sc[:] = 0.



    #Loop through angle maps
    #noweights = cn*0. + 1.
    maps = [Maps.sc,Maps.s2,Maps.c2]
    funcs = [lambda p0:  np.sin(p0)*np.cos(p0),
             lambda p0:  np.sin(p0)*np.sin(p0),
             lambda p0:  np.cos(p0)*np.cos(p0)]

    for i,m in enumerate(maps):
        Maps.sw[:] = 0.
        Maps.hw[:] = 0.
        
        Maps.sw[:],Maps.hw[:] = fBinning.bin_pix(funcs[i](phi),p,cn,int(bl),Maps.sw,Maps.hw)
        
        if i == 0:
            gd = np.where(Maps.hw != 0)[0]

        m[gd] = Maps.sw[gd]/Maps.hw[gd]
        

    Maps.DetA[:] = (Maps.c2 * Maps.s2 - Maps.sc**2 )
    bd = np.where(np.abs(Maps.DetA[p]) < 0.00)[0]
    
    cn[bd] = 1e24
    


def BinMapPol(x,bl,p,phi,cn,Maps):

    Maps.q[:] = 0. #Qmap
    Maps.u[:] = 0. #Umap
    Maps.vs[:]  = 0.
    Maps.vc[:]  = 0.
   

    #Loop through data maps
    maps = [Maps.vc,Maps.vs]
    funcs = [lambda x0,p0:  x0*np.cos(p0),
             lambda x0,p0:  x0*np.sin(p0)]
    
    for i,m in enumerate(maps):
        Maps.sw[:] = 0.
        Maps.hw[:] = 0. 
        Maps.sw[:],Maps.hw[:] = fBinning.bin_pix(funcs[i](x,phi),p,cn,int(bl),Maps.sw,Maps.hw)

        if i == 0:
            gd = np.where(Maps.hw != 0)[0]

        m[gd] = Maps.sw[gd]/Maps.hw[gd]
        
        
    gd = np.where(np.abs(Maps.DetA) > 0.00)[0]

    Maps.q[gd]  =  Maps.vs[gd] * (Maps.c2[gd]) -\
                   Maps.vc[gd] * (Maps.sc[gd])
    Maps.q[gd] /=  Maps.DetA[gd]
 
    Maps.u[gd]  =   Maps.vs[gd] * (Maps.sc[gd]) +\
                   -Maps.vc[gd] * (Maps.s2[gd])
    Maps.u[gd] /=  Maps.DetA[gd]
    Maps.u *= -1.
        

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
    sw[:],hw[:] = fBinning.bin_pix_ext(np.squeeze(a).astype('d'),p.astype('i'),cn.astype('d'),sw.astype('d'),hw.astype('d'),int(bl))
            
    #Sum up on root node:
    comm.Reduce(sw  ,swroot  ,op=MPI.SUM,root=0)
    comm.Reduce(hw  ,hwroot  ,op=MPI.SUM,root=0)
    
    if rank==0:
        meh = (hwroot>0) & (np.isnan(hwroot) == 0) & (np.isinf(hwroot) ==0)
        m[meh] = swroot[meh]/hwroot[meh]
                

    # Then broadcast it back out to all nodes
    m[:] = comm.bcast(m,root=0)


def BinMap_with_ext(tod,a,bl,p,cn,m,sw=[None],hw=[None],swroot=None,hwroot=None,comm=None):
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
    sw[:],hw[:] = fBinning.bin_pix_with_ext(tod.astype('d'),np.squeeze(a).astype('d'),p.astype('i'),cn.astype('d'),sw.astype('d'),hw.astype('d'),int(bl))
            
    #Sum up on root node:
    comm.Reduce(sw  ,swroot  ,op=MPI.SUM,root=0)
    comm.Reduce(hw  ,hwroot  ,op=MPI.SUM,root=0)
    
    if rank==0:
        meh = (hwroot>0) & (np.isnan(hwroot) == 0) & (np.isinf(hwroot) ==0)
        m[meh] = swroot[meh]/hwroot[meh]
                

    # Then broadcast it back out to all nodes
    m[:] = comm.bcast(m,root=0)


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
