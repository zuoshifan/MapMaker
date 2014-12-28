import numpy as np
from mpi4py import MPI
from matplotlib import pyplot
import math

#------------------------------------------------------------------------#
# Some MPI utilities

def MPI_sum(comm,x):
    # Sum values from all nodes and broadcast back the result

    size = comm.Get_size()
    rank = comm.Get_rank()

    if not isinstance(x,(int,float)):
        try:
            xsum = math.fsum(x)
        except TypeError:
            xsum = math.fsum(x.tolist())
    else:
        xsum = x

    vals = comm.gather(xsum,root=0)

    if rank==0:
        s = math.fsum(vals)
    else:
        s = None
            
    s = comm.bcast(s,root=0)

    return s

def MPI_concat(comm,x):
    # Concatenate values from all nodes and broadcast back the result

    size = comm.Get_size()
    rank = comm.Get_rank()

    vals = comm.gather(x,root=0)
    
    if rank==0:
        s = np.concatenate(vals)
    else:
        s = None
            
    s = np.array(comm.bcast(s,root=0))

    return s

def MPI_len(comm,x):
    # Concatenate values from all nodes and broadcast back the result
    size = comm.Get_size()
    rank = comm.Get_rank()


    datasize = x.size
    vals = comm.gather(datasize,root=0)
    
    s=0
    if rank==0:
        print 'MPI_len: ', rank,vals
        s = math.fsum(vals)
        #for val in vals:
        #    s += val
        #    print s,val
    else:
        s = None
            
    s = comm.bcast(s,root=0)

    return s
