'''
Generalised CGM implementation:

Solves: x for a linear equation Ax = b.

Inputs:

A function which describes Ax <-- Should return a column vector.
An array which describes b.
An initial guess for x.

'''

import numpy as np
import MPI_tools
from mpi4py import MPI

def CGM(x0, b, AXClass, maxiter=200, comm = MPI.COMM_WORLD):

    rank = comm.rank

    dsize = b.size
    xsize = len(x0)
    #Make b a column vector:
    b = np.reshape(b,(b.size,1))
    x0 = np.reshape(x0,(xsize,1))

    #Initial guess at Ax:
    Ax = AXClass.AXFunc(x0)


    #MAKE COLUMN VECTORS:

    xi = np.reshape(x0,(xsize,1))
    r  = np.reshape(b - Ax,(dsize,1))
    d  = np.reshape(b - Ax,(dsize,1))

    #Initial Threshold:
    del0 = r.T.dot(r)
    
    for i in range(maxiter):

        Ad = AXClass.AXFunc(d) #A.dot(d)

        sumr2  = r.T.dot(r)
        sumdAd = d.T.dot(Ad)

        sumr2  = MPI_tools.MPI_sum(comm,sumr2)
        sumdAd = MPI_tools.MPI_sum(comm,sumdAd)

        alpha = sumr2/sumdAd
        
        xi = xi + alpha*d
        if np.mod(i+1,200) == 0:
            rnew = np.reshape(b - AXClass.AXFunc(xi),(dsize,1)) # b - A.dot(x)
        else:
            rnew = r - alpha*Ad

        sumnewr2 = rnew.T.dot(rnew) 
        sumnewr2 = MPI_tools.MPI_sum(comm,sumnewr2)

        beta = sumnewr2/sumr2

        d = rnew + beta*d
        r = rnew*1.

        lim = rnew.T.dot(rnew)
        lim = MPI_tools.MPI_sum(comm,lim)

        asum = MPI_tools.MPI_sum(comm,xi)

        if rank == 0:
            print 1e-15*del0/lim,asum,'iteration: ', i

        if  lim < 1e-15*del0:
            if rank == 0:
                print 'Limit reached'

            break

    return xi
