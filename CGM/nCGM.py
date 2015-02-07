'''
Generalised CGM implementation:

Solves: x for a linear equation Ax = b.

Inputs:

A function which describes Ax <-- Modify this vector -in place-.
An array which describes b <-- Returns a column vector.
An initial guess for x.

'''

import numpy as np
from mpi4py import MPI

from MapMaker.Tools import MPI_tools

def CGM(x0, bFunc, AXFunc, args=(None,), maxiter=200, comm = MPI.COMM_WORLD, Verbose=False):
    '''
    Returns x for a linear system Ax = b where A is a symmetric, positive-definite matrix.

    Arguments
    x0  -- Initial guess at x
    bFunc -- Function describing how to generate b.
    AXFunc -- Function describing how to generate Ax.
    
    Keyword Arguments
    args -- Optional extra arguments for bFunc and AXFunc
    maxiter -- Maximum iterations of CGM loop.
    comm -- MPI.COMM_WORLD
    Verbose -- Print out distance to solution for each CGM iteration.
    '''

    rank = comm.rank

    #First, determine value of b:
    b = bFunc(x0,*args,comm=comm)


    dsize = len(b)
    xsize = len(x0)

    #Ensure b and x0 are column vectors:
    b  = np.reshape(b,(b.size,1))
    x0 = np.reshape(x0,(xsize,1))

    #Initial guess at Ax:
    Ax = np.zeros((dsize,1)) #Define for inital guess
    Ad = np.zeros((dsize,1)) #Define for use in the CGM loop
    AXFunc(x0,Ax,*args,comm=comm)


    #MAKE COLUMN VECTORS:

    xi = np.reshape(x0,(xsize,1))
    r  = np.reshape(b - Ax,(dsize,1))
    d  = np.reshape(b - Ax,(dsize,1))

    #Initial Threshold:
    del0 = r.T.dot(r)
    del0  = MPI_tools.MPI_sum(comm,del0)

    for i in range(maxiter):

        #Generate a new conjugate search vector Ad using d:
        AXFunc(d,Ad,*args,comm=comm)

        #Calculate search vector:
        sumr2  = r.T.dot(r)
        sumdAd = d.T.dot(Ad)
        sumr2  = MPI_tools.MPI_sum(comm,sumr2)
        sumdAd = MPI_tools.MPI_sum(comm,sumdAd)

        alpha = sumr2/sumdAd
        
        xi = xi + alpha*d
        if np.mod(i+1,200) == 0:
            AXFunc(xi,Ad,*args,comm=comm)
            rnew = np.reshape(b - Ad,(dsize,1))
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

        if (rank == 0) & (Verbose):
            print 'iteration: ', i, 1e-15*del0/lim

        if  lim < 1e-15*del0:
            if rank == 0:
                print 'Limit reached'

            break

    return np.squeeze(xi)
