'''
Generalised CGM implementation:

Solves: x for a linear equation Ax = b.

Inputs:

A function which describes Ax <-- Modify this vector -in place-.
An array which describes b <-- Returns a column vector.
An initial guess for x.

'''

import numpy as np
from matplotlib import pyplot

try:
    from mpi4py import MPI
    f_found=True
    from ..Tools import MPI_tools
except ImportError:
    f_found=False

def CGM(x0, bFunc, AXFunc, args=(None,), maxiter=200, comm = None, Verbose=False):
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

    # Switch on MPI 
    if f_found:
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        rank = comm.Get_rank()
    else:
        rank = 0
        pass

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
    if f_found:
        del0  = MPI_tools.MPI_sum(comm,del0)
    else:
        del0 = np.sum(del0)


    lastLim = 1e32 #Initial value
    for i in range(maxiter):

        #Generate a new conjugate search vector Ad using d:
        AXFunc(d,Ad,*args,comm=comm)


        #Calculate search vector:
        sumr2  = r.T.dot(r)
        sumdAd = d.T.dot(Ad)        
        if f_found:
            sumr2  = MPI_tools.MPI_sum(comm,sumr2)
            sumdAd = MPI_tools.MPI_sum(comm,sumdAd)
        else:
            sumr2 = np.sum(sumr2)
            sumdAd = np.sum(sumdAd)

        alpha = sumr2/sumdAd
        
        xi = xi + alpha*d

        ##test = np.squeeze(xi) - args[-1]
        #print np.mean(np.sqrt(args[3]))*3.
        #pyplot.hist(test,bins=40)
        #pyplot.show()

        if np.mod(i+1,200) == 0:
            AXFunc(xi,Ad,*args,comm=comm)
            rnew = np.reshape(b - Ad,(dsize,1))
        else:
            rnew = r - alpha*Ad

        sumnewr2 = rnew.T.dot(rnew) 
        if f_found:
            sumnewr2 = MPI_tools.MPI_sum(comm,sumnewr2)
        else:
            sumnewr2 = np.sum(sumnewr2)

        beta = sumnewr2/sumr2

        d = rnew + beta*d
        r = rnew*1.

        #AXFunc(d,Ad,*args,comm=comm)

        lim = rnew.T.dot(rnew)
        if f_found:
            lim = MPI_tools.MPI_sum(comm,lim)
        else:
            lim = np.sum(lim)

        if f_found:
            asum = MPI_tools.MPI_sum(comm,xi)
        else:
            asum = np.sum(xi)

        if (rank == 0) & (Verbose):
            print 'iteration: ', i, 1e-11*del0/lim


        if lim < lastLim:
            xbest = np.copy(xi)
            lastLim = lim

        if  lim < 1e-11*del0:
            if rank == 0:
                print 'Limit reached'

            break
        comm.barrier()




    return np.squeeze(xbest)
