import numpy as np
from mpi4py import MPI
import MPI_tools

import Binning
import time

class MapMaking:

    def __init__(self,npix,bl,comm=MPI.COMM_WORLD):
        '''
        '''
        self.bl = bl
        self.npix = npix

        self.comm = comm
        
        self.sw  = np.zeros(npix,dtype='float64')
        self.hw  = np.zeros(npix,dtype='float64')
        self.m   = np.zeros(npix,dtype='float64')
        self.nmap= np.zeros(npix,dtype='float64')
        
        self.jk1 = np.zeros(npix,dtype='float64')
        self.jk2 = np.zeros(npix,dtype='float64')
        self.hits= np.zeros(npix,dtype='float64')

        if comm.Get_rank() == 0:
            self.swmap = np.zeros(npix)
            self.hwmap = np.zeros(npix)
            self.hitmap= np.zeros(npix)
        else:
            self.swmap = None
            self.hwmap = None
            self.hitmap= None
            

    def SortPixels(self,p):
        '''
        SortPixels(p)

        p = time ordered pixel number array

        produces:

        self.ips : index of sorted pixel number array
        self.upix: unique pixel number array
        self.nupix: array of hits per unique pixel number
        '''

        #SORT TOD IN TERMS OF PIXEL NUMBER
        self.ips = p.argsort()  

        #GET THE NUMBER OF UNIQUE PIXELS AND UNIQUE PIXEL VALUES FOR MAP
        self.upix,self.iupix = np.unique(p[self.ips],return_index=1)
        nupix = np.bincount(p[self.ips])
        self.nupix = nupix[np.squeeze(np.where( (nupix > 0)))]




    def MakeMap(self,x,p,C_N):
        '''
        MakeMap(comm,x,p,C_N)
        '''

        size = self.comm.Get_size()
        rank = self.comm.Get_rank()
        
        self.m[:]    = 0.
        self.hits[:] = 0.
        self.nmap[:] = 0.

        self.sw[:]   = 0.
        self.hw[:]   = 0.        

        #LOOP THROUGH EACH MAP PIXEL TO PRODUCE MAP (AND WEIGHT MAP)

        self.sw[:],self.hw[:],self.hits[:] = Binning.bin_pix_hits(x,p,C_N,self.sw,self.hw,self.hits,int(self.bl))

        start = time.time()
        #Gather hits maps back to root to make full hits map
        #self.swmaps[:,:]  = self.comm.gather(self.sw  ,root=0)
        #self.hwmaps[:,:]  = self.comm.gather(self.hw  ,root=0)
        #self.hitmaps[:,:] = self.comm.gather(self.hits,root=0)
        
        #Sum up on root node:
        self.comm.Reduce(self.sw  ,self.swmap  ,op=MPI.SUM,root=0)
        self.comm.Reduce(self.hw  ,self.hwmap  ,op=MPI.SUM,root=0)
        self.comm.Reduce(self.hits,self.hitmap,op=MPI.SUM,root=0)

        print 'Processing:',time.time()-start

        if rank==0:
            #self.sw[:] = 0.
            #self.hw[:] = 0.
            #self.hits[:] = 0.
            #for map in self.swmaps:
            #    self.sw   += map

            #for map in self.hwmaps:
            #    self.hw   += map

            #for map in self.hitmaps:
            #    self.hits += map

            print self.swmap.shape
            print self.m.shape
            meh = np.where((self.hwmap>0) & (np.isnan(self.hwmap) == 0) & (np.isinf(self.hwmap) ==0) )[0]
            self.m[meh] = self.swmap[meh]/self.hwmap[meh]
            
            #del swmaps, hwmaps
    

        # Then broadcast it back out to all nodes
        self.m[:]    = self.comm.bcast(self.m,root=0)
        self.hits[:] = self.comm.bcast(self.hitmap,root=0)

    def MakeJackKnives(self,x,p,C_N):
        '''
        MakeMap(comm,x,p,C_N)
        '''

        if self.ips != None:

            size = self.comm.Get_size()
            rank = self.comm.Get_rank()

            self.m    *= 0.
            self.hits *= 0.
            self.nmap *= 0.

            self.sw   *= 0.
            self.hw   *= 0.

            self.jk1 *= 0.
            self.jk2 *= 0.

            #LOOP THROUGH EACH MAP PIXEL TO PRODUCE MAP (AND WEIGHT MAP)
            Binning.bin_pix_hits(x,self.upix,self.nupix,self.ips,C_N,self.sw,self.hw,self.nmap,self.hits,int(self.bl))

            
            #Gather hits maps back to root to make full hits map
            swmaps  = self.comm.gather(self.sw  ,root=0)
            hwmaps  = self.comm.gather(self.hw  ,root=0)

            hitmaps = self.comm.gather(self.hits,root=0)

            if rank==0:
                self.sw   *= 0.
                self.hw   *= 0.
                self.jk1  *= 0.
                self.jk2  *= 0.

                
                for imap in range(len(swmaps[0:size/2])):
                    self.sw   += swmaps[imap]

                for imap in range(len(hwmaps[0:size/2])):
                    self.hw   += hwmaps[imap]
              
                meh = np.where(self.hw>0)[0]
                self.jk1[meh] = self.sw[meh]/self.hw[meh]

                self.sw   *= 0.
                self.hw   *= 0.

                for imap in range(len(swmaps[size/2:])):
                    self.sw   += swmaps[imap]

                for imap in range(len(hwmaps[size/2:])):
                    self.hw   += hwmaps[imap]
              
                meh = np.where(self.hw>0)[0]
                self.jk2[meh] = self.sw[meh]/self.hw[meh]

                
                del swmaps, hwmaps
    

            # Then broadcast it back out to all nodes
            self.jk1 = self.comm.bcast(self.jk1,root=0)
            self.jk2 = self.comm.bcast(self.jk2,root=0)

        else:

            print 'RUN SORTPIXELS FIRST!'
