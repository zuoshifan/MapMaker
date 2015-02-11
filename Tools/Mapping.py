import numpy as np

class MapsClass:

    def __init__(self,npix,rank=0):

        self.m = np.zeros(npix)
        self.sw = np.zeros(npix)
        self.hw = np.zeros(npix)
        self.hits = np.zeros(npix)

        self.q = np.zeros(npix)
        self.u = np.zeros(npix)

        self.c2 = np.zeros(npix)
        self.s2 = np.zeros(npix)
        self.sc = np.zeros(npix)

        self.vs = np.zeros(npix)
        self.vc = np.zeros(npix)
        self.DetA = np.zeros(npix)

        
        if rank == 0:
            self.swroot = np.zeros(npix)
            self.hwroot = np.zeros(npix)
            self.hitmap = np.zeros(npix)            
        else:
            self.swroot = None
            self.hwroot = None
            self.hitmap = None       

    def GoodPixels(self,pix):
        '''
        Save unique pixel indexes in map
        
        '''

        self.gd = np.unique(np.sort(pix))
        self.gd = self.gd.astype('i')
