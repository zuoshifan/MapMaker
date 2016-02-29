import numpy as np

class MapsClass:

    def __init__(self,npix,rank=0):

        self.m = np.zeros(npix, dtype='f')
        self.sw = np.zeros(npix, dtype='f')
        self.hw = np.zeros(npix, dtype='f')
        self.hits = np.zeros(npix, dtype='f')
        
        if rank == 0:
            self.swroot = np.zeros(npix, dtype='f')
            self.hwroot = np.zeros(npix, dtype='f')
            self.hitmap = np.zeros(npix, dtype='f')            
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
