import numpy as np

class MapsClass:

    def __init__(self,npix,rank=0):

        self.m = np.zeros(npix)
        self.sw = np.zeros(npix)
        self.hw = np.zeros(npix)
        self.hits = np.zeros(npix)
        
        if rank == 0:
            self.swroot = np.zeros(npix)
            self.hwroot = np.zeros(npix)
        else:
            self.swroot = None
            self.hwroot = None
