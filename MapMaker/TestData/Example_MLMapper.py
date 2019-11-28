#First import all the standard modules:
import numpy as np
import pyfits
import healpy as hp
from matplotlib import pyplot

#Now import the MLMapper:
from MapMaker.MLMapper import Control

#Open the test data:
hdu = pyfits.open('TestData.fits')
tod = hdu[1].data['TOD'][0,:,0]
pix = hdu[1].data['PIX'][0,:,0]

#Number of pixels:
npix = 786432

#The control function returns a class containing several maps.
Maps = Control.MLMapper(tod,pix,npix,Verbose=True)

#Plot the ML map:
hp.mollview(Maps.m,min=-10,max=20,title='ML Map')
pyplot.show()