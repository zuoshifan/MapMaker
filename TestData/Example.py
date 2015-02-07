#First import all the standard modules:
import numpy as np
import pyfits
import healpy as hp
from matplotlib import pyplot

#Now import the Destriper:
from MapMaker.Destriper import Control 

#Open the test data:
hdu = pyfits.open('TestData.fits')
tod = hdu[1].data['TOD'][0,:,0]
pix = hdu[1].data['PIX'][0,:,0]

#Baseline length:
bl = 250

#Number of pixels:
npix = 786432

#The control function returns a class containing several maps.
Maps = Control.Destriper(tod,bl,pix,npix,Verbose=True)

#Generate a naive binned map too
Meds = Control.Destriper(tod,tod.size,pix,npix,Medians=True,Verbose=True)


#Plot the Destriped map:
hp.mollview(Maps.m,min=-10,max=20,title='Destriped Map')
hp.mollview(Meds.m,min=-10,max=20,title='Naive Map')
pyplot.show()
