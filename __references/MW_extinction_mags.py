# calc signal to noise
# max for the calcium H&K
import sys, os
import matplotlib
import numpy as np
import matplotlib.pylab as plt
from scipy import interpolate
import pandas as pd
from datetime import date
from astropy.io import fits

font = {'size'   : 14}
matplotlib.rc('font', **font)

plt.ion()

from matplotlib.ticker import (AutoMinorLocator)

plt.ion()
font = {'size'   : 14}
matplotlib.rc('font', **font)
plt.rcParams['font.size'] = '14'
plt.rcParams['font.family'] = 'sans'
plt.rcParams['axes.linewidth'] = '1.3'
fontname = 'DIN Condensed'



if __name__=='__main__':
	#########
	# plot
	filename = './data/milkyway_dense_001.fits' # Average for Milky Way using R_V=5.0.
	f = fits.open(filename) 
	data = f[1].data
	wv =  1/data['WAVELENGTH'] # loaded is 1/microns
	av_BV = data['Av/E(B-V)']	

	plt.figure()
	plt.plot(wv*1000,av_BV)	
	y=[980,1100]
	J=[1170,1327]
	H=[1490,1780]
	K=[1990,2460]
	plt.fill_between(y,0,11,alpha=0.4)
	plt.fill_between(J,0,11,alpha=0.4)
	plt.fill_between(H,0,11,alpha=0.4)
	plt.fill_between(K,0,11,alpha=0.4)
	plt.xlabel('Wavelength (nm)')
	plt.ylabel('B-V extinction')


