# calc tt rate vs mag for different SNR goals for tracking camera
import sys, os, matplotlib
import numpy as np
import matplotlib.pylab as plt
from scipy import interpolate

font = {'size'   : 14}
matplotlib.rc('font', **font)

sys.path.append('/Users/ashbake/Documents/Research/ToolBox/specsim/')

from specsim.objects import load_object
from specsim.load_inputs import fill_data

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
	#load inputs
	configfile = './hispec_tracking_0.3s_texp.cfg'
	so         = load_object(configfile)
	so.track.texp = 0.3
	so.stel.teff  = 2300
	bands         = ['JHgap','H']
	
	cload      = fill_data(so,track_on=True)

	mags        = np.arange(4,16,1)
	snrs        = np.zeros((len(mags),len(bands)))
	for i, mag in enumerate(mags):
		so.stel.mag = mag
		for j,band in enumerate(bands):
			so.track.band = band
			cload       = fill_data(so,track_on=True)
			snrs[i,j]    = so.track.snr

	#plot
	fig, ax = plt.subplots(1,1,figsize=(8,6))
	ax.plot(mags,snrs[:,0],'-o',lw=2,label=bands[0])
	ax.plot(mags,snrs[:,1],'-o',lw=2,label=bands[1])
	ax.set_xlabel('%s Magnitude'%so.filt.band)
	ax.set_ylabel('SNR')
	ax.legend()
	ax.grid()
	plt.title('Tracking SNR in 0.3 sec')
	plt.savefig('./plots/tracking_snr_vs_mag_%s.png'%so.track.band,dpi=300,bbox_inches='tight')

	

