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

def compute_rate_to_snr_cap(target_snr,so):
	"""
	compute the frame rate that still allows a given SNR
	
	inputs:
		target_snr: float
			goal SNR to reach
		so: object
			specsim object
	
	outputs:
		rate: float
			frame rate that can be run the given SNR
	"""
	if so.track.nphot/so.track.npix > so.track.saturation:
		return np.nan # error if saturated bc etc hack wont work

	# estimate etc from shot noise alone
	etc_shot = so.track.texp * (target_snr/np.sqrt(so.track.nphot))**2

	# use new shot noise etc estimate to scale to etc with noise
	# estimate etc from snr
	new_flux = so.track.nphot * (etc_shot/so.track.texp) # counts at new etc_shot texp
	new_noise = np.sqrt(so.track.noise**2 - so.track.nphot  + new_flux) # remove old shot, add new shot
	
	etc_scaled_snr = new_flux/ new_noise
	etc   = etc_shot * (target_snr/etc_scaled_snr)**2

	rate = 1/etc
	
	return rate # Hz

if __name__=='__main__':
	#load inputs
	configfile = './hispec_tracking.cfg'
	so         = load_object(configfile)
	cload      = fill_data(so,track_on=True)

	target_snrs = [20,50,100,500]
	mags        = np.arange(4,16,1)
	rates       = np.zeros((len(mags),len(target_snrs)))
	for i, mag in enumerate(mags):
		so.stel.mag = mag
		cload       = fill_data(so,track_on=True)
		for j, target_snr in enumerate(target_snrs):
			rates[i,j] = compute_rate_to_snr_cap(target_snr,so)


	fig, ax = plt.subplots(1,1,figsize=(8,6))
	ax.plot(mags,rates[:,0],label='SNR=20',lw=2)
	ax.plot(mags,rates[:,1],label='SNR=50',lw=2)
	ax.plot(mags,rates[:,2],label='SNR=100',lw=2)
	ax.plot(mags,rates[:,3],label='SNR=500',lw=2)
	ax.axhline(500,ls='--',color='k')
	ax.set_xlabel('%s Magnitude'%so.filt.band)
	ax.set_ylabel('Frame Rate (Hz)')
	ax.legend()
	ax.set_yscale('log')
	ax.grid()
	ax.set_ylim(0.13,1.6*10**5)
	ax.set_title('%s Tracking, T=%s'%(so.track.band,so.stel.teff))
	plt.savefig('./plots/tracking_framerate_vs_mag.png',dpi=300,bbox_inches='tight')

	

