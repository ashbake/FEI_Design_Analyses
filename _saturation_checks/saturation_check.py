# goal is to see where tracking camera saturations
# and to inform ND filter needs
import sys, os
import matplotlib
import numpy as np
import matplotlib.pylab as plt
from scipy import interpolate
import pandas as pd
from datetime import date

font = {'size'   : 14}
matplotlib.rc('font', **font)

sys.path.append('/Users/ashbake/Documents/Research/ToolBox/specsim/')

from specsim.objects import load_object
from specsim.load_inputs import fill_data, _get_band_mag
from specsim import wfe_tools
plt.ion()

from matplotlib.ticker import (AutoMinorLocator)

plt.ion()
font = {'size'   : 14}
matplotlib.rc('font', **font)
plt.rcParams['font.size'] = '14'
plt.rcParams['font.family'] = 'sans'
plt.rcParams['axes.linewidth'] = '1.3'
fontname = 'DIN Condensed'


def saturation_per_mag(so):
	"""
	Per stellar mag, check where the 
	tracking camera saturates. Set exp
	time to 0.1 - allows for moderate exposure
	time in cases people dont want to track fast

	Inputs:
	-------
	so - storage object for specsim

	Outputs:
	-------
	mag_arr - magnitude array
	saturation_arr - flags if saturated
	"""
	mag_arr        = np.arange(3,15)
	exp_times      = np.array([1e-3,1e-2,1e-1,1,10])
	saturation_arr = np.zeros((len(mag_arr),len(exp_times)))
	nphot          = np.zeros((len(mag_arr),len(exp_times)))
	npix           = np.zeros((len(mag_arr)))
	for i,mag in enumerate(mag_arr):
		cload.set_mag(so,mag,trackonly=True)
		npix[i] = so.track.npix
		for j, exp_time in enumerate(exp_times):
			so.track.texp=exp_time
			cload.tracking(so)
			saturation_arr[i,j] = so.track.saturation_flag
			nphot[i,j] = so.track.nphot_nocap # will cap when saturated

	return mag_arr,exp_times,saturation_arr, nphot, npix

def plot_saturation_limits(so, mags, texps, sat, nphot, npix, saturation_level=40000):
	"""
	plot saturation limits

	inputs:
	-------
	so    - storage object for specsim run (uses stel.teff, so.filt.band)
	mags  - magnitudes, array of length mags
	texps - exposure times, array of length texps
	sat   - saturation flags, array of length mags x texps
	nphot - number of photons, 2D array of length mags x texps
	npix  - number of pixels, array of length mags
	saturation_level - level at which saturation occurs, default 40000

	outputs:
	--------
	plot of saturation limits saved to /plots/ folder
	"""
	# find magnitude at which things saturate
	ODs = [0,1,2,4,6]
	sat_mags = np.zeros((len(ODs),len(texps)))
	for i,OD in enumerate(ODs):
		for j,texp in enumerate(texps):
			# find index where saturated
			#isat = np.where(np.diff(sat[:,j])==-1)[0]
			bool_saturated = 10**(-1 * OD) * nphot[:,j]/npix > saturation_level
			isat = np.where(np.diff(bool_saturated))[0] # picks out index where becomes saturated
			if len(isat)==0: # if no or all saturated..
				if bool_saturated[0]==True: 
					sat_mags[i,j] = np.max(mags) 
				else:
					sat_mags[i,j] = np.min(mags) #never 
			else: # otherwise just save the mag of saturation
				sat_mags[i,j] = mags[isat]
	# PLOT
	fig, axs = plt.subplots(1,1,figsize=(7,6),sharex=True,sharey=True)
	for i,OD in enumerate(ODs):
		p = axs.semilogx(texps,sat_mags[i],label='OD%s'%OD)
		axs.fill_between(texps,0,sat_mags[i],alpha=0.3,color=p[0].get_color())

	plt.legend()
	plt.xlabel('Exposure Time (s) ')
	plt.ylabel('%s Magnitude'%so.filt.band)
	plt.ylim(np.max(mags),np.min(mags))
	plt.xlim(np.min(texps),np.max(texps))
	plt.title(str(int(so.stel.teff)) + 'K')
	plt.savefig('./plots/ATC_saturation_OD_%s.png'%so.stel.teff)


if __name__=='__main__':
	#load inputs
	configfile = '../hispec_tracking.cfg'
	so         = load_object(configfile)
	so.stel.teff = 5800
	so.stel.mag = 8
	so.track.texp = 0.05
	cload      = fill_data(so,track_on=True)

	mags, texps, sat, nphot,npix = saturation_per_mag(so)
	
	plot_saturation_limits(so, mags, texps, sat, nphot, npix, saturation_level=40000)


