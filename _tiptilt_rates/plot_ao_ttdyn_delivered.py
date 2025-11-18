# plot tt rms vs magnitude for different AO modes
import sys, os, matplotlib
import numpy as np
import matplotlib.pylab as plt
from scipy import interpolate

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


def load_AO(so):
	"""
	load the TT rms data
	"""
	data       = wfe_tools.load_WFE(so.ao.ho_wfe_set, so.ao.ttdynamic_set, so.obs.zenith_angle, so.tel.seeing_set)
	ao_modes   = np.array(list(data.keys()))
	tt_wfes = {}
	mags        = np.arange(4,16,1)

	for ao_mode in ao_modes:
		# get magnitude in band the AO mode is defined in 
		wfe_mag  = _get_band_mag(so,so.stel.vraw, so.stel.sraw, so.stel.model, so.stel.stel_file,'Johnson',data[ao_mode]['band'],so.stel.factor_0)
		color  = wfe_mag - so.stel.mag
		# interpolate over WFEs and sample HO and TT at correct mag
		#f_howfe    = interpolate.interp1d(data[ao_mode]['ho_mag'],data[ao_mode]['ho_wfe'], bounds_error=False,fill_value=10000)
		f_ttwfe    = interpolate.interp1d(data[ao_mode]['tt_mag'],data[ao_mode]['tt_wfe'], bounds_error=False,fill_value=10000)
		#ho_wfe  = float(f_howfe(mags + color))
		tt_wfe  = f_ttwfe(mags + color)
		tt_wfes[ao_mode] = tt_wfe

	return mags, tt_wfes

def plot_ao_tt():
	"""
	plot the TT rms vs magnitude for different AO modes
	"""
	mags, tt_wfes = load_AO(so)

	fig, ax = plt.subplots(1,1,figsize=(8,6))
	for ao_mode in tt_wfes.keys():
		ax.plot(mags,tt_wfes[ao_mode],label=ao_mode,lw=2)
	ax.set_xlabel('%s Magnitude'%so.filt.band)
	ax.set_ylabel('TT WFE (mas)')
	ax.axhline(6.5,ls='--',c='k',label='6.5 mas')
	ax.legend(fontsize=10)
	ax.grid()
	ax.set_ylim(0.1,30)
	ax.set_title('TT WFE vs Magnitude')
	plt.savefig('./plots/tt_wfe_vs_mag.png',dpi=300,bbox_inches='tight')

if __name__=='__main__':
	#load inputs
	configfile = './hispec_tracking.cfg'
	so         = load_object(configfile)
	cload      = fill_data(so,track_on=True)
	
	plot_ao_tt()
