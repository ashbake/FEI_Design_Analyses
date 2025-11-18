# example tracking camera usage
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
from specsim import obs_tools
from specsim import functions
plt.ion()

from matplotlib.ticker import (AutoMinorLocator)

plt.ion()
font = {'size'   : 14}
matplotlib.rc('font', **font)
plt.rcParams['font.size'] = '14'
plt.rcParams['font.family'] = 'sans'
plt.rcParams['axes.linewidth'] = '1.3'
fontname = 'DIN Condensed'


def open_filter_files():
	"""
	Open various Files from vendors and plot
	to evaluate
	"""
	path = '/Users/ashbake/Documents/Research/Projects/HISPEC/Throughput_Budget/inputs/quoted_tests/'
	filename = 'materion_J_it0.TXT'
	l,f2,f1=np.loadtxt(path + filename).T
	plt.figure()
	plt.plot(l,f1)
	plt.plot(l,f2)
	plt.plot(so.track.xtransmit,so.track.ytransmit*100*so.track.bandpass)


def plot_KPIC_FilterProfiles(so):
	"""
	"""
	# pick filters to load
	filter_names = ['JHgapKPIC','Hkpic']
	amp_modes    = [0.9, 0.8] # mods to amplitude from DS PyPO modes
	filters = {}
	for i,filter_name in enumerate(filter_names):
		filters[filter_name],_ =  obs_tools.get_tracking_band(so.stel.v,band=filter_name)
		filters[filter_name] *= amp_modes[i]

	# integrate 
	jhgap_int = functions.integrate(so.stel.v,filters[filter_names[0]]*so.track.signal_spec)
	h_int     = functions.integrate(so.stel.v,filters[filter_names[1]]*so.track.signal_spec)

	print(jhgap_int/h_int)
	print(h_int/jhgap_int)

	# plot
	plt.figure()
	plt.plot(so.stel.v,filters[filter_names[0]]*so.track.signal_spec,label=filter_names[0])
	plt.plot(so.stel.v,filters[filter_names[1]]*so.track.signal_spec,label=filter_names[1])
	plt.xlabel('Wavelength [nm]')
	plt.ylabel('Flux (a.u.)')
	plt.title('H to 1450nm Ratio: %s'%(round(h_int/jhgap_int,1)))
	plt.grid()
	plt.subplots_adjust(left=0.2,bottom=0.15)
	plt.xlim(1400,1700)
	plt.savefig('KPIC_JHgap_test_filterFluxRatio.png')

if __name__=='__main__':
	#load inputs
	configfile = './hispec_tracking.cfg'
	so         = load_object(configfile)
	cload      = fill_data(so,track_on=True)

	plot_KPIC_FilterProfiles(so)
