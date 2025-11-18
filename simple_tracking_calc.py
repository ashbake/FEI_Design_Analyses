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
plt.ion()

from matplotlib.ticker import (AutoMinorLocator)

plt.ion()
font = {'size'   : 14}
matplotlib.rc('font', **font)
plt.rcParams['font.size'] = '14'
plt.rcParams['font.family'] = 'sans'
plt.rcParams['axes.linewidth'] = '1.3'
fontname = 'DIN Condensed'

def get_track_req(so):
	# compute requirement. requirement is 0.2lambda/D in y band
	yband_wavelength       = 1020 # nm, center of y band
	tracking_requirement_arcsec = 206265 * 0.2 * yband_wavelength / (so.inst.tel_diam*10**9) 
	tracking_requirement_pixel  = tracking_requirement_arcsec/so.track.platescale

	return tracking_requirement_pixel


if __name__=='__main__':
	#load inputs
	configfile = './hispec_tracking.cfg'
	so         = load_object(configfile)
	cload      = fill_data(so,track_on=True)

	print(so.track.centroid_err)

	print(get_track_req(so))

