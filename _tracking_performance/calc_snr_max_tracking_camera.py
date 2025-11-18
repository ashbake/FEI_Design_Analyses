# calc HISPEC FEI ATC SNR
import sys, os
import matplotlib
import numpy as np
import matplotlib.pylab as plt
from scipy import interpolate
from datetime import date

font = {'size'   : 14}
matplotlib.rc('font', **font)

sys.path.append('/Users/ashbake/Documents/Research/ToolBox/specsim/')

from specsim.objects import load_object
from specsim.load_inputs import fill_data, load_filter,get_band_mag
from specsim.functions import *

plt.ion()

from matplotlib.ticker import (AutoMinorLocator)

plt.ion()
font = {'size'   : 14}
matplotlib.rc('font', **font)
plt.rcParams['font.size'] = '14'
plt.rcParams['font.family'] = 'sans'
plt.rcParams['axes.linewidth'] = '1.3'
fontname = 'DIN Condensed'

############## RUN FUNCTIONS
def get_track_req(so):
	# compute requirement. requirement is 0.2lambda/D in y band
	yband_wavelength       = 1020 # nm, center of y band
	tracking_requirement_arcsec = 206265 * 0.2 * yband_wavelength / (so.inst.tel_diam*10**9) 
	tracking_requirement_pixel  = tracking_requirement_arcsec/so.track.platescale

	return tracking_requirement_pixel

def run_mags_exptimes(so,magarr,exptimes,output_path = './output/centroid_arrs/'):
	"""
	compute centroid error for a range of magnitudes and exposure times

	inputs:
	-------
	so: object
		config object
	magarr: array
		magnitude array
	exptimes: array	
		exposure times in seconds
	output_path: str
		path to save output
	
	outputs:
	--------
	none (saves files to output_path)
	"""
	centroid  = np.zeros((len(magarr),len(exptimes)))
	fwhms     = np.zeros((len(magarr)))
	signal    = np.zeros((len(magarr),len(exptimes)))
	noise     = np.zeros((len(magarr),len(exptimes)))
	strehl    = np.zeros((len(magarr),len(exptimes)))
	skyinstbkg       = np.zeros((len(magarr),len(exptimes)))
	for i,mag in enumerate(magarr): # this is the magnitude in filter band
		cload.set_filter_band_mag(so,'H','2mass',mag,trackonly=True)
		fwhms[i] = so.track.fwhm
		if np.isnan(so.track.noise) or np.isinf(so.track.fwhm): continue
		for j, texp in enumerate(exptimes):
			cload.set_tracking_band_texp(so,so.track.band,texp)
			centroid[i,j]   = float(so.track.centroid_err)
			signal[i,j]     = float(so.track.signal)
			noise[i,j]      = float(so.track.noise)
			strehl[i,j]     = float(so.track.strehl)
			skyinstbkg[i,j] = float(so.track.inst_bg_ph + so.track.sky_bg_ph)
	
	# save centroid array for stellar temp, magntiude, ao mode, camera
	savename = 'cendata_%s_ao_%s_band_%s_teff_%s_%smag_fieldr_%sarcsec' %(so.track.camera,so.ao.mode,so.track.band,so.stel.teff,so.filt.band,so.track.field_r)
	np.save(output_path + 'centroid_%s'%savename,centroid)
	np.save(output_path + 'fwhm_%s'%savename,fwhms)
	np.save(output_path + 'snr_%s'%savename,signal/noise)
	np.save(output_path + 'signal_%s'%savename,signal)
	np.save(output_path + 'noise_%s'%savename,noise)
	np.save(output_path + 'strehl_%s'%savename,strehl)
	np.save(output_path + 'magarr',magarr)
	np.save(output_path + 'exptimes',exptimes)
	np.save(output_path + 'skyinstbg_%s'%savename,skyinstbkg)

	plt.close('all')

def run_all(so,path='./output/centroid_arrs/',field_radius=0):
	"""
	step through all params and run tracking cam for each

	inputs:
	-------	
	so: object
		config object
	path: str
		path to save output
	field_radius: float
		radius of field in arcsec
	
	outputs:
	--------
	none
	"""
	# make folder if doesnt exist
	if not os.path.isdir(path):
		os.makedirs(path)
	else:
		text = ''
		while text!='y' and text!='y':
			text = input ("WARNING, folder %s already exists. Want to still proceed? y/n: "%path)
			if text=='y': pass
			if text=='n': pass
			if (text!='y') and (text!='n'): print('input must be y/n')
	
	# Define parameters to step through
	teffs = np.array([1000,1500,2300,3000,3600,4200,5800])
	track_bands = ['JHgap','Jplus','Hplus','JHplus']
	exptimes = np.array([0.001,0.01,0.1,1,10])
	magarr   = np.arange(4,17,1) # keep in H band?

	# loop through all parameters
	for ii,teff in enumerate(teffs):
		so.stel.teff=teff
		so.track.field_r = field_radius
		cload      = fill_data(so)
		#plot_tools.plot_tracking_bands(so)
		#plt.close()
		for jj, track_band in enumerate(track_bands):
			# reset tracking band
			cload.set_tracking_band_texp(so,track_band,1)
			run_mags_exptimes(so,magarr,exptimes,output_path=path) # this saves files to specified path



if __name__=='__main__':
	#load inputs
	configfile = './hispec_tracking.cfg'
	so         = load_object(configfile)
	cload      = fill_data(so,track_on=True)

	run_all(so,path='./output/centroid_arrs_rad3/',field_radius=3)

	track_requirement = get_track_req(so)
	print("REQ:%0.2f CEN_ERR:%0.2f"%(track_requirement,so.track.centroid_err))



