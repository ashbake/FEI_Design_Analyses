# example tracking camera usage
# plot total throughput for back
# injection paths
# restuls summarized here:
# https://caltech.sharepoint.com/sites/coo/hispec/_layouts/15/doc.aspx?sourcedoc={7916470a-5781-4fa3-be73-b2d2d80c2a55}&action=edit
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

from matplotlib.ticker import (AutoMinorLocator)

plt.ion()
font = {'size'   : 14}
matplotlib.rc('font', **font)
plt.rcParams['font.size'] = '14'
plt.rcParams['font.family'] = 'sans'
plt.rcParams['axes.linewidth'] = '1.3'
fontname = 'DIN Condensed'

path = './data/asahi_roundthree/'


def load_phoenix(teff, x):
	# load phoenix spectrum
	stelpath = '/Users/ashbake/Documents/Research/_DATA/phoenix/'
	stelname = 'lte%s-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'%teff
	f = fits.open(stelpath + stelname)
	spec = f[0].data / (1e8) # ergs/s/cm2/cm to ergs/s/cm2/Angstrom for conversion
	f.close()
	
	wave_file = os.path.join(stelpath + 'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits') #assume wave in same folder
	f = fits.open(wave_file)
	lam = f[0].data # angstroms
	f.close()
	
	spec *=  5.03*10**7 * lam # phot/cm2/s/angstrom
	spec /=np.max(spec) # relative units
	lam/=10 # nm

	# reinterpolate spec onto x
	fspec = interpolate.interp1d(lam,spec,bounds_error=False,fill_value=-1)
	df_star  = fspec(x)

	return df_star

def load_sonora(teff,x):
	"""
	load sonora model file
	
	return subarray 
	
	wav_start, wav_end specified in nm
	
	convert s from erg/cm2/s/Hz to phot/cm2/s/nm using
	https://hea-www.harvard.edu/~pgreen/figs/Conversions.pdf

	wavelenght loaded is microns high to low
	"""
	stelname = '/Users/ashbake/Documents/Research/_DATA/sonora/' + 'sp_t%sg%snc_m0.0' %(str(int(teff)),'316')
	f = np.loadtxt(stelname,skiprows=2)

	lam  = 10000* f[:,0][::-1] #microns to angstroms, needed for conversiosn
	spec = f[:,1][::-1] # erg/cm2/s/Hz
	
	spec *= 3e18/(lam**2)# convert spec to erg/cm2/s/angstrom
	
	conversion_factor = 5.03*10**7 * lam #lam in angstrom here
	spec *= conversion_factor # phot/cm2/s/angstrom
	
	lam,spec =  lam/10.0,spec/np.max(spec)
	
	fspec = interpolate.interp1d(lam,spec,bounds_error=False,fill_value=-1)
	df_star  = fspec(x)

	return df_star

def load_asahi_profiles():
	"""
	"""
	# Asahi dichroic  - JH gap
	filename_jhgap = 'Final transmission model(JH gap dichroic).xlsx'
	xl_jhgap = pd.ExcelFile(path + filename_jhgap)
	first_sheet = xl_jhgap.sheet_names[0]
	df_jhgap = xl_jhgap.parse(first_sheet).values.T

	# Asahi dichroic  - H+JH gap
	filename_hgap = 'Final transmission model(H + JH gap dichroic).xlsx'
	xl_hgap = pd.ExcelFile(path + filename_hgap)
	first_sheet = xl_hgap.sheet_names[0]
	df_hgap = xl_hgap.parse(first_sheet).values.T

	# Asahi dichroic  - J+JH gap
	filename_jgap = 'Final transmission model(J + JH gap dichroic).xlsx'
	xl_jgap = pd.ExcelFile(path + filename_jgap)
	first_sheet = xl_jgap.sheet_names[0]
	df_jgap = xl_jgap.parse(first_sheet).values.T

	# Asahi dichroic  - J+H gap
	filename_jh = 'Final transmission model(J + H dichroic).xlsx'
	xl_jh = pd.ExcelFile(path + filename_jh)
	first_sheet = xl_jh.sheet_names[0]
	df_jh = xl_jh.parse(first_sheet).values.T

	# Asahi dichroic  - CSD
	filename_csd = 'Final transmission model(Channel splitting dichroic).xlsx'
	xl_csd = pd.ExcelFile(path + filename_csd)
	first_sheet = xl_csd.sheet_names[0]
	df_csd = xl_csd.parse(first_sheet).values.T

	# AR coatings
	filename_ar_atc = 'AR_ATC Dichroic.xlsx'
	xl_ar_atc   = pd.ExcelFile(path + filename_ar_atc)
	first_sheet = xl_ar_atc.sheet_names[0]
	df_ar_atc   = xl_ar_atc.parse(first_sheet).values.T

	# AR coatings
	filename_ar_csd = 'AR_Channel splitting dichroic.xlsx'
	xl_ar_csd   = pd.ExcelFile(path + filename_ar_csd)
	first_sheet = xl_ar_csd.sheet_names[0]
	df_ar_csd   = xl_ar_csd.parse(first_sheet).values.T

	filename_block = 'Transmission model(HISPEC FEI blocking filter_IR grade fused silica)-f6_0deg.xlsx'
	xl_ar_block   = pd.ExcelFile(path + filename_block)
	first_sheet   = xl_ar_block.sheet_names[0]
	df_ar_block   = xl_ar_block.parse(first_sheet).values.T
	f = interpolate.interp1d(df_ar_block[0],df_ar_block[1],bounds_error=False,fill_value=-1)
	df_block       = [df_ar_csd[0],f(df_ar_csd[0])]

	return df_jhgap, df_hgap, df_jgap, df_jh, df_csd, df_ar_atc, df_ar_csd, df_block


def compare_throughput_out_of_band(teff,df_atc,df_csd, df_ar_atc, df_ar_csd, df_block):
	"""
	compare out of band flux to in band flux
	for a hotter star with more flux in the blue
	"""
	x = df_atc[0]
	y = (100-df_atc[1])/100 * df_block[1]/100

	# load AO
	xao,yao = np.loadtxt('/Users/ashbake/Documents/Research/Projects/HISPEC/_DATA/throughput/hispec_trackingcamera.csv',delimiter=',').T
	f       = interpolate.interp1d(xao*1000,yao,bounds_error=False,fill_value=-1)
	f_ao    = f(x)

	# load star
	if int(teff) >=2300:
		df_star = load_phoenix(teff,x)
	else:
		df_star = load_sonora(teff,x)

	# plot
	fig, ax = plt.subplots(1,sharex=True,figsize=(10,7))
	ax.plot(x,y,label='ATC Dichroic') # times spectrum
	ax.plot(x,f_ao,label='AO Throughput') # times spectrum
	ax.plot(x,df_star,label='Stellar Flux Density') # times spectrum

	ax.fill_between([980,1100],-0.1,1,color='gray',alpha=0.2)
	ax.fill_between([1170,1327],-0.1,1,color='gray',alpha=0.2)
	ax.fill_between([1490,1780],-0.1,1,color='gray',alpha=0.2)
	ax.fill_between([1990,2460],-0.1,1,color='gray',alpha=0.2)
	ax.set_ylim(-0.1,1)
	ax.set_xlabel('Wavelength [nm]')
	ax.set_ylabel('Throughput')
	ax.legend()

	# integrate between 950 and 1300
	iout  = np.where((x>950) & (x<1300))[0]
	iout2 = np.where((x>1498) & (x<2500))[0]
	iin   = np.where((x>1330) & (x<1488))[0]

	out_flux  = np.trapz(y[iout] * df_star[iout] * f_ao[iout], x[iout])
	out2_flux = np.trapz(y[iout2] * df_star[iout2]* f_ao[iout2], x[iout2])
	in_flux   = np.trapz(y[iin] * df_star[iin]* f_ao[iin], x[iin])

	ratio = (out_flux +out2_flux)/in_flux
	ax.set_title('Star Teff=%sK\nOut to In Band Ratio:%spct'%(teff,np.round(100*ratio,1)))

	# take ratio
	print(ratio)
	return ax

if __name__=='__main__':
	atc_dichroics = {}
	atc_dichroics['jhgap'], atc_dichroics['hgap'], atc_dichroics['jgap'], atc_dichroics['jh'], df_csd, df_ar_atc, df_ar_csd, df_block = load_asahi_profiles()
	band = 'jhgap'
	teff='01000' # 5 digits leading 0
	ax = compare_throughput_out_of_band(teff,atc_dichroics[band],df_csd, df_ar_atc, df_ar_csd, df_block)




