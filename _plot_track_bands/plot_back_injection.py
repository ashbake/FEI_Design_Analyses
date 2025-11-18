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

datapath = './data/asahi_roundthree/'



def load_asahi_profiles(path):
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


def plot_red_backinjection(atc_dichroic,df_csd, df_ar_atc, df_ar_csd, df_block,savename='test'):
	"""
	plto back injection
	"""	
	# CSD is in transmission
	# also initate total throughput seen by ATC
	csd_t = (100 - df_ar_csd[1])/100 * df_csd[1].copy()/100

	# dichroic is in reflection - passes through AR coating twice
	dichroic_r = (100-df_ar_atc[1])**2 * (100-atc_dichroic[1])/100/100/100

	# assume retroreflector is achromatic

	# dichroic is in transmission, AR in 
	dichroic_t = (100-df_ar_atc[1]) * atc_dichroic[1]/100/100

	# filters - assume no filter

	# cold window
	# pass

	# blocking filter in transmission
	block_t = df_block[1]/100

	# total
	total = csd_t * dichroic_r * dichroic_t * block_t

	# get ghosts
	ghost_atc = get_red_ghost_ATC(atc_dichroic,df_csd, df_ar_atc, df_ar_csd, df_block)
	ghost_csd = get_red_ghost_CSD(atc_dichroic,df_csd, df_ar_atc, df_ar_csd, df_block)

	fig, axs = plt.subplots(2,figsize=(10,8),sharex=True)
	ax,ax2 = axs
	ax.plot(df_csd[0], csd_t,label='CSD trans.')
	ax.plot(df_csd[0], dichroic_r,label='ATC refl.')
	ax.plot(df_csd[0], dichroic_t,label='ATC trans.')
	ax.plot(df_csd[0], block_t,label='Cold Filter')
	ax.plot(df_csd[0], df_ar_csd[1]/100,label='CSD AR coating')

	ax.plot(df_csd[0], total,'k',label='Red arm')

	ax.fill_between([980,1100],-0.1,1,color='gray',alpha=0.2)
	ax.fill_between([1170,1327],-0.1,1,color='gray',alpha=0.2)
	ax.fill_between([1490,1780],-0.1,1,color='gray',alpha=0.2)
	ax.fill_between([1990,2460],-0.1,1,color='gray',alpha=0.2)
	
	ax2.plot(df_csd[0], total,'k')
	
	# plot ghosts
	ax2.plot(df_csd[0], ghost_atc,'gray',label='ATC ghost')
	ax2.plot(df_csd[0], ghost_csd,'gray',label='CSD ghost', ls='--')

	ax2.fill_between([980,1100],-0.1,1,color='gray',alpha=0.2)
	ax2.fill_between([1170,1327],-0.1,1,color='gray',alpha=0.2)
	ax2.fill_between([1490,1780],-0.1,1,color='gray',alpha=0.2)
	ax2.fill_between([1990,2460],-0.1,1,color='gray',alpha=0.2)

	ax2.legend(fontsize=10,loc=7)
	ax.legend(fontsize=10,loc=7)
	ax2.set_xlabel('Wavelength [nm]')
	ax2.set_ylabel('Transmission')
	ax.set_ylabel('Transmission')
	ax.set_title('Red Back Injection Throughput\n%s'%savename)
	ax2.set_title('Total Transmission')
	ax.grid()
	ax2.grid()
	plt.savefig('plots_back_injection/Red_BackInjection_trans_%s.png'%savename)

	############
	#
	#
	# plot ratio	
	############
	fig, axs = plt.subplots(2,sharex=True,figsize=(10,7))
	ax,ax2=axs
	ax.plot(df_csd[0], total,'k',label='Red arm total')
	ax.plot(df_csd[0], ghost_atc,'gray',label='Ghost ATC')
	ax.plot(df_csd[0], ghost_csd,'gray',ls='--',label='Ghost CSD')

	ax2.fill_between(df_csd[0], total/ghost_atc,facecolor='steelblue', alpha=0.5,label='ATC Ghost ratio')
	#ax2twin = ax2.twinx()
	ax2.semilogy(df_csd[0], total/ghost_csd,'gray', label='CSD Ghost ratio')
	
	ax.fill_between([980,1100],-0.1,1,color='gray',alpha=0.2)
	ax.fill_between([1170,1327],-0.1,1,color='gray',alpha=0.2)
	ax.fill_between([1490,1780],-0.1,1,color='gray',alpha=0.2)
	ax.fill_between([1990,2460],-0.1,1,color='gray',alpha=0.2)
	ax2.fill_between([980,1100],-0.1,80,color='gray',alpha=0.2)
	ax2.fill_between([1170,1327],-0.1,80,color='gray',alpha=0.2)
	ax2.fill_between([1490,1780],-0.1,80,color='gray',alpha=0.2)
	ax2.fill_between([1990,2460],-0.1,80,color='gray',alpha=0.2)

	ax.set_xlim(950,1800)
	ax.set_ylim(-.1,0.3)


	aerodiodes = [1410, 1430, 1450, 1470, 1480, 1490,1510,1525,1530]
	for aerodiode in aerodiodes:
		ax2.axvline(aerodiode,color='m',lw=2,ls='--')

	ax.legend(fontsize=10,loc='best')	
	ax2.legend(fontsize=10,loc='best')
	ax.set_title('Red Ghosts (%s)'%savename)
	ax2.set_title('Red Beacon to ATC/CSD ghost flux')
	ax2.set_xlabel('Wavelength [nm]')
	ax2.set_ylabel('Flux Ratio')
	ax.set_ylabel('Transmission')
	plt.savefig('plots_back_injection/Red_BackInjection_ghostratio_%s.png'%savename)

	return total

def get_red_ghost_ATC(atc_dichroic,df_csd, df_ar_atc, df_ar_csd, df_block):
	"""
	plot ghosts from back injection
	ATC ghost - this is if the beam reflects off the ATC AR coating
	"""	
	# CSD is in transmission, AR in transmission
	# also initate total throughput seen by ATC
	csd_t = (100 - df_ar_csd[1])/100 * df_csd[1].copy()/100

	# dichroic is in reflection - ghosts reflects off AR coating
	dichroic_r = (df_ar_atc[1])/100

	# assume retroreflector is achromatic

	# dichroic is in transmission including AR coating
	dichroic_t = (100-df_ar_atc[1]) * atc_dichroic[1]/100/100

	# filters - assume no filter

	# cold window
	# pass

	# blocking filter in transmission
	block_t = df_block[1]/100

	# total
	total_ghost = csd_t * dichroic_r * dichroic_t * block_t

	return total_ghost


def get_red_ghost_CSD(atc_dichroic,df_csd, df_ar_atc, df_ar_csd, df_block):
	"""
	plot ghosts from back injection
	CSD ghost - this is if the beam does a back and forth in the CSD
	"""	
	# AR transmission, CSD reflection, then AR reflect, CSD transmission
	# also initate total throughput seen by ATC
	csd_1 = (100 - df_ar_csd[1])/100 * (100 - df_csd[1])/100
	csd_2 = (df_ar_csd[1])/100 * (df_csd[1])/100

	# dichroic is in reflection - passes through AR coating twice
	dichroic_r = (100-df_ar_atc[1])**2 * (100-atc_dichroic[1])/100/100/100

	# assume retroreflector is achromatic

	# dichroic is in transmission, AR in 
	dichroic_t = (100-df_ar_atc[1]) * atc_dichroic[1]/100/100

	# filters - assume no filter

	# cold window
	# pass

	# blocking filter in transmission
	block_t = df_block[1]/100

	# total
	total_ghost = csd_1 * csd_2 * dichroic_r * dichroic_t * block_t

	return total_ghost


def plot_blue_backinjection(atc_dichroic,df_csd, df_ar_atc, df_ar_csd, df_block,savename='test'):
	"""

	"""	
	# CSD is in reflection
	# blue doesnt see AR coating
	csd_r = (1 - df_csd[1].copy()/100)

	# dichroic is in reflection - passes through AR coating twice
	dichroic_r = (100-df_ar_atc[1])**2 * (100-atc_dichroic[1])/100/100/100

	# assume retroreflector is achromatic

	# dichroic is in transmission
	dichroic_t = (100-df_ar_atc[1]) * atc_dichroic[1]/100/100

	# filters - assume no filter

	# cold window
	# pass

	# blocking filter in transmission
	block_t = df_block[1]/100

	# total
	total = csd_r * dichroic_r * dichroic_t * block_t
	
	# load ghost totals
	total_ghost_atc = get_blue_ATCghost(atc_dichroic,df_csd, df_ar_atc, df_ar_csd, df_block)
	total_ghost_csd = get_blue_CSDghost(atc_dichroic,df_csd, df_ar_atc, df_ar_csd, df_block)

	# 
	fig, axs = plt.subplots(2,sharex=True,figsize=(10,7))
	ax,ax2=axs
	ax.plot(df_csd[0], csd_r,label='CSD refl.')
	ax.plot(df_csd[0], dichroic_r,label='ATC refl.')
	ax.plot(df_csd[0], dichroic_t,label='ATC trans.')
	ax.plot(df_csd[0], block_t,label='Cold Filter')
	ax.plot(df_csd[0], total,'k',label='Blue arm total')
	ax.plot(df_csd[0],(df_ar_atc[1])/100,label='AR Reflection')

	ax2.plot(df_csd[0], total,'k',label='Blue arm total')
	ax2.plot(df_csd[0], total_ghost_atc,'gray',label='Ghost ATC')
	ax2.plot(df_csd[0], total_ghost_csd,'gray',ls='--',label='Ghost CSD')
	
	ax.legend(loc=4,fontsize=12)
	ax.fill_between([980,1100],-0.1,1,color='gray',alpha=0.2)
	ax.fill_between([1170,1327],-0.1,1,color='gray',alpha=0.2)
	ax.fill_between([1490,1780],-0.1,1,color='gray',alpha=0.2)
	ax.fill_between([1990,2460],-0.1,1,color='gray',alpha=0.2)

	ax2.set_xlabel('Wavelength [nm]')
	ax2.set_ylabel('Transmission')
	ax2.legend(fontsize=10)

	ax2.fill_between([980,1100],-0.1,1,color='gray',alpha=0.2)
	ax2.fill_between([1170,1327],-0.1,1,color='gray',alpha=0.2)
	ax2.fill_between([1490,1780],-0.1,1,color='gray',alpha=0.2)
	ax2.fill_between([1990,2460],-0.1,1,color='gray',alpha=0.2)

	ax.set_ylabel('Transmission')
	ax.grid()
	ax2.grid()
	ax2.set_title('Total Transmission')

	ax.set_title('Blue Back Injection Throughput\n%s'%savename)
	plt.savefig('plots_back_injection/Blue_BackInjection_trans_%s.png'%savename)

	# plot ghost ratio
	fig, axs = plt.subplots(2,sharex=True,figsize=(10,7))
	ax,ax2=axs
	ax.plot(df_csd[0], total,'k',label='Blue arm total')
	ax.plot(df_csd[0], total_ghost_atc,'gray',label='Ghost ATC')
	ax.plot(df_csd[0], total_ghost_csd,'gray',ls='--',label='Ghost CSD')

	ax2.fill_between(df_csd[0], total/total_ghost_atc,facecolor='steelblue', alpha=0.5,label='ATC Ghost ratio')
	ax2.semilogy(df_csd[0], total/total_ghost_csd,'gray', label='CSD Ghost ratio')
	ax.fill_between([980,1100],-0.1,1,color='gray',alpha=0.2)
	ax.fill_between([1170,1327],-0.1,1,color='gray',alpha=0.2)
	ax.fill_between([1490,1780],-0.1,1,color='gray',alpha=0.2)
	ax.fill_between([1990,2460],-0.1,1,color='gray',alpha=0.2)
	ax2.fill_between([980,1100],-0.1,80,color='gray',alpha=0.2)
	ax2.fill_between([1170,1327],-0.1,80,color='gray',alpha=0.2)
	ax2.fill_between([1490,1780],-0.1,80,color='gray',alpha=0.2)
	ax2.fill_between([1990,2460],-0.1,80,color='gray',alpha=0.2)
	aerodiodes = [1270, 1290, 1310, 1330, 1350, 1370, 1390]
	for aerodiode in aerodiodes:
		ax2.axvline(aerodiode,color='m',lw=2,ls='--')

	ax.set_xlim(950,2800)
	ax.set_ylim(-.1,0.4)

	ax.legend(fontsize=10,loc='best')
	ax2.legend(fontsize=10,loc='best')
	ax.set_title('Blue Ghosts(%s)'%savename)
	ax2.set_title('Blue Beacon to ATC/CSD ghost ratio')
	ax2.set_xlabel('Wavelength [nm]')
	ax2.set_ylabel('Flux Ratio')
	ax.set_ylabel('Transmission')
	plt.savefig('plots_back_injection/Blue_BackInjection_ghostratio_%s.png'%savename)

	return total

def get_blue_ATCghost(atc_dichroic,df_csd, df_ar_atc, df_ar_csd, df_block):
	"""
	ghost from reflecting off ATC
	"""
	# CSD is in reflection
	# also initate total throughput seen by ATC
	csd_r = (1 - df_csd[1].copy()/100)

	# ghost reflects off AR coating
	dichroic_r = df_ar_atc[1] /100

	# assume retroreflector is achromatic

	# dichroic is in transmission, one pass through AR coating
	dichroic_t = (100-df_ar_atc[1]) * atc_dichroic[1]/100/100

	# filters - assume no filter

	# cold window
	# pass

	# blocking filter in transmission
	block_t = df_block[1]/100

	# total
	total_ghost = csd_r * dichroic_r * dichroic_t * block_t
	
	return total_ghost


def get_blue_CSDghost(atc_dichroic,df_csd, df_ar_atc, df_ar_csd, df_block):
	"""

	"""

	# ghost transmits through CSD, then reflects off AR coating
	# then passes 
	csd_r =  df_csd[1].copy()/100 * (df_ar_csd[1])/100 * (1 - df_csd[1].copy()/100)

	# dichroic is in reflection - passes through AR coating twice
	dichroic_r = (100-df_ar_atc[1])**2 * (100-atc_dichroic[1])/100/100/100

	# assume retroreflector is achromatic

	# dichroic is in transmission
	dichroic_t = (100-df_ar_atc[1]) * atc_dichroic[1]/100/100

	# filters - assume no filter

	# cold window
	# pass

	# blocking filter in transmission
	block_t = df_block[1]/100

	# total
	total_ghost2 = csd_r * dichroic_r * dichroic_t * block_t
	
	return total_ghost2


def compare_throughput_out_of_band():
	"""
	compare out of band flux to in band flux
	for a hotter star with more flux in the blue
	"""
	df_jhgap, _, _, _, _, _, _, df_block = load_asahi_profiles()

	x = df_jhgap[0]
	y = (100-df_jhgap[1])/100 * df_block[1]/100


	# load phoenix spectrum
	from astropy.io import fits
	stelpath = '/Users/ashbake/Documents/Research/_DATA/phoenix/'
	stelname = 'lte06600-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
	f = fits.open(stelpath + stelname)
	spec = f[0].data / (1e8) # ergs/s/cm2/cm to ergs/s/cm2/Angstrom for conversion
	f.close()
	
	wave_file = os.path.join(stelpath + 'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits') #assume wave in same folder
	f = fits.open(wave_file)
	lam = f[0].data # angstroms
	f.close()
	
	spec *=  5.03*10**7 * lam # phot/cm2/s/angstrom
	lam/=10 # nm

	# reinterpolate spec onto x

	# plot
	fig, ax = plt.subplots(1,sharex=True,figsize=(10,7))
	ax.plot(x,y) # times spectrum

	ax.fill_between([980,1100],-0.1,1,color='gray',alpha=0.2)
	ax.fill_between([1170,1327],-0.1,1,color='gray',alpha=0.2)
	ax.fill_between([1490,1780],-0.1,1,color='gray',alpha=0.2)
	ax.fill_between([1990,2460],-0.1,1,color='gray',alpha=0.2)

	# integrate between 950 and 1300

	# integrate between 1550 and 1800

	# take ratio
	# try diff stel temp


if __name__=='__main__':
	atc_dichroics = {}
	atc_dichroics['jhgap'], atc_dichroics['hgap'], atc_dichroics['jgap'], atc_dichroics['jh'], df_csd, df_ar_atc, df_ar_csd, df_block = load_asahi_profiles(datapath)

	# SAVE DATA TO DICTIONARY
	atcfilts = ['jhgap','jgap','hgap','jh']
	totals = {}
	totals['red'] = {}
	totals['blue'] = {}
	for atcfilt in atcfilts:
		totals['red'][atcfilt]  = plot_red_backinjection(atc_dichroics[atcfilt],df_csd, df_ar_atc, df_ar_csd, df_block,savename=atcfilt)
		totals['blue'][atcfilt] = plot_blue_backinjection(atc_dichroics[atcfilt],df_csd, df_ar_atc, df_ar_csd, df_block,savename=atcfilt)

	# PLOT ALL THROUGHPUTS ON ONE GRAPH
	fig, axs = plt.subplots(2,1,figsize=(7,6))
	for atcfilt in atcfilts:
		p=axs[1].plot(atc_dichroics[atcfilt][0], totals['red'][atcfilt], label=atcfilt)
		axs[0].plot(atc_dichroics[atcfilt][0], totals['blue'][atcfilt], label=atcfilt,c=p[0].get_color())
		np.savetxt('plots_back_injection/Red_Backinjection_throughputs_%s.txt'%atcfilt,
			 		np.vstack((atc_dichroics[atcfilt][0], totals['red'][atcfilt])).T)
		np.savetxt('plots_back_injection/Blue_Backinjection_throughputs_%s.txt'%atcfilt,
			 		np.vstack((atc_dichroics[atcfilt][0], totals['blue'][atcfilt])).T)

	axs[0].legend()
	axs[1].set_title('Red Back Injection Throughput')
	axs[0].set_title('Blue Back Injection Throughput')
	axs[1].set_xlabel('Wavelength [nm]')
	axs[0].set_ylabel('Throughput')
	axs[1].set_ylabel('Throughput')

	axs[0].set_xlim(1330, 1485)
	axs[1].set_xlim(1400, 1470)
	plt.subplots_adjust(hspace=0.4)
	axs[0].set_yscale("log")  
	axs[1].set_yscale("log")

	#axs[0].set_ylim(1-5)
	