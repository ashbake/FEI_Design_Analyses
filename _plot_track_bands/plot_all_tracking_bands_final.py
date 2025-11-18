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

datapath = './data/'#'/Users/ashbake/Documents/Research/Projects/HISPEC/Throughput_Budget/inputs/quoted_tests/'

def load_asahi_single(filename):
	"""
	Load excel file in asahi format. Assumes contents in first sheet
	
	inputs
	------
		Filename: string
			name of file to open
	returns
	-------
	df - tuple
		x and y in tuple of loaded file
	"""
	xl = pd.ExcelFile(filename)
	first_sheet = xl.sheet_names[0]
	df = xl.parse(first_sheet).values.T

	return df

def load_asahi_data():
	"""
	"""
	path = datapath + 'asahi_final/'

	# Asahi dichroic  - JH gap
	df_jhgap = load_asahi_single(path+'Final transmission model(JH gap dichroic)-3.xlsx')
	df_hgap = load_asahi_single(path+'Final transmission model(H + JH gap dichroic)-2.xlsx')
	df_jgap =  load_asahi_single(path+'Final transmission model(J + JH gap dichroic)-2.xlsx')
	df_jh = load_asahi_single(path+'Final transmission model(J + H dichroic)-2.xlsx')
	df_csd = load_asahi_single(path+'Final transmission model(Channel splitting dichroic).xlsx')

	# Asahi blocking filter
	df_block = load_asahi_single(path+'Transmission model(HISPEC FEI blocking filter_IR grade fused silica)-3.xlsx')

	# Asahi AR coatings
	ar_csd = load_asahi_single(path+'AR_Channel splitting dichroic.xlsx')
	ar_atc = load_asahi_single(path+'AR_ATC Dichroic.xlsx')

	return df_jhgap, df_hgap, df_jgap, df_jh, df_csd, df_block, ar_csd, ar_atc


def plot_asahi_profiles():
	df_jhgap, df_hgap, df_jgap, df_jh, df_csd, df_block, ar_csd, ar_atc  = load_asahi_data() 

	fig, ax = plt.subplots(1,figsize=(9,6))
	ax.plot(df_csd[0],df_csd[1]/100,'k',label='CSD')
	ax.plot(ar_csd[0],ar_csd[1]/100,'k--',lw=.5,label='CSD')
	ax.plot(df_jgap[0], df_jgap[1]/100+4.4,'r',label='J + Gap')
	ax.plot(df_jhgap[0],df_jhgap[1]/100+3.3,'r',label='JH Gap')
	ax.plot(df_hgap[0],df_hgap[1]/100+2.2,'r',label='H + Gap')
	ax.plot(df_jh[0], df_jh[1]/100+1.1,'r',label='J + H')
	ax.plot(ar_atc[0], ar_atc[1]/100+1.1,'r--',lw=.5,label='AR ATC')
	ax.plot(ar_atc[0], ar_atc[1]/100+2.2,'r--',lw=.5)
	ax.plot(ar_atc[0], ar_atc[1]/100+3.3,'r--',lw=.5)
	ax.plot(ar_atc[0], ar_atc[1]/100+4.4,'r--',lw=.5)
	ax.plot(df_block[0], df_block[1]/100+5.5,'c')

	ax.set_xlabel('Wavelength [nm]')
	ax.set_ylabel('Transmission + offset')

	ax.fill_between(so.inst.y,-2.5,10,color='gray',alpha=0.2)
	ax.fill_between(so.inst.J,-2.5,10,color='gray',alpha=0.2)
	ax.fill_between(so.inst.H,-2.5,10,color='gray',alpha=0.2)
	ax.fill_between(so.inst.K,-2.5,10,color='gray',alpha=0.2)

	ax.set_xlim(900,2470)
	ax.set_ylim(-2.2,6.8)

	plt.subplots_adjust(bottom=0.15,left=0.15)

	ax.plot(so.tel.v,so.tel.s -1.1,'b',alpha=0.3)
	ax.plot(so.obs.v,so.obs.s/np.max(so.obs.s) -2.2,'g',alpha=0.3)

	plt.savefig('asahi_profiles_final.png',dpi=350)



if __name__=='__main__':
	#load inputs
	configfile = '../hispec_tracking.cfg'
	so         = load_object(configfile)
	cload      = fill_data(so,track_on=True)

	plot_asahi_profiles()
