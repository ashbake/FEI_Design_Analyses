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


def load_asahi_profiles():
	"""
	"""
	path = datapath + 'asahi_roundone/'

	# Asahi CSD
	filename_csd = 'AsahiTransmissionmodel(Channelsplittingdichroic).xlsx'
	xl_csd = pd.ExcelFile(path + filename_csd)
	first_sheet = xl_csd.sheet_names[0]
	df_csd = xl_csd.parse(first_sheet).values.T

	# Asahi dichroic  - JH gap
	filename_jhgap = 'AsahiTransmissionmodel(JHgapdichroic).xlsx'
	xl_jhgap = pd.ExcelFile(path + filename_jhgap)
	first_sheet = xl_jhgap.sheet_names[0]
	df_jhgap = xl_jhgap.parse(first_sheet).values.T

	# Asahi dichroic  - H+JH gap
	filename_hgap = 'AsahiTransmissionmodel(H+JHgapdichroic).xlsx'
	xl_hgap = pd.ExcelFile(path + filename_hgap)
	first_sheet = xl_hgap.sheet_names[0]
	df_hgap = xl_hgap.parse(first_sheet).values.T

	# Asahi dichroic  - J+JH gap
	filename_jgap = 'AsahiTransmissionmodel(J+JHgapdichroic).xlsx'
	xl_jgap = pd.ExcelFile(path + filename_jgap)
	first_sheet = xl_jgap.sheet_names[0]
	df_jgap = xl_jgap.parse(first_sheet).values.T

	# Asahi dichroic  - J+H gap
	filename_jh = 'AsahiTransmissionmodel(J+Hdichroic).xlsx'
	xl_jh = pd.ExcelFile(path + filename_jh)
	first_sheet = xl_jh.sheet_names[0]
	df_jh = xl_jh.parse(first_sheet).values.T

	return df_csd, df_jhgap, df_hgap, df_jgap, df_jh

def load_asahi_j_h():
	filename_h = 'AsahiTransmissionmodel(Hbandfilter).xlsx'
	filename_j = 'AsahiTransmissionmodel(Jbandfilter).xlsx'
	path = datapath + 'asahi_roundone/'

	# Asahi dichroic  - H+JH gap
	xl_h = pd.ExcelFile(path + filename_h)
	first_sheet = xl_h.sheet_names[0]
	df_h = xl_h.parse(first_sheet).values.T

	# Asahi dichroic  - J+JH gap
	xl = pd.ExcelFile(path + filename_j)
	first_sheet = xl.sheet_names[0]
	df_j = xl.parse(first_sheet).values.T

def load_asahi_round2():
	"""
	"""
	path = datapath + 'asahi_roundtwo/'

	# Asahi dichroic  - JH gap
	filename_jhgap = 'Modified transmission model(JH gap dichroic).xlsx'
	xl_jhgap = pd.ExcelFile(path + filename_jhgap)
	first_sheet = xl_jhgap.sheet_names[0]
	df_jhgap = xl_jhgap.parse(first_sheet).values.T

	# Asahi dichroic  - H+JH gap
	filename_hgap = 'Modified transmission model(H + JH gap dichroic).xlsx'
	xl_hgap = pd.ExcelFile(path + filename_hgap)
	first_sheet = xl_hgap.sheet_names[0]
	df_hgap = xl_hgap.parse(first_sheet).values.T

	# Asahi dichroic  - J+JH gap
	filename_jgap = 'Modified transmission model(J + JH gap dichroic).xlsx'
	xl_jgap = pd.ExcelFile(path + filename_jgap)
	first_sheet = xl_jgap.sheet_names[0]
	df_jgap = xl_jgap.parse(first_sheet).values.T

	# Asahi dichroic  - J+H gap
	filename_jh = 'Modified transmission model(J + H dichroic).xlsx'
	xl_jh = pd.ExcelFile(path + filename_jh)
	first_sheet = xl_jh.sheet_names[0]
	df_jh = xl_jh.parse(first_sheet).values.T

	return df_jhgap, df_hgap, df_jgap, df_jh


def load_asahi_round3():
	"""
	"""
	path = datapath + 'asahi_roundthree/'

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

	return df_jhgap, df_hgap, df_jgap, df_jh, df_csd

def load_materion():
	path = datapath + 'materion/'

	filename = 'materionCSDJY_HKDBS.csv'
	lcsd,_,Tcsd=np.loadtxt(path + filename,delimiter=',',skiprows=3).T

	filename = 'materionJHGapDBS.csv'
	ljhgap, _,Tjhgap = np.loadtxt(path + filename,delimiter=',',skiprows=3).T
	
	filename = 'materionJ+JHGapDBS.csv'
	ljgap, _,Tjgap = np.loadtxt(path + filename,delimiter=',',skiprows=3).T
	
	filename = 'materionH+JHGapDBS.csv'
	lhgap, _,Thgap = np.loadtxt(path + filename,delimiter=',',skiprows=3).T

	filename = 'materionJ+HDBS.csv'
	ljh, _,Tjh = np.loadtxt(path + filename,delimiter=',',skiprows=3).T

	filename='materionFEIBlockingFilter.csv'
	lbl, _,Tbl= np.loadtxt(path + filename,delimiter=',',skiprows=3).T

	# these aren't final optimized files and
	# generally look worse than asahi ones
	# so i never got around to making a plotting script 
	return [lcsd,Tcsd], [ljhgap,Tjhgap],[lhgap,Thgap], [ljgap,Tjgap], [ljh,Tjh], [lbl,Tbl]

def plot_asahi_profiles():
	df_csd, df_jhgap, df_hgap, df_jgap, df_jh = load_asahi_profiles()
	df_jhgap2, df_hgap2, df_jgap2, df_jh2     = load_asahi_round2()
	df_csd3, df_jhgap3, df_hgap3, df_jgap3, df_jh3, df_bl3 = load_materion()
	df_jhgap3, df_hgap3, df_jgap3, df_jh3, df_csd3  = load_asahi_round3() 

	fig, ax = plt.subplots(1,figsize=(7,5))
	ax.plot(df_csd[0],df_csd[1]/100,'k',label='CSD')
	ax.plot(df_jgap[0], df_jgap[1]/100+4.4,'r',label='J + Gap')
	ax.plot(df_jhgap[0],df_jhgap[1]/100+3.3,'r',label='JH Gap')
	ax.plot(df_hgap[0],df_hgap[1]/100+2.2,'r',label='H + Gap')
	ax.plot(df_jh[0], df_jh[1]/100+1.1,'r',label='J + H')

	ax.plot(df_jgap2[0], df_jgap2[1]/100+4.4,'orange',label='J + Gap')
	ax.plot(df_jhgap2[0],df_jhgap2[1]/100+3.3,'orange',label='JH Gap')
	ax.plot(df_hgap2[0],df_hgap2[1]/100+2.2,'orange',label='H + Gap')
	ax.plot(df_jh2[0], df_jh2[1]/100+1.1,'orange',label='J + H')

	ax.plot(df_csd3[0],df_csd3[1]/100,'g',label='CSD')
	ax.plot(df_jgap3[0], df_jgap3[1]/100+4.4,'g',label='J + Gap')
	ax.plot(df_jhgap3[0],df_jhgap3[1]/100+3.3,'g',label='JH Gap')
	ax.plot(df_hgap3[0],df_hgap3[1]/100+2.2,'g',label='H + Gap')
	ax.plot(df_jh3[0], df_jh3[1]/100+1.1,'g',label='J + H')
	
	#ax.plot(df_j[0], df_j[1]/100+2.2,'b',label='J Filter')
	#ax.plot(df_h[0], df_h[1]/100+1.1,'b',label='H Filter')
	#ax.plot(df_jgap[0], 1-df_jgap[1]/100+2.2,'r--',label='J + Gap')
	#ax.plot(df_hgap[0],1-df_hgap[1]/100+1.1,'r--',label='H + Gap')

	ax.set_xlabel('Wavelength [nm]')
	ax.set_ylabel('Transmission + offset')

	ax.fill_between(so.inst.y,-2.5,10,color='gray',alpha=0.2)
	ax.fill_between(so.inst.J,-2.5,10,color='gray',alpha=0.2)
	ax.fill_between(so.inst.H,-2.5,10,color='gray',alpha=0.2)
	ax.fill_between(so.inst.K,-2.5,10,color='gray',alpha=0.2)

	ax.set_xlim(950,2470)
	ax.set_ylim(-2.2,8.5)

	plt.subplots_adjust(bottom=0.15,left=0.15)

	ax.plot(so.tel.v,so.tel.s -1.1,'b',alpha=0.3)
	ax.plot(so.obs.v,so.obs.s/np.max(so.obs.s) -2.2,'g',alpha=0.3)

	plt.savefig('asahi_profiles.png')

def plot_total_throughput():
	"""
	"""
	# load throughput of tracking camera without filters to assess cutoff
	xao,yao = np.loadtxt('/Users/ashbake/Documents/Research/Projects/HISPEC/_DATA/throughput/hispec_trackingcamera.csv',delimiter=',').T
	df_csd, df_jhgap, df_hgap, df_jgap, df_jh = load_asahi_profiles()
	
	f= interpolate.interp1d(xao*1000,yao,bounds_error=False,fill_value=-1)

	fig, ax = plt.subplots(1,figsize=(7,5))
	ax.plot(df_jhgap[0], 1- df_jhgap[1]/100,label='Dichroic only')
	ax.plot(xao*1000,yao,label='Pre Dichroic Tracking\n Camera Throughput')
	ax.plot(df_jhgap[0], (1 - df_jhgap[1]/100)* f(df_jhgap[0]),lw=2,label='Total Throughput')
	ax.set_xlabel('Wavelength [nm]')
	ax.set_ylabel('Transmission')
	ax.legend(fontsize=10,loc=1)
	ax.set_ylim(-0.1,1.1)
	plt.savefig('total_throughput_jhgap.png')


	fig, ax = plt.subplots(1,figsize=(7,5))
	ax.plot(df_hgap[0], 1- df_hgap[1]/100,label='Dichroic only')
	ax.plot(xao*1000,yao,label='Pre Dichroic Tracking\n Camera Throughput')
	ax.plot(df_hgap[0], (1 - df_hgap[1]/100)* f(df_hgap[0]),lw=2,label='Total Throughput')
	ax.set_xlabel('Wavelength [nm]')
	ax.set_ylabel('Transmission')
	ax.legend(fontsize=10,loc=1)
	ax.set_ylim(-0.1,1.1)
	plt.savefig('total_throughput_hgap.png')


	fig, ax = plt.subplots(1,figsize=(7,5))
	ax.plot(df_jgap[0], 1- df_jgap[1]/100,label='Dichroic only')
	ax.plot(xao*1000,yao,label='Pre Dichroic Tracking\n Camera Throughput')
	ax.plot(df_jgap[0], (1 - df_jgap[1]/100)* f(df_jgap[0]),lw=2,label='Total Throughput')
	ax.set_xlabel('Wavelength [nm]')
	ax.set_ylabel('Transmission')
	ax.legend(fontsize=10,loc=1)
	ax.set_ylim(-0.1,1.1)
	plt.savefig('total_throughput_jgap.png')

	fig, ax = plt.subplots(1,figsize=(7,5))
	ax.plot(df_jh[0], 1- df_jh[1]/100,label='Dichroic only')
	ax.plot(xao*1000,yao,label='Pre Dichroic Tracking\n Camera Throughput')
	ax.plot(df_jh[0], (1 - df_jh[1]/100)* f(df_jh[0]),lw=2,label='Total Throughput')
	ax.set_xlabel('Wavelength [nm]')
	ax.set_ylabel('Transmission')
	ax.legend(fontsize=10,loc=1)
	ax.set_ylim(-0.1,1.1)
	plt.savefig('total_throughput_jh.png')


def plot_asahi_profiles_final():
	df_csd, df_jhgap, df_hgap, df_jgap, df_jh = load_asahi_profiles()
	#df_jhgap2, df_hgap2, df_jgap2, df_jh2     = load_asahi_round2()
	#df_csd3, df_jhgap3, df_hgap3, df_jgap3, df_jh3, df_bl3 = load_materion()
	df_jhgap3, df_hgap3, df_jgap3, df_jh3, df_csd3  = load_asahi_round3() 

	fig, ax = plt.subplots(1,figsize=(7,5))
	ax.plot(df_csd[0],df_csd[1]/100,'k',label='CSD')

	ax.plot(df_jgap3[0], df_jgap3[1]/100+4.4,'g',label='J + Gap')
	ax.plot(df_jhgap3[0],df_jhgap3[1]/100+3.3,'g',label='JH Gap')
	ax.plot(df_hgap3[0],df_hgap3[1]/100+2.2,'g',label='H + Gap')
	ax.plot(df_jh3[0], df_jh3[1]/100+1.1,'g',label='J + H')
	
	#ax.plot(df_j[0], df_j[1]/100+2.2,'b',label='J Filter')
	#ax.plot(df_h[0], df_h[1]/100+1.1,'b',label='H Filter')
	#ax.plot(df_jgap[0], 1-df_jgap[1]/100+2.2,'r--',label='J + Gap')
	#ax.plot(df_hgap[0],1-df_hgap[1]/100+1.1,'r--',label='H + Gap')

	ax.set_xlabel('Wavelength [nm]')
	ax.set_ylabel('Transmission + offset')

	ax.fill_between(so.inst.y,-2.5,10,color='gray',alpha=0.2)
	ax.fill_between(so.inst.J,-2.5,10,color='gray',alpha=0.2)
	ax.fill_between(so.inst.H,-2.5,10,color='gray',alpha=0.2)
	ax.fill_between(so.inst.K,-2.5,10,color='gray',alpha=0.2)

	ax.set_xlim(950,2470)
	ax.set_ylim(-2.2,5.5)

	plt.subplots_adjust(bottom=0.15,left=0.15)

	ax.plot(so.tel.v,so.tel.s -1.1,'b',alpha=0.3)
	ax.plot(so.obs.v,so.obs.s/np.max(so.obs.s) -2.2,'g',alpha=0.3)

	# Spec sheet profiles
	# jh gap
	# ref 1347.7 nm to 1453.9 
	#tran 980 nm to 1312.7 nm and 1505.7 to 2500 nm
	ax.plot([1347.7, 1453.9],[3.3,3.3],'k',lw=2)
	ax.plot([980, 1312.7],[3.3,3.3],'r',lw=2)
	ax.plot([1505.7, 2500],[3.3,3.3],'r',lw=2)

	# J+jhgap  1167.9 nm - 1447.9 , 980 to 1132.1nm and 1519.4 nm - 2500 nm
	ax.plot([1167.9, 1447.9],[4.4,4.4],'k',lw=2)
	ax.plot([980, 1132.1],[4.4,4.4],'r',lw=2)
	ax.plot([1519.4 , 2500],[4.4,4.4],'r',lw=2)

	# H+jhgap  1355.4 nm - 1783.5 ,1305.9 nm and 1856.8
	ax.plot([1355.4, 1783.5 ],[2.2,2.2],'k',lw=2)
	ax.plot([980, 1305.9],[2.2,2.2],'r',lw=2)
	ax.plot([1856.8 , 2500],[2.2,2.2],'r',lw=2)

	# H+J  
	ax.plot([1165.3,1782.8],[1.1,1.1],'k',lw=2)
	ax.plot([980,1132.8],[1.1,1.1],'r',lw=2)
	ax.plot([1855.7 , 2500],[1.1,1.1],'r',lw=2)

	plt.savefig('asahi_profiles_final.png')




if __name__=='__main__':
	#load inputs
	configfile = '../hispec_tracking.cfg'
	so         = load_object(configfile)
	cload      = fill_data(so,track_on=True)

	plot_asahi_profiles()
