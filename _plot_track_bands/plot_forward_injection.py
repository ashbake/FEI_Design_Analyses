# example tracking camera usage
import sys, os
import matplotlib
import numpy as np
import matplotlib.pylab as plt
from scipy import interpolate
import pandas as pd
from datetime import date

from matplotlib.ticker import (MultipleLocator)


font = {'size'   : 14}
matplotlib.rc('font', **font)

sys.path.append('/Users/ashbake/Documents/Research/ToolBox/specsim/')

from specsim import throughput_tools
#plt.ion()

from matplotlib.ticker import (AutoMinorLocator)

plt.ion()
font = {'size'   : 14}
matplotlib.rc('font', **font)
plt.rcParams['font.size'] = '14'
plt.rcParams['font.family'] = 'sans'
plt.rcParams['axes.linewidth'] = '1.3'
fontname = 'DIN Condensed'

datapath = './data/'#'/Users/ashbake/Documents/Research/Projects/HISPEC/Throughput_Budget/inputs/quoted_tests/'
datapath = '/Users/ashbake/Documents/Research/Projects/HISPEC/Tests/FEI_TrackingCamera/_plot_track_bands/data/'

def load_asahi_profiles():
	"""
	round 3 is the final ones
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

	#blocking
	filename_block = 'Transmission model(HISPEC FEI blocking filter_IR grade fused silica)-f6_0deg.xlsx'
	xl_block   = pd.ExcelFile(path + filename_block)
	first_sheet   = xl_block.sheet_names[0]
	df_block   = xl_block.parse(first_sheet).values.T
	f = interpolate.interp1d(df_block[0],df_block[1],bounds_error=False,fill_value=-1)
	df_block       = [df_csd[0],f(df_csd[0])]

	return df_jhgap, df_hgap, df_jgap, df_jh, df_csd, df_block

def plot_asahi_profiles():
	df_jhgap3, df_hgap3, df_jgap3, df_jh3, df_csd3  = load_asahi_round3() 

	fig, ax = plt.subplots(1,figsize=(7,5))

	ax.plot(df_csd3[0],df_csd3[1]/100,'g',label='CSD')
	ax.plot(df_jgap3[0], df_jgap3[1]/100+4.4,'g',label='J + Gap')
	ax.plot(df_jhgap3[0],df_jhgap3[1]/100+3.3,'g',label='JH Gap')
	ax.plot(df_hgap3[0],df_hgap3[1]/100+2.2,'g',label='H + Gap')
	ax.plot(df_jh3[0], df_jh3[1]/100+1.1,'g',label='J + H')
	

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

	plt.savefig('plots_tracking_bands/final_filter_profiles.png')

def plot_total_throughput(ploton=False):
	"""
	inputs:
	-------
		ploton: Bool
			will plot all plots if true, otherwise returns throughput data
	"""
	savepath = './plots_forward_throughput/'
	
	# load throughput of tracking camera without filters to assess cutoff
	df_jhgap, df_hgap, df_jgap, df_jh, df_csd, df_block = load_asahi_profiles()
	
	# load throughputs to spectrograph, no PL, no piaa boost
	transmission_path = '/Users/ashbake/Documents/Research/Projects/HISPEC/_DATA/throughput/hispec_subsystems/'
	coupling,_   = throughput_tools.pick_coupling_rounded(transmission_path, df_jhgap[0], 0, 0.5,lo_wfe=0,pl_on=0,piaa_boost=1)
	_, data = throughput_tools.get_base_throughput(df_jhgap[0],ploton=False,datapath='/Users/ashbake/Documents/Research/Projects/HISPEC/_DATA/throughput/hispec_subsystems/')

	# total atc is 1-trans times ao,fei cam throughput times blocking filter
	total_atc_jhgap = (1 - df_jhgap[1]/100)* data['red']['feicom'] *df_block[1]/100
	total_atc_hgap  = (1 - df_hgap[1]/100) * data['red']['feicom'] *df_block[1]/100
	total_atc_jgap  = (1 - df_jgap[1]/100) * data['red']['feicom'] *df_block[1]/100
	total_atc_jh    = (1 - df_jh[1]/100)   * data['red']['feicom'] *df_block[1]/100

	# total csd is FEI optics times ATC dichroic times CSD times coupling, fib, and spectrograph
	total_csd_red_jhgap  = df_jhgap[1]/100* data['red']['feicom'] *df_csd[1]/100 * data['red']['feired'] #* coupling* data['red']['rspec']* data['red']['fibred']
	total_csd_blue_jhgap = df_jhgap[1]/100* data['red']['feicom'] *(1-df_csd[1]/100)* data['blue']['feiblue'] #* coupling* data['blue']['bspec']* data['blue']['fibblue']
	total_csd_red_jgap   = df_jgap[1]/100* data['red']['feicom'] *df_csd[1]/100 * data['red']['feired']# * coupling* data['red']['rspec']* data['red']['fibred']
	total_csd_blue_jgap  = df_jgap[1]/100* data['red']['feicom'] *(1-df_csd[1]/100)* data['blue']['feiblue']# * coupling* data['blue']['bspec']* data['blue']['fibblue']
	total_csd_red_jhgap  = df_jhgap[1]/100* data['red']['feicom'] *df_csd[1]/100* data['red']['feired'] #* coupling* data['red']['rspec']* data['red']['fibred']
	total_csd_blue_jhgap = df_jhgap[1]/100* data['red']['feicom'] *(1-df_csd[1]/100)* data['blue']['feiblue'] #* coupling* data['blue']['bspec']* data['blue']['fibblue']
	total_csd_red_hgap   = df_hgap[1]/100* data['red']['feicom'] *df_csd[1]/100* data['red']['feired']#* coupling* data['red']['rspec']* data['red']['fibred']
	total_csd_blue_hgap  = df_hgap[1]/100* data['red']['feicom'] *(1-df_csd[1]/100)* data['blue']['feiblue'] #* coupling* data['blue']['bspec']* data['blue']['fibblue']
	total_csd_red_jh     = df_jh[1]/100* data['red']['feicom'] *df_csd[1]/100* data['red']['feired']# * coupling* data['red']['rspec']* data['red']['fibred']
	total_csd_blue_jh    = df_jh[1]/100* data['red']['feicom'] *(1-df_csd[1]/100)* data['blue']['feiblue'] #* coupling* data['blue']['bspec']* data['blue']['fibblue']

	if not ploton:
		return df_jhgap[0],total_atc_jhgap, total_atc_hgap, total_atc_jgap,total_atc_jh,\
		 			 total_csd_red_jhgap, total_csd_red_hgap, total_csd_red_jgap, total_csd_red_jh,\
					 total_csd_blue_jhgap, total_csd_blue_hgap, total_csd_blue_jgap, total_csd_blue_jh
	else:
		# ATCs
		fig, ax = plt.subplots(1,figsize=(7,5))
		ax.plot(df_jhgap[0], 1- df_jhgap[1]/100,label='Dichroic only')
		ax.plot(df_block[0], df_block[1]/100,label='Cold Filter')
		#ax.plot(xao*1000,yao,label='Pre Dichroic Tracking\n Camera Throughput')
		total_atc_jhgap = (1 - df_jhgap[1]/100)* f(df_jhgap[0]) *df_block[1]/100
		ax.plot(df_jhgap[0], total_atc_jhgap ,lw=2,label='Total Throughput')
		ax.set_xlabel('Wavelength [nm]')
		ax.set_ylabel('Transmission')
		ax.legend(fontsize=10,loc=1)
		ax.set_ylim(-0.1,1.1)
		ax.set_title('JH gap')
		plt.savefig(savepath + 'total_throughput_jhgap.png')


		fig, ax = plt.subplots(1,figsize=(7,5))
		ax.plot(df_hgap[0], 1- df_hgap[1]/100,label='Dichroic only')
		ax.plot(df_block[0], df_block[1]/100,label='Cold Filter')
		#ax.plot(xao*1000,yao,label='Pre Dichroic Tracking\n Camera Throughput')
		total_atc_hgap = (1 - df_hgap[1]/100)* f(df_hgap[0])*df_block[1]/100
		ax.plot(df_hgap[0], total_atc_hgap,lw=2,label='Total Throughput')
		ax.set_xlabel('Wavelength [nm]')
		ax.set_ylabel('Transmission')
		ax.legend(fontsize=10,loc=1)
		ax.set_ylim(-0.1,1.1)
		ax.set_title('H + JH gap')
		plt.savefig(savepath + 'total_throughput_hgap.png')


		fig, ax = plt.subplots(1,figsize=(7,5))
		ax.plot(df_jgap[0], 1- df_jgap[1]/100,label='Dichroic only')
		ax.plot(df_block[0], df_block[1]/100,label='Cold Filter')
		#ax.plot(xao*1000,yao,label='Pre Dichroic Tracking\n Camera Throughput')
		total_atc_jgap = (1 - df_jgap[1]/100)* f(df_jgap[0])*df_block[1]/100
		ax.plot(df_jgap[0], total_atc_jgap,lw=2,label='Total Throughput')
		ax.set_xlabel('Wavelength [nm]')
		ax.set_ylabel('Transmission')
		ax.legend(fontsize=10,loc=1)
		ax.set_ylim(-0.1,1.1)
		ax.set_title('J + JH gap')
		plt.savefig(savepath + 'total_throughput_jgap.png')

		fig, ax = plt.subplots(1,figsize=(7,5))
		ax.plot(df_jh[0], 1- df_jh[1]/100,label='Dichroic only')
		ax.plot(df_block[0], df_block[1]/100,label='Cold Filter')
		#ax.plot(xao*1000,yao,label='Pre Dichroic Tracking\n Camera Throughput')
		total_atc_jh = (1 - df_jh[1]/100)* f(df_jh[0])*df_block[1]/100
		ax.plot(df_jh[0], total_atc_jh,lw=2,label='Total Throughput')
		ax.set_xlabel('Wavelength [nm]')
		ax.set_ylabel('Transmission')
		ax.legend(fontsize=10,loc=1)
		ax.set_ylim(-0.1,1.1)
		ax.set_title('J+H')
		plt.savefig(savepath + 'total_throughput_jh.png')

		# CSD
		fig, ax = plt.subplots(1,figsize=(7,5))
		ax.plot(df_jhgap[0], df_jhgap[1]/100,label='Dichroic Transmission')
		ax.plot(df_csd[0], df_csd[1]/100,label='CSD')
		#ax.plot(xao*1000,yao,label='Pre Dichroic Tracking\n Camera Throughput')
		ax.plot(df_jhgap[0], total_csd_red_jhgap,'r',lw=2,label='To RSPEC')
		ax.plot(df_jhgap[0], total_csd_blue_jhgap ,'b',lw=2,label='To BSPEC')
		ax.plot(df_jhgap[0], total_atc_jhgap,'k--',label='ATC transmission')
		ax.set_xlabel('Wavelength [nm]')
		ax.set_ylabel('Transmission')
		ax.legend(fontsize=10,loc=1)
		ax.set_ylim(-0.1,1.1)
		ax.set_title('Red Arm')

		ax.fill_between([980,1327],-0.1,1,color='gray',alpha=0.2)
		ax.fill_between([1490,2460],-0.1,1,color='gray',alpha=0.2)

		ax.axvline(1330,color='m',ls='--',lw=2)
		ax.axvline(1470,color='m',ls='--',lw=2)
		plt.savefig(savepath + 'total_throughput_csd.png')

		##########
		##########
		# JH gap
		fig, ax = plt.subplots(1,figsize=(10,5))
		ax.plot(df_jhgap[0], total_csd_red_jhgap,'r',lw=2,label='To RSPEC')
		ax.plot(df_jhgap[0], total_csd_blue_jhgap ,'b',lw=2,label='To BSPEC')
		ax.plot(df_jhgap[0], total_atc_jhgap,'k--',label='ATC transmission')
		ax.set_xlabel('Wavelength [nm]')
		ax.set_ylabel('Transmission')
		ax.legend(fontsize=10,loc=1)
		ax.set_ylim(-0.1,1.1)
		ax.set_title('JH gap')

		ax.fill_between([980,1327],-0.1,1,color='gray',alpha=0.2)
		ax.fill_between([1490,2460],-0.1,1,color='gray',alpha=0.2)
		ax.grid()
		ax.xaxis.grid(True, which='minor')
		ax.yaxis.grid(True, which='minor')
		ax.xaxis.set_minor_locator(MultipleLocator(5))
		ax.yaxis.set_minor_locator(MultipleLocator(0.05))

		#ax.axvline(1330,color='m',ls='--',lw=2)
		#ax.axvline(1470,color='m',ls='--',lw=2)
		plt.savefig(savepath + 'total_throughput_csd_jhgap.png')

		#########
		# J
		fig, ax = plt.subplots(1,figsize=(10,5))
		ax.plot(df_jhgap[0], total_csd_red_jgap,'r',lw=2,label='To RSPEC')
		ax.plot(df_jhgap[0], total_csd_blue_jgap ,'b',lw=2,label='To BSPEC')
		ax.plot(df_jhgap[0], total_atc_jgap,'k--',label='ATC transmission')
		ax.set_xlabel('Wavelength [nm]')
		ax.set_ylabel('Transmission')
		ax.legend(fontsize=10,loc=1)
		ax.set_ylim(-0.1,1.1)
		ax.set_title('J + JH Gap')

		ax.fill_between([980,1327],-0.1,1,color='gray',alpha=0.2)
		ax.fill_between([1490,2460],-0.1,1,color='gray',alpha=0.2)
		ax.grid()
		ax.xaxis.grid(True, which='minor')
		ax.yaxis.grid(True, which='minor')
		ax.xaxis.set_minor_locator(MultipleLocator(5))
		ax.yaxis.set_minor_locator(MultipleLocator(0.05))

		#ax.axvline(1330,color='m',ls='--',lw=2)
		#ax.axvline(1470,color='m',ls='--',lw=2)
		plt.savefig(savepath + 'total_throughput_csd_jgap.png')


		## plot ratio
		fig, ax = plt.subplots(1,figsize=(10,5))
		ax.plot(df_jhgap[0], total_csd_red_jgap/total_atc_jgap,c='r')
		ax.plot(df_jhgap[0], total_csd_blue_jgap/total_atc_jgap,c='b')
		ax.set_ylim(0,50)


		#########
		# H
		fig, ax = plt.subplots(1,figsize=(10,5))
		#ax,ax2 = axs
		ax.plot(df_jhgap[0], total_csd_red_hgap,'r',lw=2,label='To RSPEC')
		ax.plot(df_jhgap[0], total_csd_blue_hgap ,'b',lw=2,label='To BSPEC')
		ax.plot(df_jhgap[0], total_atc_hgap,'k--',label='ATC transmission')
		#ax2.plot(df_jhgap[0], total_atc_hgap * total_csd_blue ,'b--',label='ATC transmission')
		#ax2.plot(df_jhgap[0], total_atc_hgap * total_csd_red ,'r--',label='ATC transmission')
		ax.set_xlabel('Wavelength [nm]')
		ax.set_ylabel('Transmission')
		ax.legend(fontsize=10,loc=1)
		ax.set_ylim(-0.1,1.1)
		ax.set_title('H + JH Gap')

		ax.fill_between([980,1327],-0.1,1,color='gray',alpha=0.2)
		ax.fill_between([1490,2460],-0.1,1,color='gray',alpha=0.2)
		ax.grid()
		ax.xaxis.grid(True, which='minor')
		ax.yaxis.grid(True, which='minor')
		ax.xaxis.set_minor_locator(MultipleLocator(5))
		ax.yaxis.set_minor_locator(MultipleLocator(0.05))

		#ax.axvline(1330,color='m',ls='--',lw=2)
		#ax.axvline(1470,color='m',ls='--',lw=2)
		plt.savefig(savepath + 'total_throughput_csd_Hgap.png')


		#########
		# J+H
		fig, ax = plt.subplots(1,figsize=(10,5))
		ax.plot(df_jhgap[0], total_csd_red_jh,'r',lw=2,label='To RSPEC')
		ax.plot(df_jhgap[0], total_csd_blue_jh ,'b',lw=2,label='To BSPEC')
		ax.plot(df_jhgap[0], total_atc_jh,'k--',label='ATC transmission')
		ax.set_xlabel('Wavelength [nm]')
		ax.set_ylabel('Transmission')
		ax.legend(fontsize=10,loc=1)
		ax.set_ylim(-0.1,1.1)
		ax.set_title('J+H')

		ax.fill_between([980,1327],-0.1,1,color='gray',alpha=0.2)
		ax.fill_between([1490,2460],-0.1,1,color='gray',alpha=0.2)
		ax.grid()
		ax.xaxis.grid(True, which='minor')
		ax.yaxis.grid(True, which='minor')
		ax.xaxis.set_minor_locator(MultipleLocator(5))
		ax.yaxis.set_minor_locator(MultipleLocator(0.05))

		#ax.axvline(1330,color='m',ls='--',lw=2)
		#ax.axvline(1470,color='m',ls='--',lw=2)
		plt.savefig(savepath + 'total_throughput_csd_jh.png')
		
		fig, ax = plt.subplots(1,figsize=(10,5))
		ax.plot(df_jhgap[0], total_csd_red_jgap/total_atc_jgap,c='r')
		ax.plot(df_jhgap[0], total_csd_blue_jgap/total_atc_jgap,c='b')
		ax.set_ylim(0,50)


if __name__=='__main__':
	#load inputs
	#configfile = '../hispec_tracking.cfg'
	#so         = load_object(configfile)
	#cload      = fill_data(so,track_on=True)

	x,total_atc_jhgap, total_atc_hgap, total_atc_jgap,total_atc_jh,\
		 			 total_csd_red_jhgap, total_csd_red_hgap, total_csd_red_jgap, total_csd_red_jh,\
					 total_csd_blue_jhgap, total_csd_blue_hgap, total_csd_blue_jgap, total_csd_blue_jh = plot_total_throughput(ploton=False)

	
	# plot on one graph
	fig, axs = plt.subplots(3,1,figsize=(7,6),sharex=True)
	p=axs[0].plot(x, total_atc_jhgap, label='JH gap')
	axs[0].plot(x,  total_atc_hgap, label='H gap')
	axs[0].plot(x,  total_atc_jgap, label='J gap')
	axs[0].plot(x,  total_atc_jh, label='J+H')

	p=axs[1].plot(x, total_csd_red_jhgap)
	axs[1].plot(x,  total_csd_red_hgap)
	axs[1].plot(x,  total_csd_red_jgap)
	axs[1].plot(x,  total_csd_red_jh)

	p=axs[2].plot(x, total_csd_blue_jhgap)
	axs[2].plot(x,  total_csd_blue_hgap)
	axs[2].plot(x,  total_csd_blue_jgap)
	axs[2].plot(x,  total_csd_blue_jh)

	axs[0].legend()
	axs[0].set_title('ATC Forward Throughput')
	axs[1].set_title('RSPEC Total Throughput')
	axs[2].set_title('BSPEC Total Throughput')
	axs[2].set_xlabel('Wavelength [nm]')
	axs[0].set_ylabel('Throughput')
	axs[1].set_ylabel('Throughput')
	axs[2].set_ylabel('Throughput')

	#axs[0].set_xlim(1330, 1485)
	#axs[1].set_xlim(1400, 1470)
	plt.subplots_adjust(hspace=0.14)
	for i in np.arange(3):
		axs[i].set_yscale('log')
	
	# save data
	all_out =  plot_total_throughput(ploton=False)
	hdr = 'x,total_atc_jhgap, total_atc_hgap, total_atc_jgap,total_atc_jh,\
		 			 total_csd_red_jhgap, total_csd_red_hgap, total_csd_red_jgap, total_csd_red_jh,\
					 total_csd_blue_jhgap, total_csd_blue_hgap, total_csd_blue_jgap, total_csd_blue_jh\n'
	np.savetxt('all_throughputs.csv',np.array(all_out).T,delimiter=',',header=hdr)

