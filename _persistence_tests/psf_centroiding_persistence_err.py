# calc the positional error for a
# gaussian PSF with noise or no noise
# as a function of persistence level
#
# see reference https://www.astro.ucla.edu/~ghezgroup/meeting_papers/fritz_limiting.pdf
# for guidance on theory
import sys, os
import matplotlib
import numpy as np
import matplotlib.pylab as plt
from scipy import interpolate
import scipy
import pandas as pd
from datetime import date

from photutils.aperture import ApertureStats, CircularAperture
from photutils.centroids import (centroid_1dg, centroid_2dg,
                                 centroid_com, centroid_quadratic)


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


def twoD_Gaussian(x,y, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    """
    """
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    return offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))

def add_noise(psf_im):
	"""
	given 2D image, produce a photon + read noise + dark frame for it
	"""
	seed       = 42
	rs         = np.random.RandomState(seed)
	shot_noise = rs.poisson(psf_im)


	read_noise = 12 # electrons - cds noise
	electrons_out = np.round(rs.normal(scale=read_noise, size=shot_noise.shape) + shot_noise)

	snr_peak  = np.max(psf_im) /np.sqrt(np.max(psf_im) + read_noise**2)
	npix = 12
	snr_total = np.sum(psf_im) /np.sqrt(np.sum(psf_im) + npix * read_noise**2)

	return electrons_out, snr_peak, snr_total

def get_peak_cen(data,ploton=False):
    """
    convolve image with  gaussian of expected fwhm to smooth out for peak finding
    """
    xtarr,ytarr = np.meshgrid(np.arange(30), np.arange(30)) # x y arrays corresponding to pixel position of 10x10 sub box
    gaus = twoD_Gaussian(xtarr,ytarr, 1,15, 15, 4/2.355, 4/2.355, 0, 0)
    data_conv  = scipy.signal.convolve2d(data, gaus, mode='same', boundary='fill', fillvalue=0)

    trim       = 2
    sub_data   = data_conv[trim:-1*trim, trim:-1*trim]

    # fit for centroid
    x1,y1 = centroid_2dg(data)
    x2,y2 = centroid_1dg(data)
    x3,y3 = centroid_quadratic(data)
    
    #xcen,ycen = x3+trim, y3+trim
    if ploton:
        plt.figure()
        plt.imshow(data)
        plt.plot(xcen,ycen,'kx')
        #plt.plot(x1,y1,'kx')
        
    return x3,y3#xcen,ycen

def calc_persistence_bias(snr,contrast,offset,noise_on=False):
	"""
	make guassian 2d psf and convolve with optical radius for science beam
	and for ghost beam

	Add the two images, add noise if turned on, then compute peak center error
	"""
	# make 2D gaussian and fit center
	pixel_pitch = 18#um/pix
	nim       = 30
	sci_rms   = 5 / pixel_pitch # in pixels

	x,y       = np.meshgrid(np.arange(nim), np.arange(nim)) # x y arrays corresponding to pixel position of 10x10 sub box
	amplitude = snr**2  # matters if add noise
	xo,yo     = nim//2, nim//2
	sigma_x, sigma_y = 4.1/2.355,4.1/2.355 # fwhm is 4 pixels in K band for hispec
	theta     = 0

	# perfect PSF no convolution, persistent spot amplitude is lower by contrast
	# sci psf at 0, 
	perfect_psf   = twoD_Gaussian(x,y, amplitude, xo, yo, sigma_x, sigma_y, theta, 0)
	perfect_psf_g = twoD_Gaussian(x,y, contrast*amplitude, xo+offset, yo, sigma_x, sigma_y, theta, 0)

	#  convolve with spot size, persistent spot has same rms and image qualities
	sci_spot    = twoD_Gaussian(x,y, 1,nim//2, nim//2, sci_rms, sci_rms, 0, 0)
	psf         = scipy.signal.convolve2d(perfect_psf, sci_spot, mode='same', boundary='fill', fillvalue=0)
	psf_g       = scipy.signal.convolve2d(perfect_psf_g, sci_spot, mode='same', boundary='fill', fillvalue=0)

	if noise_on:
		psf,snr_peak, snr_total = add_noise(psf)
		psf_g,_,_ = add_noise(psf_g)
		snr_peak  = round(snr_peak)
		snr_total = round(snr_total)
	else:
		snr_peak =np.inf
		snr_total=np.inf		

	xcen,ycen      = get_peak_cen(psf,ploton=False)
	xcen_g, ycen_g = get_peak_cen(psf+psf_g,ploton=False)

	return xcen - xcen_g, psf,psf_g, snr_peak, snr_total


def calc_peak_to_total_flux_ratio(amplitude):
	"""
	If peak flux is X what is total flux for a 4.1 px FWHM gaussian?

	answer is 0.05 (amp is 1/20th of total flux)
	"""
	# make 2D gaussian and fit center
	pixel_pitch = 18 #um/pix
	nim       = 30   # size of image

	x,y       = np.meshgrid(np.arange(nim), np.arange(nim)) # x y arrays corresponding to pixel position of 10x10 sub box
	#amplitude = peaksnr**2  # matters if add noise
	xo,yo     = nim//2, nim//2
	fwhm = 4.1 # pixels  - size of PSF
	sigma_x, sigma_y = fwhm/2.355,fwhm/2.355 # fwhm is 4 pixels in K band for hispec
	theta, offset    = 0, 0

	# perfect PSF no convolution
	perfect_psf   = twoD_Gaussian(x,y, amplitude, xo, yo, sigma_x, sigma_y, theta, offset)

	# total flux
	total_flux = np.sum(perfect_psf)

	print(amplitude/total_flux)


def get_track_req():
	# compute requirement. requirement is 0.2lambda/D in y band
	fratio = 35
	D      = 10
	pixel_pitch = 18
	platescale_arcsec_um = 206265 / fratio / (D * 10**6) #arc/um
	platescale_arcsec_pix = platescale_arcsec_um * pixel_pitch

	yband_wavelength       = 1020 # nm, center of y band
	tracking_requirement_arcsec = 206265 * 0.2 * yband_wavelength / (D*10**9) 
	tracking_requirement_pixel  = tracking_requirement_arcsec/platescale_arcsec_pix

	return tracking_requirement_pixel

def get_track_exp(snr,fwhm):
	"""
	get tracking expectation
	"""
	return 1/np.pi * fwhm/snr

if __name__=='__main__':
	#load inputs
	offsets   = [0.1,0.5, 1,1.5, 2, 3,5] # in case the psf is offset
	contrasts = [1e-5,1e-4,1e-3,1e-2,1e-1] 
	snr       = 20
	noise_on  = True

	# Compute offsets and plot!
	pix_offsets = np.zeros((len(offsets),len(contrasts)))
	
	plt.figure()
	for i, offset in enumerate(offsets):
		for j, contrast in enumerate(contrasts):
			pix_offsets[i,j], psf,psf_g, snr_peak,snr_total = calc_persistence_bias(snr,\
												contrast=contrast,\
												offset=offset,\
												noise_on=noise_on)

		plt.loglog(contrasts, np.abs(pix_offsets[i,:]),'-o',label='Offset: %spx'%offset)

	plt.xlabel('Persistent Flux to Science Ratio')
	plt.ylabel('Centroiding Error [pix]')
	track_req = get_track_req()
	try:
		track_expect = get_track_exp(snr_peak,fwhm=4.1)
		track_expect_total = get_track_exp(snr_total,fwhm=4.1)
		#plt.axhline(track_expect,c='g',ls='--',label='expectation (peak)')	
		plt.axhline(track_expect_total,c='g',ls='-.',label='expectation (tot)')	
	except TypeError:
		pass
	plt.axhline(track_req,c='k',ls='--',label='requirement')

	plt.legend(fontsize=10)
	plt.subplots_adjust(bottom=0.15,left=0.2)
	plt.grid()
	plt.ylim(1e-7,1)
	plt.title('Impact on Tracking for Persistent\n PSFs Ranging in Offset and Contrast',fontsize=14)
	plt.savefig('ghost_PSFcentroid_error_noise_%s_peaksnr_%s.png'%(noise_on,snr_peak))

	#
	# plot image of psf and ghost
	offset, psf,psf_g, snr_peak,snr_total = calc_persistence_bias(snr,\
										contrast=contrast,\
										offset=offset,\
										noise_on=noise_on)

	plt.figure()
	# plt.imshow(psf)
	plt.imshow(psf+psf_g)
	if not noise_on:
		plt.contour(psf,colors='r',alpha=0.5)
		plt.contour(psf_g,colors='w',alpha=0.2)
	if noise_on: plt.title('SNR=%s'%snr)
	if not noise_on: plt.title('offset=%spx'%round(np.abs(offset),2))
	plt.savefig('example_psf_snr_%s_noise_%s.png'%(snr,noise_on))




