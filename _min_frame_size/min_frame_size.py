# calc centroiding ability vs frame size
# for HISPEC ATC to guide minimum possible
# frame size in pixels
# This simulates ATC image based on user defined SNR
# and includes WFE modeled with zernikes
# environment: poppy
# ref: https://github.com/danechever/VFN-Simulations/blob/master/KPICVFN_TrackCam_Simulator/DEMO_KPIC_image_generator.m
import sys, os
import matplotlib
import numpy as np
import matplotlib.pylab as plt
from scipy import interpolate
import scipy
from datetime import date

from photutils.aperture import ApertureStats, CircularAperture
from photutils.centroids import (centroid_1dg, centroid_2dg,
                                 centroid_com, centroid_quadratic)

import poppy 

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
	given 2D image, produce a photon + read noise for it

	inputs
	-------
	psf_im - 2D image of psf

	outputs:
	--------
	electrons_out - 2d image of psf with noise added
	snr           - estimated snr frame of image
	"""
	# create random noise state
	seed          = 42
	rs            = np.random.RandomState(seed)
	
	# pull from noise the science frame with poisson noise
	science_frame = psf_im + rs.poisson(psf_im)

	# make a bias frame
	read_noise = 12 # electrons - cds noise
	bias_frame= rs.normal(scale=read_noise, size=psf_im.shape)

	# make a dark frame
	dark_noise = 0.01 # electrons - cds noise
	dark_frame = rs.poisson(scale=dark_noise, size=psf_im.shape)

	# final frame
	final_frame = science_frame + bias_frame + dark_frame
	
	snr = np.max(psf_im) /np.sqrt(np.max(psf_im) + read_noise**2)

	return electrons_out, snr

def make_ATC_image(nim=30,peaksnr=100,noise_on=False):
	"""
	inputs
	-----
	nim - size of image frame in pixels (default=30)
	nolls - 
	outputs -
	-------
	simulated image of a psf 
	"""
	# make 2D gaussian and fit center
	pixel_pitch = 18#um/pix
	#nim       = 30  #size of image in pixels
	sci_rms   = 5 / pixel_pitch # in pixels, std of rays for science spot is 5um

	x,y       = np.meshgrid(np.arange(nim), np.arange(nim)) # x y arrays corresponding to pixel position of 10x10 sub box
	amplitude = peaksnr**2  # matters if add noise
	xo,yo     = nim//2, nim//2
	sigma_x, sigma_y = 4.1/2.355,4.1/2.355 # fwhm is 4 pixels in K band for hispec
	theta, offset    = 0, 0

	# perfect PSF no convolution
	perfect_psf   = twoD_Gaussian(x,y, amplitude, xo, yo, sigma_x, sigma_y, theta, offset)

	# science PSF convolve with spot size
	sci_spot    = twoD_Gaussian(x,y, 1,nim//2, nim//2, sci_rms, sci_rms, 0, 0)
	psf         = scipy.signal.convolve2d(perfect_psf, sci_spot, mode='same', boundary='fill', fillvalue=0)
	
	if noise_on:
		psf,snr = add_noise(psf)
		snr = round(snr)
	else:
		snr='inf'

	return psf, snr

def get_peak_cen(data,ploton=False):
    """
    convolve image with  gaussian of expected fwhm to smooth out for peak finding
    
    inputs:
    ------
    data - 2d image of psf
    	must be bigger than 30x30

    ouputs:
    -------
    x3, y3 - x and y positions of psf centroid in pixels

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
    
    #xcen,ycen = x1+trim, y1+trim
    if ploton:
        plt.figure()
        plt.imshow(data)
        #plt.plot(xcen,ycen,'kx')
        plt.plot(x3,y3,'kx')
        
    return x3,y3#xcen,ycen

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

def get_track_expectation(snr,fwhm):
	"""
	get tracking expectation
	
	inputs:
	------
	snr - snr of psf
	fwhm - fwhm of psf

	outputs
	-------
	returns 1/pi * fwhm/snr (theoretical centroiding ability)
	"""
	return 1/np.pi * fwhm/snr

def getPSF_mft(Epup, lambdas, foc, coordsPP, coordsFP):
	"""
	%   Generates the point spread function (E-field itself) 
	%   - uses falco's propcustom_mft_PtoF() function 
	%   - Returns the wavelength-averaged PSF (Broadband intensity)
	%   - Also returns a cube with the E-field at each wavelength
	%
	%   Inputs:
	%       Epup: the entrance pupil field (NxNxnumWavelengths)
	%       lambdas: wavelengths (meters)
	%       foc: focal length of focusing optic in meters
	%       coordsPP: Pupil-plane coordinates structure
	%               - Must contain: N (num pix), dx (pix resolution in m)
	%       coordsFP: Focal-plane coordinates structure
	%               - Must contain: N (num pix), dx (pix resolution in m)
	%   Outputs:
	%       iPSF: 2D array with point spread function (normalized irradiance)
	%       PSF:  cube with complex-valued electric field in focal plane
	"""
    numWavelengths = len(lambdas)
    ncoords = coordsFP.N
    dx      = coordsFP.dx
    PSF     = np.zeros(ncoords,ncoords,numWavelengths);
    
    for ch in np.arange(len(numWavelengths)):
        #% Get PSF at each wavelength (channel)
        PSF[:,:,ch] = propcustom_mft_PtoF(Epup[:,:,ch],foc, lambdas[ch], coordsPP.dx, coordsFP.dx, ncoords, coordsFP.dx, ncoords);


    iPSF = mean(abs(PSF).^2,3);
    iPSF = iPSF/max(iPSF(:));

def gen_PSF():
	"""
	generate PSF from zernike coeffs - def
	"""
	nolls  = [5,8,9] + list(np.arange(50,150))
	coeffs = [-0.08, 0.05, 0.09] +  list(np.random.normal(0.0, 0.04, len(nolls)-3))
	zernike_indices = [poppy.zernike.noll_indices(nolls[i]) for i in np.arange(len(nolls))]
	zernikes   = poppy.zernike.zernike(1,1,100)
	


if __name__=='__main__':
	Tint = 0.028
	nolls = [5, 8, 9, 50:150]
	coeffs = [-0.08, 0.05, 0.09, 0.04*np.random.normal(1,len(nolls)-3)]
	nframes = 9
	# Sample tilts from 2 gaussian distributions with different parameters
	sig_x = 10;     #% [mas] rms X-jitter 
	sig_y = 10;     #% [mas] rms Y-jitter
	mu_x  = 0;      #% [mas] average X position
	mu_y  = 0;      #% [mas] average Y position
	tilts = [np.random.normal(mu_x, sig_x, nframes), np.random.normal(mu_y,sig_y,nframes)];
	# Generate images
	cred2sim   = setUp_CRED2ImSim(Npup, Nim, charge, lam0, starHmag);
	noisy, psf = CRED2ImSim_getNoisyIm(cred2sim, Tint, nolls, coeffs, tilts, True);


	psf,snr    = make_ATC_image(nim=30)
	xcen,ycen  = get_peak_cen(psf,ploton=True)



