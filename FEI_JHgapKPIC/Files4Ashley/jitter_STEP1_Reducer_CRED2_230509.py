##########################################
# Script to process CRED2 image cubes and determine the PSF position
# in each each frame of the cube. This is to reduce the data once so that
# we can then plot at will without having to re-run reduction every time.
#
# The PSF finding is done with the JW1DG algo which has proven most accurate
# in multiple simulations
#
# The output of this script can be plotted with the CRED2_Jitter_Plotter.py
# script.
#
#
# NOTE: this copy is based off the version in the AlgoDev/TMP_2ShareWGary/230507 directory
#   - This version uses the new "COPYTAGS" header keyword to deal with tags
# NOTE: The MaskRad in PSF_Algos has been changed to 4 for this analysis since that seems
#       to help with the VFN-mode PSF finding.
##########################################
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import time
from glob import glob
import os

#-- KPIC Imports
import PSF_Algos

startTime = time.time()

#-- Filenames and stuff
# Cred2 Stuff (modify paths as needed)
fl_creddir = 'rawdata/CRED2/'
fl_bkgddir = 'rawdata/CRED2/'

specify_bkgd = True   # Flag to denote if specific background should be used. If so, provide bkgd file below
#(OPTIONAL) Pull a more recent BG file
fl_bkgd = fl_bkgddir+'cred2_0068.fits'

#-- Option to save results or not
isSave = True
isSavePlot = True   # Option to save plot rather than show it during reduction
savePath = fl_creddir+'/reduced/'
if isSave and not os.path.exists(savePath):
    # Create save directory if it doesn't exist
    os.makedirs(savePath)


#-- String identifier to use with glob to identify files to reduce
globelement = 'CL_Test*.fits'

#-- Analysis constants
fluxwin = 5.5       # window size for photometric calc
plate_scale = 7.240      # mas/pix
lamOD = 2.2/(10e6) * 206265*1000            # Lambda over D of Keck [mas]

####### Main code ######
#-- Define a function to load data
#def loader(fl_num, pth=fl_creddir, pre=fl_credpre):
#    flnm = pth+pre+f'{fl_num:04d}'
#    with fits.open(flnm+'.fits') as f:
#        data = f[0].data
#        head = f[0].header
#    return (data, head, flnm)


if not isSavePlot:
    plt.ion()

#-- Find files to analyze
flnms = glob(fl_creddir+globelement)

#-- Exclude files that were already analyzed
for flnm in flnms:
    print(flnm)

#-- Iterate through data doing analysis
for ind, flnm in enumerate(flnms):
    fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(14, 4.4))
    flnm_short = os.path.basename(flnm)
    print('Frame %d of %d'%(ind+1, len(flnms)))
    print('\t'+flnm_short)
    # Load the data frame
    with fits.open(flnm, memmap=True) as f:
        data = f[0].data
        hdr = f[0].header
        tims = f[1].data

    if not specify_bkgd:
        # Find the corresponding background
        fl_bkgd = os.path.basename(hdr['BIAS'])
        fl_bkgd = fl_bkgddir + fl_bkgd
    # Pull the background
    with fits.open(fl_bkgd, memmap=False) as f:
        bkgd = f[0].data

    if specify_bkgd and (bkgd.ndim == 3):
        # A non-standard background was loaded, condense temporal dimension 
        bkgd = np.median(bkgd, axis=1)
    
    if bkgd.shape != data.shape:
        print('Background shape does not match data shape; cropping background')
        cropcen = [hdr['CROPWIN0'], hdr['CROPWIN1']]
        # Crop background to match cropping of raw data
        rawdatsize = data.shape[-2:]
        bkgd = bkgd[cropcen[0]-rawdatsize[0]//2:cropcen[0]+rawdatsize[0]//2+1, cropcen[1]-rawdatsize[1]//2:cropcen[1]+rawdatsize[1]//2+1]
    # Subtract background
    data = data - bkgd
    if hdr['COPYTAGS']:
        # Deal with bad pixels via matching to adjacent values
        data[:,0,:4] = np.expand_dims(data[:,0,4], axis=1)

    # Perform fit on the first frame to find the rough PSF position
    im = data[0]
    _, _, _, cent_x, cent_y, fwhm_x, fwhm_y = PSF_Algos.get_psf_params(im.copy(), 5)
    # Rename coordinate system since x,y always messes me up... who works in XY when you can work in Row,Col?!?!
    colcen, rowcen = cent_x, cent_y
    fwhm_c, fwhm_r = fwhm_x, fwhm_y
    # Display results
    # Full image for reference
    ax1.imshow(im)
    ax1.set_title("Original")
    # Image cropped around found PSF
    ax2.imshow(im)
    # Add crosshair on found-PSF
    ax2.plot(colcen, rowcen, 'r+')  
    ax2.set_title("Fit - Coords (%0.2f, %0.2f) - Size (%0.2f, %0.2f)"%(colcen,rowcen,fwhm_c,fwhm_r))
    # Actually crop the image display now
    cropwin = 8
    ax2.set_xlim([colcen-cropwin, colcen+cropwin+1]); ax2.set_ylim([rowcen-cropwin,rowcen+cropwin+1])
   
    # Preallocate matrix for final results
        # Rows = frames in data, 
        # Col1 = rowCenter, Col2 = colCenter, Col3 = rowSigma, Col4 = colSigma, Col5 = encircled Flux
    Coords = np.full((data.shape[0], 5), np.nan) 
    # Get mask size fora full, well-centered mask
    masksize = PSF_Algos.gen_circ_mask(data[0].copy(), center=np.array(data[0].shape[::-1])//2, radius=fluxwin).sum()
    # Iterate through frames finding PSF and save coords
    for find in range(data.shape[0]):
        if (find%(data.shape[0]//100)) == 0:
            print('    Analyzing %d of %d'%(find+1, data.shape[0]))
        # Perform fit
        try:
            _, _, _, cent_x, cent_y, fwhm_x, fwhm_y = PSF_Algos.get_psf_params(data[find,:,:].copy(), 5)
            # swap coords again cuz no like x,y 
            colcen, rowcen = cent_x, cent_y     
            fwhm_c, fwhm_r = fwhm_x, fwhm_y
            # get encircled flux 
                # Use integer center coords so that mask size is constant between frames
            mask = PSF_Algos.gen_circ_mask(data[find].copy(), center=(int(colcen), int(rowcen)), radius=fluxwin)
            flx = (mask*data[find]).sum() 
            # if PSF at edge, don't compute flux since sum skewed by missing pixels w.r.t frames where PSF fully in frame
                # (old edge check was to size of array vs. size of mask radius)
            # if np.array( [np.array([rowcen,colcen]) < fluxwin , (data[find].shape-np.array([rowcen,colcen])) < fluxwin] ).any():
                # (new edge check just makes sure mask is full size for all frames)
            if mask.sum() != masksize:
                flx = -9999.
            # (old square calculation method below for reference)
            # flx = data[find,int(rowcen)-fluxwin:int(rowcen)+fluxwin, int(colcen)-fluxwin:int(colcen)+fluxwin].sum()
        except KeyboardInterrupt as e:
            print('keyboard interrupt detected')
            raise e 
        except Exception as e:
            print('--- Fit failed on %d:%d'%(ind,find))
            # Set the values to be obviosly wrong (using default "wrong" value from PSF fitting lib)
            res = PSF_Algos.INVALID_PSF
            # swap coords again cuz no like x,y 
            colcen, rowcen = res[3], res[4]   
            fwhm_c, fwhm_r = res[5], res[6]
            flx = -9999.
        # Save coords 
        Coords[find, :] = np.array([rowcen, colcen, fwhm_r, fwhm_c, flx]).T

    # Make a copy of the coords result to be able to save "raw" coords later
    Coords_orig = Coords.copy()

    # Mask out bad coords entries (from where fit failed)
    badentries = np.where(Coords[:,0] == PSF_Algos.INVALID_PSF[4])
    Coords[badentries,:2] = np.nan
    

    # Display results summary
    #print(Coords[:,0], Coords[:,1])
    std_vals = (np.nanstd(Coords[:,0]*plate_scale), np.nanstd(Coords[:,1]*plate_scale))
    print('-- STDs: row=%f [mas]    col=%f [mas]'%std_vals)

    #-- Plot 2D displacement
    markerSize = 8
    ax3.scatter(Coords[:,1]-np.nanmean(Coords[:,1]), Coords[:,0]-np.nanmean(Coords[:,0]),markerSize,c=np.linspace(0,1,len(Coords[:,1])))
    ax3.set_xlabel('Displacement [pix]')
    ax3.set_ylabel('Displacement [pix]')
    ax3.set_title('PSF Positions (relative to mean)')
    #- Axis limits (fixed axis size)
    axlim = 5
    ax3.set_xlim([-axlim, axlim]); ax3.set_ylim([-axlim,axlim])
    ax3.set_aspect('equal')
    
    #-- Add title to main figure panel marking which frame this is and some analysis properties
    tit = flnm_short
    tit += '\nSTDC = {:0.2f}, STDR = {:0.2f} [pix]'.format(std_vals[0],std_vals[1])
    tit += '\nPSF Finding failed on %d frames'%badentries[0].shape
    plt.suptitle(tit)   

    plt.tight_layout()
    
    executionTime = (time.time() - startTime)
    print('Total Runtime so far: %f'%executionTime)
    
    if isSave: 
        #-- Save results into new fits file
        # Add analysis properties to the fits header
        hdr.append(('BKG4PROC', os.path.basename(fl_bkgd), 'File used by STEP1_Reducer in analysis'))
        fitsfl = fits.HDUList([fits.PrimaryHDU(data=Coords_orig, header=hdr), fits.ImageHDU(data=tims)])
        fitsfl.writeto(os.path.splitext(flnm)[0]+'_proc.fits')

    if isSavePlot:
        #-- Save plot 
        plt.savefig(savePath+os.path.splitext(flnm_short)[0]+'step1plot.pdf')
        plt.close()
        
