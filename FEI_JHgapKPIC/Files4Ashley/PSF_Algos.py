#from photutils import centroids
from scipy import ndimage, optimize
import numpy as np
#from PSF_finder_ORIG import PSF_Finder


##### Constants needed for code
INVALID_PSF = [False, -9999, -9999, -9999, -9999, -9999, -9999] 
MASK_RAD = 5.5

#-- Functions pulled from kpicctrl_Star_Tracker script
    # Modified minimally only so that they are standalone in this library as opposed
    # to relying on shms and class variables as they do in the original script.
def gen_circ_mask(im, center=None, radius = None):
    """Function to generate a circular binary mask around a point

    Args:
        im = the image (used only for determining mask size)
        center = (optional) center-point for mask (in X/Y pix coords, not row/col)
                    if None, center will be center of image
        radius = (optional) radius for mask in pix 
                    if None, radius is largest inscribed circle
    Returns:
        ndarray = boolean mask of True inside circle and False outside
    """
    # TODO::: change orientation of "center" to be row/col not X/Y
    
    # get image size
    h,w = im.shape
    
    # Deal with center and radius inputs
    if center is None:
        # use the middle of the image
        center = (int(w/2), int(h/2))
    if radius is None:
        # use largest inscribed radius
        radius = min(center[0], center[1], w-center[0], h-center[1])

    # Generate coordinate system around center
    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    # Create binary mask
    mask = dist_from_center <= radius
    return mask

def get_psf_params(im, algo_int):
    """Method to calculate and return psf parameters
    
    Args:
        im = the image to interpret
        algo_int = int representing which algo to use
    Returns:
        bool  = whether a psf was actually found
        float = the peak
        float = snr of the fit
        float = the x (width/col) coordinate of the center of the psf
        float = the y (height/row) coordinate of the center of the psf
        float = the fwhm of the psf in x
        float = the fwhm of the psf in y
    """
    # Apply selected algo
    if algo_int == 0:
        # Original KPIC Algorithm (PSF_Finder)
        func = get_psf_parameters_old
        isMasked = False
    elif algo_int == 1:
        # SciPy Center of Mass, Masked
        func = get_psf_parameters_scipyCOMmasked
        isMasked = True
    elif algo_int == 2:
        # Photutils Qaudratic Fit, Unmasked
        func = get_psf_parameters_photQDR
        isMasked = False            
    elif algo_int == 3:
        # Photutils Quadratic Fit, Masked
        func = get_psf_parameters_photQDRmasked
        isMasked = True
    elif algo_int == 4:
        # Photutils 1-D Gaussian Masked
        func = get_psf_parameters_phot1DGmasked
        isMasked = True    
    elif algo_int == 5:
        # Hacked Jason Method
        func = get_psf_parameters_JW1DG
        isMasked=False
    else:
        raise ValueError('Error in get_psf_params: algo index {} not recognized'.format(algo_int))
    
    # Deal with masking when necessary
    if isMasked:
        # Get radius for mask
        mask_rad = MASK_RAD
        # Find brightest pixel to use as center
        rind, cind = np.unravel_index(im.argmax(), im.shape)
        # generate mask
        mask = gen_circ_mask(im, center=(cind, rind), radius=mask_rad)
        # Apply mask
        im *= mask.astype(int)

    return func(im)

global last_win
last_win = np.array([0, 0, 0, 0], dtype=np.uint16)  # provide initial condition for window size
def get_psf_parameters_old(im):
    """Old method for getting the PSF
    
    Args:
        im = the image to interpret
    Returns:
        bool  = whether a psf was actually found
        float = the peak
        float = snr of the fit (to be implemented)
        float = the x (width/col) coordinate of the center of the psf
        float = the y (height/row) coordinate of the center of the psf
        float = the fwhm of the psf in x
        float = the fwhm of the psf in y
    """
    global last_win

    # crop image
    win = last_win
    try:
        im_full = im
        # window of 0, 0, 0, 0 means we want full frame
        assert list(win) != [0, 0, 0, 0]
        im_sub = im[win[2]:win[3], win[0]:win[1]]
    except: pass

    # try to find PSF (error, continue if invalid)
    try:
        ImP, _, _ = PSF_Finder(im_sub, 
                inputs = {"Radius":min([win[3]-win[2], win[1]-win[0]])}, method="GaussFit")
        used_full = False
    except:
        # in this case, we were subwindowing, so try with full frame
        # set window to full frame if it's not already
        if not (list(win) == [0, 0, 0, 0]):
            last_win[:] = 0
        if im_full is not None:
            try: 
                ImP, _, _ = PSF_Finder(im_full, method="GaussFit")
                used_full = True
            except:
                # indicate that PSF is not valid 
                return INVALID_PSF 
        else:
            # indicate that PSF is not valid in shm
            return INVALID_PSF 

    # if we got here, fit succeeded
    # x/y returned here are row/col
    cent_x = ImP['Center_Y'][0]
    # add back in subwindow
    if not used_full:
        cent_x += win[0]
    cent_y = ImP['Center_X'][0]
    # add back in subwindow
    if not used_full:
        cent_y += win[2]
    # get correct cropping (based on size of original image provided, not cred2 cropping since this is a sim...)
    lb = 0
    rb = im.shape[1] 
    ub = 0
    bb = im.shape[0]
    # for whether the psf is valid, check that we're at least 1 std away from the edge
    psf_v = cent_x >= ImP["Sigma_Y"][0]
    psf_v = psf_v and cent_x <= (rb - lb) - ImP["Sigma_Y"][0]
    psf_v = psf_v and cent_y >= ImP["Sigma_X"][0]
    psf_v = psf_v and cent_y <= (bb - ub) - ImP["Sigma_X"][0]

    # compute new subwindow
    fwhm_x = ImP['Sigma_Y'][0]*2.355
    fwhm_y = ImP['Sigma_X'][0]*2.355
    fwhm = np.max([fwhm_x, fwhm_y], 0)
    sub_win = [np.max([cent_x - 2*fwhm, 0], 0), 
           np.min([cent_x + 2*fwhm, im_full.shape[1]], 0),
           np.max([cent_y - 2*fwhm, 0], 0),
           np.min([cent_y + 2*fwhm, im_full.shape[0]], 0)]
    last_win = np.array(sub_win, dtype=last_win.dtype)

    # return
    return psf_v, ImP['Amplitude'][0], 0, cent_x, cent_y, fwhm_x, fwhm_y 

def get_psf_parameters_scipyCOMmasked(im):
    """Scipy Center of Mass Algorithm for finding PSF
    
    Args:
        im = the image to interpret (pre-masked)
    Returns:
        bool  = whether a psf was actually found
        float = the peak
        float = snr of the fit (to be implemented)
        float = the x (width/col) coordinate of the center of the psf
        float = the y (height/row) coordinate of the center of the psf
        float = the fwhm of the psf in x
        float = the fwhm of the psf in y
    """

    # try to find PSF (error, continue if invalid)
    try: 
        cent_row, cent_col = ndimage.center_of_mass(im) 
    except:
        # indicate that PSF is not valid 
        return INVALID_PSF

    # if we got here, PSF finding succeeded
    
    # Get validity of PSF finding (make sure we weren't too close to edge)
    # get correct cropping (based on size of original image provided, not cred2 cropping since this is a sim...)
    lb = 0
    rb = im.shape[1] 
    ub = 0
    bb = im.shape[0]
    # Get the mask radius
    mask_rad = MASK_RAD
    # for whether the psf is valid, check mask window fits within frame
    psf_v = cent_col >= mask_rad
    psf_v = psf_v and cent_col <= (rb - lb) - mask_rad
    psf_v = psf_v and cent_row >= mask_rad
    psf_v = psf_v and cent_row <= (bb - ub) - mask_rad

    # compute new subwindow
        # ONLY KEEPING THIS IN SO THAT THE RUNTIME IS ROUGHLY CORRECT, this algo doesn't actually use the subwindowing anyway
    #TODO::: Genuinely compute the FWHM
    fwhm_x = mask_rad
    fwhm_y = mask_rad
    fwhm = np.max([fwhm_x, fwhm_y], 0)
    sub_win = [np.max([cent_col - 3*fwhm, 0], 0), 
           np.min([cent_col + 3*fwhm, im.shape[1]], 0),
           np.max([cent_row - 3*fwhm, 0], 0),
           np.min([cent_row + 3*fwhm, im.shape[0]], 0)]

    # return
    return psf_v, im.max(), 0, cent_col, cent_row, fwhm_x, fwhm_y 

def get_psf_parameters_photQDR(im):
    """Photutils Quadratic Fitting algorithm for finding PSF
    
    Args:
        im = the image to interpret 
    Returns:
        bool  = whether a psf was actually found
        float = the peak
        float = snr of the fit (to be implemented)
        float = the x (width/col) coordinate of the center of the psf
        float = the y (height/row) coordinate of the center of the psf
        float = the fwhm of the psf in x
        float = the fwhm of the psf in y
    """

    # try to find PSF (error, continue if invalid)
    try: 
        cent_col, cent_row = centroids.centroid_quadratic(im)
    except:
        # indicate that PSF is not valid
        return INVALID_PSF

    # if we got here, PSF finding succeeded
    
    # Get validity of PSF finding (make sure we weren't too close to edge)
    # get correct cropping (based on size of original image provided, not cred2 cropping since this is a sim...)
    lb = 0
    rb = im.shape[1] 
    ub = 0
    bb = im.shape[0]
    # Get the mask radius
    mask_rad = MASK_RAD
    # for whether the psf is valid, check mask window fits within frame
    psf_v = cent_col >= mask_rad
    psf_v = psf_v and cent_col <= (rb - lb) - mask_rad
    psf_v = psf_v and cent_row >= mask_rad
    psf_v = psf_v and cent_row <= (bb - ub) - mask_rad

    # compute new subwindow
        # ONLY KEEPING THIS IN SO THAT THE RUNTIME IS ROUGHLY CORRECT, this algo doesn't actually use the subwindowing anyway
    #TODO::: Genuinely compute the FWHM
    fwhm_x = mask_rad
    fwhm_y = mask_rad
    fwhm = np.max([fwhm_x, fwhm_y], 0)
    sub_win = [np.max([cent_col - 3*fwhm, 0], 0), 
           np.min([cent_col + 3*fwhm, im.shape[1]], 0),
           np.max([cent_row - 3*fwhm, 0], 0),
           np.min([cent_row + 3*fwhm, im.shape[0]], 0)]

    # return
    return psf_v, im.max(), 0, cent_col, cent_row, fwhm_x, fwhm_y 

def get_psf_parameters_photQDRmasked(im):
    """Photutils Quadratic Fitting algorithm (with mask) for finding PSF
    
    Args:
        im = the image to interpret (pre-masked) 
    Returns:
        bool  = whether a psf was actually found
        float = the peak
        float = snr of the fit (to be implemented)
        float = the x (width/col) coordinate of the center of the psf
        float = the y (height/row) coordinate of the center of the psf
        float = the fwhm of the psf in x
        float = the fwhm of the psf in y
    """

    # try to find PSF (error, continue if invalid)
    try: 
        cent_col, cent_row = centroids.centroid_quadratic(im)
    except:
        # indicate that PSF is not valid in shm
        return INVALID_PSF

    # if we got here, PSF finding succeeded
    
    # Get validity of PSF finding (make sure we weren't too close to edge)
    # get correct cropping (based on size of original image provided, not cred2 cropping since this is a sim...)
    lb = 0
    rb = im.shape[1] 
    ub = 0
    bb = im.shape[0]
    # Get the mask radius
    mask_rad = MASK_RAD
    # for whether the psf is valid, check mask window fits within frame
    psf_v = cent_col >= mask_rad
    psf_v = psf_v and cent_col <= (rb - lb) - mask_rad
    psf_v = psf_v and cent_row >= mask_rad
    psf_v = psf_v and cent_row <= (bb - ub) - mask_rad

    # compute new subwindow
        # ONLY KEEPING THIS IN SO THAT THE RUNTIME IS ROUGHLY CORRECT, this algo doesn't actually use the subwindowing anyway
    #TODO::: Genuinely compute the FWHM
    fwhm_x = mask_rad
    fwhm_y = mask_rad
    fwhm = np.max([fwhm_x, fwhm_y], 0)
    sub_win = [np.max([cent_col - 3*fwhm, 0], 0), 
           np.min([cent_col + 3*fwhm, im.shape[1]], 0),
           np.max([cent_row - 3*fwhm, 0], 0),
           np.min([cent_row + 3*fwhm, im.shape[0]], 0)]

    # return
    return psf_v, im.max(), 0, cent_col, cent_row, fwhm_x, fwhm_y 

def get_psf_parameters_phot1DGmasked(im):
    """Photutils 1-D Gaussian Fitting algorithm (with mask) for finding PSF
    
    Args:
        im = the image to interpret (pre-masked) 
    Returns:
        bool  = whether a psf was actually found
        float = the peak
        float = snr of the fit (to be implemented)
        float = the x (width/col) coordinate of the center of the psf
        float = the y (height/row) coordinate of the center of the psf
        float = the fwhm of the psf in x
        float = the fwhm of the psf in y
    """

    # try to find PSF (error, continue if invalid)
    try: 
        cent_col, cent_row = centroids.centroid_1dg(im)
    except:
        # indicate that PSF is not valid in shm
        return INVALID_PSF

    # if we got here, PSF finding succeeded
    
    # Get validity of PSF finding (make sure we weren't too close to edge)
    # get correct cropping (based on size of original image provided, not cred2 cropping since this is a sim...)
    lb = 0
    rb = im.shape[1] 
    ub = 0
    bb = im.shape[0]
    # Get the mask radius
    mask_rad = MASK_RAD
    # for whether the psf is valid, check mask window fits within frame
    psf_v = cent_col >= mask_rad
    psf_v = psf_v and cent_col <= (rb - lb) - mask_rad
    psf_v = psf_v and cent_row >= mask_rad
    psf_v = psf_v and cent_row <= (bb - ub) - mask_rad

    # compute new subwindow
        # ONLY KEEPING THIS IN SO THAT THE RUNTIME IS ROUGHLY CORRECT, this algo doesn't actually use the subwindowing anyway
    #TODO::: Genuinely compute the FWHM
    fwhm_x = mask_rad
    fwhm_y = mask_rad
    fwhm = np.max([fwhm_x, fwhm_y], 0)
    sub_win = [np.max([cent_col - 3*fwhm, 0], 0), 
           np.min([cent_col + 3*fwhm, im.shape[1]], 0),
           np.max([cent_row - 3*fwhm, 0], 0),
           np.min([cent_row + 3*fwhm, im.shape[0]], 0)]

    # return
    return psf_v, im.max(), 0, cent_col, cent_row, fwhm_x, fwhm_y 

# Jason's Method (Hacked for simplicity)
def gauss2d(x0, y0, peak, sigma):
    """
    2d symmetric guassian function for guassfit2d
    from pyklip.fakes

    Args:
        x0,y0: center of gaussian
        peak: peak amplitude of guassian
        sigma: stddev in both x and y directions
    """
    sigma *= 1.0
    return lambda y,x: peak*np.exp( -(((x-x0)/sigma)**2+((y-y0)/sigma)**2)/2)
    
def gaussfit2d(frame, xguess, yguess, searchrad=5, guessfwhm=3, guesspeak=1):
    """
    Fits a 2d gaussian to the data at point (xguess, yguess)
    from pyklip.fakes
    
    Args:
        frame: the data - Array of size (y,x)
        xguess,yguess: location to fit the 2d guassian to (should be pretty accurate)
        searchrad: 1/2 the length of the box used for the fit
        guessfwhm: approximate fwhm to fit to
        guesspeak: approximate flux

    Returns:
        peakflux: the peakflux of the gaussian
        fwhm: fwhm of the PFS in pixels
        xfit: x position (only chagned if refinefit is True)
        yfit: y position (only chagned if refinefit is True)
    """
    if not isinstance(searchrad, int):
        raise ValueError("searchrad needs to be an integer")

    x0 = np.rint(xguess).astype(int)
    y0 = np.rint(yguess).astype(int)
    #construct our searchbox
    fitbox = np.copy(frame[y0-searchrad:y0+searchrad+1, x0-searchrad:x0+searchrad+1])

    #fit a least squares gaussian 
    #construct the residual to the fit
    errorfunction = lambda p: np.ravel(gauss2d(*p)(*np.indices(fitbox.shape)) - fitbox)

    #do a least squares fit. Note that we use searchrad for x and y centers since we're narrowed it to a box of size
    #(2searchrad+1,2searchrad+1)

    guess = (searchrad, searchrad, guesspeak, guessfwhm/(2 * np.sqrt(2*np.log(2))))

    p, success = optimize.leastsq(errorfunction, guess)

    xfit = p[0]
    yfit = p[1]
    peakflux = p[2]
    fwhm = p[3] * (2 * np.sqrt(2*np.log(2)))
   

    # convert xfit, yfit back to image coordinates
    xfit = xfit - searchrad + x0
    yfit = yfit - searchrad + y0

    return xfit, yfit, peakflux, fwhm

def get_psf_parameters_JW1DG(im):
    """Jason's 1-D Gaussian Fitting algorithm (with mask) for finding PSF
    
    Args:
        im = the image to interpret (pre-masked) 
    Returns:
        bool  = whether a psf was actually found
        float = the peak
        float = snr of the fit (to be implemented)
        float = the x (width/col) coordinate of the center of the psf
        float = the y (height/row) coordinate of the center of the psf
        float = the fwhm of the psf in x
        float = the fwhm of the psf in y
    """
    # Find peak to provide initial conditions
    rind, cind = np.unravel_index(im.argmax(), im.shape)
    peak = im[rind,cind]
    
    # try to find PSF (error, continue if invalid)
    try: 
        cent_col, cent_row, peak, fwhm = gaussfit2d(im, cind, rind, searchrad=int(np.ceil(MASK_RAD)), guessfwhm=4.35, guesspeak=peak)
    except Exception as e:
        # indicate that PSF is not valid in shm
        return INVALID_PSF

    # if we got here, PSF finding succeeded
    
    # Get validity of PSF finding (make sure we weren't too close to edge)
    # get correct cropping (based on size of original image provided, not cred2 cropping since this is a sim...)
    lb = 0
    rb = im.shape[1] 
    ub = 0
    bb = im.shape[0]
    # for whether the psf is valid, check that we're at least 1 std away from the edge
    psf_v = cent_col >= fwhm
    psf_v = psf_v and cent_col <= (rb - lb) - fwhm
    psf_v = psf_v and cent_row >= fwhm
    psf_v = psf_v and cent_row <= (bb - ub) - fwhm

    # compute new subwindow
        # ONLY KEEPING THIS IN SO THAT THE RUNTIME IS ROUGHLY CORRECT, this algo doesn't actually use the subwindowing anyway
    fwhm = np.max(fwhm, 0)
    sub_win = [np.max([cent_col - 2*fwhm, 0], 0), 
           np.min([cent_col + 2*fwhm, im.shape[1]], 0),
           np.max([cent_row - 2*fwhm, 0], 0),
           np.min([cent_row + 2*fwhm, im.shape[0]], 0)]

    # return
    return psf_v, peak, 0, cent_col, cent_row, fwhm, fwhm 

