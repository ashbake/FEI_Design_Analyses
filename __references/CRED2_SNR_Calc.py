# dan's code for kpic cred2 calculation for comparison
import numpy as np

H_zpt   = 1.15e-9   # [W/m2/um] zero point (Allen's Astrophysical Quantities Table 7.5)
Area    = 76        # [m2]      Keck telescope area
lam0    = 1.578     # [um]      central wavelength for tracking
dlam    = 0.15      # [um]      wavelength coverage on CRED2 (assume 1.475 - 1.70 um)       
T       = 0.3       # [e-/ph]   Throughput to CRED2 (including QE) -- assuming PyPO out but vortex IN
FWHM    = 4.15      # [px]      FWHM of PSF
H_mag   = 11        # []        H-band magnitude to assume for star
DC      = 315       # [e-/s/px] dark current (CRED2 spec)
RN      = 30        # [e-/px]   read noise (CRED2 spec)
FPS     = 27        # [Hz]      frame rate to assume for CRED2 (sets itime as 1/FPS)
e_max   = 33e3      # [e-]      well-depth for CRED2 (CRED2 spec for high gain mode

# Photon energy at central wavelength
h = 6.26068e-34     # [m2.kg/s] Planck's constant
c = 299792458       # [m/s]     speec of light
ph_energy = h * c / (lam0 * 1e-6)       # [Joules/ph: W.s/ph] 

# Compute stellar flux
Phi     = H_zpt * 10**(-H_mag/2.5)  # [W/m2/um]     flux per um in H-band
Phi     /= ph_energy                # [ph/m2/s/um]  photon flux per um in H band
Phi     *= dlam                     # [ph/m2/s]     photon flux in bandpass
Phi     *= Area                     # [ph/s]        photon flux in keck aperture

# Compute misc. terms
tau     = 1/FPS             # [s]   integration time on CRED2
psf_A   = np.pi*(FWHM/2)**2           # [px]  number of pixels covered by the PSF

# Compute signals
S_star  = 0.8*Phi * T * tau     # [e-]   electrons from star
S_DC    = DC * tau * psf_A  # [e-]   Dark Current

# Compute SNR
SNR     = S_star / np.sqrt(np.sqrt(S_star)**2 + np.sqrt(S_DC)**2 + psf_A * RN**2)   # [sqrt(e-)]    SNR on the PSF

# Compute total signal (on each pixel)
S_tot = (S_star + S_DC) / psf_A

print('Star Mag:    %0.2f'%H_mag)
print('FPS: %d  | itime: %0.2f [ms]'%(FPS, tau*1e3))
print('e-/pix: %d'%S_tot)
print('===>>> SNR: %0.1f'%SNR)
print((1/np.pi) * FWHM/SNR)
if S_tot > e_max:
    print('*** SATURATED ({:,d} e- signal > {:,d} e- well-depth)'.format(int(S_tot), int(e_max)))
