"""Show phase and amplitude at key planes in the Habex PROPER model."""
# Copyright 2020, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.
# -------------------------------------------------------------------------
# Script to test the Habex PROPER prescription in Python
#
# Written by A.J. Riggs (JPL, CIT) in February, 2020
#
# PYTHONPATH (or sys.path) must know the locations of
# - the Habex prescription 'habex'
# - the PROPER library
# - the FALCO package

import numpy as np
from astropy.io import fits  # for reading in FITS files
import matplotlib.pyplot as plt
import proper
import sys
sys.path.insert(0, '/Users/ajriggs/Repos/proper-models/simple_habex')

map_dir = '/Users/ajriggs/Documents/habex/maps/'  # Change for your computer

prescription = 'habex'

# Focal Plane values
res = 3  # pixels per lambda0/D
FOV = 30  # FOV radius
Rout = 26  # field stop radius in lambda0/D

# Pupil Plane Resolutions
NbeamFull = 62*7  # points across the beam at the pupil
gridsize = int(2**np.ceil(np.log2(NbeamFull)))  # points across the pupil array

lambda0 = 550e-9  # Central wavelength of the whole spectral bandpass [meters]
lambda_um = 550e-9 * 1e6

optval = dict([
        ('map_dir', map_dir),
        ('nout', int(2*np.ceil(0.5*(1 + res*(2*FOV))))),
        ('final_sampling_lam0', 1/res),
        ('use_field_stop', True),
        ('field_stop_radius', Rout),
        ('lambda0_um', lambda0*1e6),
        ('pupil_diam_pix', NbeamFull),
        ('normLyotDiam', 0.95),
        ('vortexCharge', 6),
        ('xoffset', 0),
        ('yoffset', 0)
            ])

# %% No Correction

# Phase Retrieval of Pupil
optval['use_pr'] = True
optval['pr_pupil_diam_pix'] = 248
optval['use_errors'] = 1  # 1 = use optical surface errors, 0 = none
optval['use_fpm'] = 0  # use focal plane mask (0 = no FPM)
optval['use_lyot_stop'] = 1  # use Lyot stop (0 = no stop)
optval['use_field_stop'] = 0  # use field stop (0 = no stop)
[Epup, sampling_m] = proper.prop_run(prescription, lambda_um, gridsize,
                                     PASSVALUE=optval, QUIET=True)

mask = np.zeros(Epup.shape)
mask[np.abs(Epup) > 0.1*np.max(np.abs(Epup))] = 1

plt.figure(1); plt.imshow(np.abs(Epup)); plt.colorbar(); plt.title('abs(Epupil)'); plt.pause(1e-2)
plt.figure(2); plt.imshow(mask*np.angle(Epup)); plt.colorbar(); plt.clim(-1, 1); plt.title('angle(Epupil)'); plt.pause(1e-2)

# PSF for normalization
optval['use_pr'] = False
optval['use_errors'] = 1;		#-- 1 = use optical surface errors, 0 = none
optval['use_fpm'] = 0;		#-- use focal plane mask (0 = no FPM)
optval['use_lyot_stop'] = 1 	#-- use Lyot stop (0 = no stop)
optval['use_field_stop'] = 0 	#-- use field stop (0 = no stop)
[EforNorm, sampling_m] = proper.prop_run(prescription, lambda_um, gridsize, PASSVALUE=optval, QUIET=True )
IforNorm = np.abs(EforNorm)**2
I00 = np.max(IforNorm)
plt.figure(3); plt.imshow(np.log10(IforNorm/I00)); plt.colorbar(); plt.clim(-7, 0); plt.title('PSF for Normalization)'); plt.pause(1e-2)

# Coronagraphic PSF
optval['use_pr'] = False
optval['use_errors'] = 1;		#-- 1 = use optical surface errors, 0 = none
optval['use_fpm'] = 1;		#-- use focal plane mask (0 = no FPM)
optval['use_lyot_stop'] = 1;	#-- use Lyot stop (0 = no stop)
optval['use_field_stop'] = 1;#1;	#-- use field stop (0 = no stop)
[Ecoro, sampling_m] = proper.prop_run(prescription, lambda_um, gridsize, PASSVALUE=optval, QUIET=True )
Inorm = (np.abs(Ecoro)**2) / I00
plt.figure(4); plt.imshow(np.log10(Inorm)); plt.title('Aberrated Coronagraphic PSF'); plt.colorbar(); plt.clim(-7, -2); plt.pause(1)


# %% With Phase Flattened at Pupil before the FPM

optval['use_dm1'] = True
optval['dm1'] = fits.getdata(optval['map_dir']+'flat_map.fits')

#--Phase Retrieval of Pupil
optval['use_pr'] = True
optval['pr_pupil_diam_pix'] = 248;
optval['use_errors'] = 1;		#-- 1 = use optical surface errors, 0 = none
optval['use_fpm'] = 0;		#-- use focal plane mask (0 = no FPM)
optval['use_lyot_stop'] = 1;	#-- use Lyot stop (0 = no stop)
optval['use_field_stop'] = 0;	#-- use field stop (0 = no stop)
[Epup, sampling_m] = proper.prop_run(prescription, lambda_um, gridsize, PASSVALUE=optval, QUIET=True )

mask = np.zeros(Epup.shape)
mask[np.abs(Epup) > 0.1*np.max(np.abs(Epup))] = 1

plt.figure(11); plt.imshow(np.abs(Epup)); plt.colorbar(); plt.title('abs(Epupil)'); plt.pause(1e-2)
plt.figure(12); plt.imshow(mask*np.angle(Epup)); plt.colorbar(); plt.clim(-1, 1); plt.title('angle(Epupil)'); plt.pause(1e-2)


# PSF for normalization
optval['use_pr'] = False
optval['use_errors'] = 1;		#-- 1 = use optical surface errors, 0 = none
optval['use_fpm'] = 0;		#-- use focal plane mask (0 = no FPM)
optval['use_lyot_stop'] = 1 	#-- use Lyot stop (0 = no stop)
optval['use_field_stop'] = 0 	#-- use field stop (0 = no stop)
[EforNorm, sampling_m] = proper.prop_run(prescription, lambda_um, gridsize, PASSVALUE=optval, QUIET=True )
IforNorm = np.abs(EforNorm)**2
I00 = np.max(IforNorm)
plt.figure(13); plt.imshow(np.log10(IforNorm/I00)); plt.colorbar(); plt.clim(-7, 0); plt.title('PSF for Normalization)'); plt.pause(1e-2)


# Coronagraphic PSF
optval['use_pr'] = False
optval['use_errors'] = 1;		#-- 1 = use optical surface errors, 0 = none
optval['use_fpm'] = 1;		#-- use focal plane mask (0 = no FPM)
optval['use_lyot_stop'] = 1;	#-- use Lyot stop (0 = no stop)
optval['use_field_stop'] = 1;	#-- use field stop (0 = no stop)
[Ecoro, sampling_m] = proper.prop_run(prescription, lambda_um, gridsize, PASSVALUE=optval, QUIET=True )
Inorm = np.abs(Ecoro)**2 / I00
plt.figure(14); plt.imshow(np.log10(Inorm)); plt.title('Coronagraphic PSF after Flattening'); plt.clim(-7, -2); plt.colorbar(); plt.pause(5)
