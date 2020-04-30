"""Run WFC with FALCO and the Habex PROPER model."""
# Copyright 2020, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.
# -------------------------------------------------------------------------
#
# Script to perform WFSC for the Habex vortex model.
#
# Written by A.J. Riggs (JPL, CIT) in February 2020.
#
# PYTHONPATH (or sys.path) must know the locations of
# - the Habex prescription 'habex.py'
# - the PROPER library
# - the FALCO package
#
# In falco_defaults_Habex_VC.py, change the value of mp.full.map_dir to be
# for your computer.
# -------------------------------------------------------------------------

# Change this directory to your own, or add your directory to PYTHONPATH
import sys
sys.path.insert(0, '/Users/ajriggs/Repos/proper-models/simple_habex')

import numpy as np
import copy
from scipy.interpolate import RectBivariateSpline
import matplotlib.pyplot as plt
from astropy.io import fits

import falco
import proper
import falco_defaults_Habex_VC as DEFAULTS
mp = DEFAULTS.mp


mp.path = falco.config.Object()

mp.path.falco = './'  # Location of FALCO. Change to be correct for your machine

# Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
# mp.path.config = './' #--Location of config files and minimal output files. Default is [mainPath filesep 'data' filesep 'brief' filesep]
# mp.path.ws = './' # (Mostly) complete workspace from end of trial. Default is [mainPath filesep 'data' filesep 'ws' filesep];


# Step 3: Overwrite default values as desired

# ##--Special Computational Settings
mp.flagPlot = True
mp.flagMultiproc = True  # whether to use multiprocessing to parallelize some large computations
# mp.Nthreads = 2  # Number of threads to use when using multiprocessing. If undefined, it is set to the number of cores 

# Record Keeping
mp.TrialNum = 1
mp.SeriesNum = 1

# DEBUGGING and DEMONSTRATION: Use minimal settings
mp.fracBW = 0.01  # fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 1  # Number of sub-bandpasses to divide the whole bandpass into for estimation and control
mp.Nwpsbp = 1  # Number of wavelengths to used to approximate an image in each sub-bandpass
mp.Nitr = 3  # Number of wavefront control iterations
# mp.flagMultiproc = False

# %% Step 3b: Obtain the phase retrieval phase.

optval = copy.copy(vars(mp.full))
optval['xoffset'] = 0
optval['yoffset'] = 0
optval['use_dm1'] = True
optval['dm1'] = fits.getdata(mp.full.map_dir+'flat_map.fits')

optval['use_pr'] = True
# optval['end_at_fpm_exit_pupil'] = True
# optval['output_field_rootname'] = [fileparts(mp.full.input_field_rootname) filesep 'fld_at_xtPup'];
optval['use_fpm'] = False
nout = mp.P1.full.Narr			# nout > pupil_diam_pix
optval['output_dim'] = nout  # Get the Input Pupil's E-field

if mp.Nsbp == 1:
    lambdaFacs = (1,)
else:
    lambdaFacs = np.linspace(1-mp.fracBW/2, 1+mp.fracBW/2, mp.Nsbp)

prescription = 'habex';

# Get the Input Pupil's E-field
Nf = nout #--N full
Nc = falco.util.ceil_even((mp.P1.compact.Nbeam/mp.P1.full.Nbeam)*Nf)
mp.P1.compact.E = np.ones((Nc, Nc, mp.Nsbp), dtype=complex)
# mp.P1.compact.E = ones(mp.P1.compact.Nbeam+2, mp.P1.compact.Nbeam+2, mp.Nsbp); %--Initialize
for si in range(mp.Nsbp):
    lambda_um = 1e6*mp.lambda0*lambdaFacs[si]
    [fldFull, sampling_m] = proper.prop_run(prescription, lambda_um, nout, PASSVALUE=optval, QUIET=True)
    if(mp.flagPlot):
        pass
        plt.figure(605); plt.imshow(np.angle(fldFull)); plt.colorbar; plt.hsv(); plt.pause(1e-2)
        plt.figure(606); plt.imshow(np.abs(fldFull)); plt.colorbar; plt.magma(); plt.pause(1e-2)

    # Downsampling for the compact model
    dxF = 1
    dxC = mp.P1.full.Nbeam/mp.P1.compact.Nbeam

    xF = np.arange(-Nf/2., Nf/2.)*dxF  # (-Nf/2:(Nf/2-1))*dxF;
    xC = np.arange(-Nc/2., Nc/2.)*dxC  # (-Nc/2:(Nc/2-1))*dxC;
    interp_spline_real = RectBivariateSpline(xF, xF, np.real(fldFull))
    interp_spline_imag = RectBivariateSpline(xF, xF, np.imag(fldFull)) 
    fldC = interp_spline_real(xC, xC) + 1j*interp_spline_imag(xC, xC)

    # Ncrop = falco.util.ceil_even(mp.P1.compact.Nbeam+1)
    # fldC = falco.util.pad_crop(fldC, (Ncrop, Ncrop))
    if(mp.flagPlot):
        plt.figure(607); plt.imshow(np.angle(fldC)); plt.colorbar; plt.hsv(); plt.pause(1e-2)
        plt.figure(608); plt.imshow(np.abs(fldC)); plt.colorbar; plt.magma(); plt.pause(1e-2)

    # Assign to initial E-field in compact model.
    mp.P1.compact.E[:, :, si] = falco.prop.relay(fldC, 1)


# %% Obtain DM1 Commands to Flatten the Wavefront Prior to WFSC

# % optval['use_pr = true;
# % optval['pr_pupil_diam_pix = mp.P1.compact.Nbeam;
# %
# % optval['use_errors = 1;		%-- 1 = use optical surface errors, 0 = none
# % % optval['use_fpm = 0;		%-- use focal plane mask (0 = no FPM)
# % % optval['use_lyot_stop = 0;	%-- use Lyot stop (0 = no stop)
# % % optval['use_field_stop = 0;	%-- use field stop (0 = no stop)
# % [Epup, sampling_m]  = habex_vortex(mp.lambda0, mp.P1.full.Narr, optval);
# %
# % mask = 0*Epup;
# % mask(abs(Epup) > 1e-1*max(abs(Epup(:)))) = 1;
# % % mask = ones(size(Epup));
# %
# % surfaceToFit = -0.5*mask.*angle(Epup)*(mp.lambda0/(2*pi));
# % 
# % figure(1); imagesc(abs(Epup)); axis xy equal tight; colorbar; 
# % % title('', 'Fontsize', 16);
# % drawnow;
# % figure(2); imagesc(surfaceToFit); axis xy equal tight; colorbar; 
# %
# %
# % mp.dm1.inf0 = fitsread(mp.dm1.inf_fn);
# % mp.dm1.dm_spacing = 400e-6;
# % mp.dm1.dx_inf0 = mp.dm1.dm_spacing/10;
# % mp.dm1.dx = mp.P2.D/mp.P2.compact.Nbeam;
# % mp.dm1.centering = 'pixel';
# % V0 = falco_fit_dm_surf(mp.dm1,surfaceToFit);
# % fitswrite(V0, [mp.full.map_dir, 'flat_map.fits'])
# % figure(3); imagesc(V0); axis xy equal tight; colorbar;


# Step 4: Generate the label associated with this trial
mp.runLabel = 'Series' + ('%04d'%(mp.SeriesNum)) + '_Trial' + ('%04d_'%(mp.TrialNum)) + mp.coro + '_' \
+ mp.whichPupil + '_' + str(np.size(mp.dm_ind)) + 'DM' + str(mp.dm1.Nact) + '_z' + str(mp.d_dm1_dm2) \
+ '_IWA' + str(mp.Fend.corr.Rin) + '_OWA' + str(mp.Fend.corr.Rout) + '_' + str(mp.Nsbp) + 'lams' \
+ str(round(1e9*mp.lambda0)) + 'nm_BW' + str(mp.fracBW*100) + '_' + mp.controller

print(mp.runLabel)


# Step 5: Perform the Wavefront Sensing and Control

out = falco.setup.flesh_out_workspace(mp)

falco.wfsc.loop(mp, out)
