#   Copyright 2019 California Institute of Technology
# ------------------------------------------------------------------

import proper
import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle
import wfirst_phaseb_proper
from wfirst_phaseb_proper import trim

def run_hlc():

    nlam = 7
    lam0 = 0.575
    bandwidth = 0.1
    minlam = lam0 * (1 - bandwidth/2)
    maxlam = lam0 * (1 + bandwidth/2)
    lam_array = np.linspace( minlam, maxlam, nlam )

    n = 256                 # output image dimension (must be power of 2)
    final_sampling = 0.1    # output sampling in lam0/D

    # compute unaberrated coronagraphic field using Dwight's DM wavefront maps

    print( "Computing unaberrated coronagraphic field using DM wavefront map..." )
    (fields, sampling) = proper.prop_run_multi('wfirst_phaseb', lam_array, n, QUIET=True, \
        PASSVALUE={'cor_type':'hlc','use_errors':0,'zindex':[4],'zval_m':[0.19e-9],
        'use_hlc_dm_patterns':1,'final_sampling_lam0':final_sampling} )
    images = np.abs(fields)**2
    image = np.sum( images, 0 ) / nlam

    # compute unaberrated coronagraphic field using DM actuator settings

    print( "Computing unaberrated coronagraphic field using DM actuator pistons..." )
    dm1 = proper.prop_fits_read( wfirst_phaseb_proper.lib_dir + '/examples/hlc_dm1.fits' )
    dm2 = proper.prop_fits_read( wfirst_phaseb_proper.lib_dir + '/examples/hlc_dm2.fits' )
    (fields, sampling) = proper.prop_run_multi('wfirst_phaseb', lam_array, n, QUIET=True, \
        PASSVALUE={'cor_type':'hlc','use_errors':0,'zindex':[4],'zval_m':[0.19e-9],
        'use_hlc_dm_patterns':0,'final_sampling_lam0':final_sampling,
        'use_dm1':1, 'dm1_m':dm1, 'use_dm2':1, 'dm2_m':dm2} )
    images = np.abs(fields)**2
    image_dm = np.sum( images, 0 ) / nlam

    # compute aberrated coronagraphic field using DM actuator settings

    print( "Computing aberrated coronagraphic field using DM actuator pistons..." )
    dm1 = proper.prop_fits_read( wfirst_phaseb_proper.lib_dir + '/examples/hlc_with_aberrations_dm1.fits' )
    dm2 = proper.prop_fits_read( wfirst_phaseb_proper.lib_dir + '/examples/hlc_with_aberrations_dm2.fits' )
    (fields, sampling) = proper.prop_run_multi('wfirst_phaseb', lam_array, n, QUIET=True, \
        PASSVALUE={'cor_type':'hlc', 'use_errors':1, 'polaxis':10, 
        'use_hlc_dm_patterns':0, 'final_sampling_lam0':final_sampling,
        'use_dm1':1, 'dm1_m':dm1, 'use_dm2':1, 'dm2_m':dm2} )
    images = np.abs(fields)**2
    image_ab = np.sum( images, 0 ) / nlam

    # move source to 7 lam/D

    print( "Computing offset source to compute NI..." )
    (fields, sampling) = proper.prop_run_multi('wfirst_phaseb_compact', lam_array, n, QUIET=True, \
        PASSVALUE={'cor_type':'hlc','source_x_offset':7.0,'use_hlc_dm_patterns':1,'final_sampling_lam0':final_sampling} )
    psfs = np.abs(fields)**2
    psf = np.sum( psfs, 0 ) / nlam

    max_psf = np.max(psf)

    ni = image / max_psf
    ni_dm = image_dm / max_psf
    ni_ab = image_ab / max_psf

    fig, ax  = plt.subplots( nrows=1, ncols=3, figsize=(11,4) )

    im = ax[0].imshow(ni, norm=LogNorm(vmin=1e-10,vmax=1e-7), cmap=plt.get_cmap('jet'))
    circ_in = Circle((n/2,n/2),3/final_sampling,edgecolor='white', facecolor='none')
    circ_out = Circle((n/2,n/2),9/final_sampling,edgecolor='white', facecolor='none')
    ax[0].add_patch(circ_in)
    ax[0].add_patch(circ_out)
    ax[0].set_title('Unaberrated,\nUsing wavefront maps')
    fig.colorbar(im, ax=ax[0], shrink=0.5) 

    im = ax[1].imshow(ni_dm, norm=LogNorm(vmin=1e-10,vmax=1e-7), cmap=plt.get_cmap('jet'))
    circ_in = Circle((n/2,n/2),3/final_sampling,edgecolor='white', facecolor='none')
    circ_out = Circle((n/2,n/2),9/final_sampling,edgecolor='white', facecolor='none')
    ax[1].add_patch(circ_in)
    ax[1].add_patch(circ_out)
    ax[1].set_title('Unaberrated,\nUsing DM piston maps')
    fig.colorbar(im, ax=ax[1], shrink=0.5) 

    im = ax[2].imshow(ni_ab, norm=LogNorm(vmin=1e-10,vmax=1e-7), cmap=plt.get_cmap('jet'))
    circ_in = Circle((n/2,n/2),3/final_sampling,edgecolor='white', facecolor='none')
    circ_out = Circle((n/2,n/2),9/final_sampling,edgecolor='white', facecolor='none')
    ax[2].add_patch(circ_in)
    ax[2].add_patch(circ_out)
    ax[2].set_title('Aberrated,\nUsing DM piston maps')
    fig.colorbar(im, ax=ax[2], shrink=0.5) 

    plt.show()
    
if __name__ == '__main__':
    run_hlc()
