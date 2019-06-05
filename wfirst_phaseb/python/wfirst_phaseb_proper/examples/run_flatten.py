import proper
import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle
import wfirst_phaseb_proper
from wfirst_phaseb_proper import trim

def run_flatten():

    nlam = 7
    lam0 = 0.575
    bandwidth = 0.1
    minlam = lam0 * (1 - bandwidth/2)
    maxlam = lam0 * (1 + bandwidth/2)
    lam_array = np.linspace( minlam, maxlam, nlam )

    n = 256
    final_sampling = 0.1

    # compute field before flattening (use HLC DM WFE maps) 

    print( "Computing field before flattening..." )
    fields, sampling = proper.prop_run_multi( 'wfirst_phaseb', lam_array, n, QUIET=True, \
        PASSVALUE={'final_sampling_lam0':final_sampling, 'use_hlc_dm_patterns':1, 'use_errors':1, 'polaxis':10} )
    images = np.abs(fields)**2
    image_before = np.sum(images,0) / nlam

    # compute field after flattening 

    dm1 = proper.prop_fits_read( wfirst_phaseb_proper.lib_dir + '/examples/errors_polaxis10_dm.fits' )

    print( "Computing field after flattening..." )
    fields, sampling = proper.prop_run_multi( 'wfirst_phaseb', lam_array, n, QUIET=True, \
        PASSVALUE={'final_sampling_lam0':final_sampling, 'use_dm1':1, 'dm1_m':dm1, 'use_hlc_dm_patterns':1, 'use_errors':1, 'polaxis':10} )
    images = np.abs(fields)**2
    image_after = np.sum(images,0) / nlam

    # move source by 7.0 lam/D

    print( "Computing offset source to compute NI..." )
    fields, sampling = proper.prop_run_multi( 'wfirst_phaseb_compact', lam_array, n, QUIET=True, \
        PASSVALUE={'source_x_offset':7.0,'final_sampling_lam0':final_sampling,'use_hlc_dm_patterns':1} )
    images = np.abs(fields)**2
    psf = np.sum(images,0) / nlam
    max_psf = np.max(psf)

    ni_before = image_before / max_psf
    ni_after = image_after / max_psf

    fig, ax = plt.subplots( nrows=1, ncols=2, figsize=(9,4) )

    im = ax[0].imshow(ni_before, norm=LogNorm(vmin=1e-7,vmax=1e-2), cmap=plt.get_cmap('jet'))
    circ = Circle((n/2,n/2),3/final_sampling,edgecolor='white', facecolor='none')
    ax[0].add_patch(circ)
    circ = Circle((n/2,n/2),9/final_sampling,edgecolor='white', facecolor='none')
    ax[0].add_patch(circ)
    ax[0].set_title('Before')
    fig.colorbar( im, ax=ax[0], shrink=0.5 ) 

    im = ax[1].imshow(ni_after, norm=LogNorm(vmin=1e-7,vmax=1e-2), cmap=plt.get_cmap('jet'))
    circ = Circle((n/2,n/2),3/final_sampling,edgecolor='white', facecolor='none')
    ax[1].add_patch(circ)
    circ = Circle((n/2,n/2),9/final_sampling,edgecolor='white', facecolor='none')
    ax[1].add_patch(circ)
    ax[1].set_title('After')
    fig.colorbar( im, ax=ax[1], shrink=0.5 ) 

    plt.show()
    
if __name__ == '__main__':
    run_flatten()
