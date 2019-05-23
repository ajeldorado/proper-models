import proper
import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
from trim import trim

def run_flatten():

    lam_array = np.array([ 0.54625, 0.5534375, 0.560625, 0.5678125, 0.575, 0.5821875, 0.589375, 0.5965625, 0.60375 ])
    nlam = 9

    final_sampling = 0.1

    # compute field before flattening 

    fields, sampling = proper.prop_run_multi( 'wfirst_phaseb', lam_array, 512, QUIET=True, \
        PASSVALUE={'final_sampling_lam0':final_sampling} )
    images = np.abs(fields)**2
    image_before = np.sum(images,0) / nlam

    # compute field after flattening 

    dm1 = proper.prop_fits_read('dm1_flatten.fits')

    fields, sampling = proper.prop_run_multi( 'wfirst_phaseb', lam_array, 512, QUIET=True, \
        PASSVALUE={'final_sampling_lam0':final_sampling,'use_dm1':1,'dm1_m':dm1} )
    images = np.abs(fields)**2
    image_after = np.sum(images,0) / nlam

    # move source by 7.0 lam/D

    fields, sampling = proper.prop_run_multi( 'wfirst_phaseb_compact', lam_array, 512, QUIET=True, \
        PASSVALUE={'source_x_offset':7.0,'final_sampling_lam0':final_sampling} )
    images = np.abs(fields)**2
    psf = np.sum(images,0) / nlam
    max_psf = np.max(psf)

    ni_before = image_before / max_psf
    ni_after = image_after / max_psf

    fig, ax = plt.subplots( nrows=1, ncols=2 )

    im = ax[0].imshow(ni_before, norm=LogNorm(vmin=1e-7,vmax=1e-2))
    ax[0].set_title('Before')
    fig.colorbar( im, ax=ax[0], shrink=0.5 ) 

    im = ax[1].imshow(ni_after, norm=LogNorm(vmin=1e-7,vmax=1e-2))
    ax[1].set_title('After')
    fig.colorbar( im, ax=ax[1], shrink=0.5 ) 

    plt.show()
    
if __name__ == '__main__':
    run_flatten()
