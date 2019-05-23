import proper
import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
from trim import trim
import time

def run_hlc_input_fields():

    n = 512                 # output image dimension
    final_sampling = 0.1    # output sampling in lam0/D
    lam = np.array([ 0.54625, 0.5534375, 0.560625, 0.5678125, 0.575, 0.5821875, 0.589375, 0.5965625, 0.60375 ])
    nlam = 9

    # compute coronagraphic field with full prescription

    print 'Computing field with full prescription'
    t1 = time.time()
    (fields, sampling) = proper.prop_run_multi('wfirst_phaseb', lam, n, QUIET=True, \
        PASSVALUE={'use_errors':1,'polaxis':-2,'zindex':[4],'zval_m':[0.19e-9],'final_sampling_lam0':final_sampling} )
    t2 = time.time()
    print 'Time (sec) for full prescription = ', t2-t1
    images = np.abs(fields)**2
    image_full = np.sum( images, 0 ) / nlam

    # write FPM exit pupil fields with full prescription

    print 'Writing FPM exit pupil fields using full prescription'
    (fields, sampling) = proper.prop_run_multi('wfirst_phaseb', lam, n, QUIET=True, \
        PASSVALUE={'use_errors':1,'polaxis':-2,'zindex':[4],'zval_m':[0.19e-9],'use_hlc_dm_patterns':0,'use_fpm':0,'end_at_fpm_exit_pupil':1, \
    'output_field_rootname':'hlc_input_field'} )

    # compute coronagraphic field with compact prescription and input fields

    print 'Computing field with compact prescription and pre-computed input fields'
    t1 = time.time()
    (fields, sampling) = proper.prop_run_multi('wfirst_phaseb_compact', lam, n, QUIET=True, \
        PASSVALUE={'input_field_rootname':'hlc_input_field','polaxis':-2,'final_sampling_lam0':final_sampling} )
    t2 = time.time()
    print 'Time (sec) for full prescription = ', t2-t1
    images = np.abs(fields)**2
    image = np.sum( images, 0 ) / nlam

    # move source to 7 lam/D

    print 'Computing offset source'
    (fields, sampling) = proper.prop_run_multi('wfirst_phaseb_compact', lam, n, QUIET=True, \
        PASSVALUE={'source_x_offset':7.0,'final_sampling_lam0':final_sampling} )
    psfs = np.abs(fields)**2
    psf = np.sum( psfs, 0 ) / nlam
    max_psf = np.max(psf)

    ni_full = image_full / max_psf
    ni = image / max_psf

    fig, ax = plt.subplots( nrows=1, ncols=2 )

    im = ax[0].imshow(ni_full, norm=LogNorm(vmin=1e-7,vmax=1e-2))
    ax[0].set_title('Full prescription')
    fig.colorbar( im, ax=ax[0], shrink=0.5 ) 

    im = ax[1].imshow(ni, norm=LogNorm(vmin=1e-7,vmax=1e-2))
    ax[1].set_title('Compact prescription')
    fig.colorbar( im, ax=ax[1], shrink=0.5 ) 

    plt.show()
    
if __name__ == '__main__':
    run_hlc_input_fields()
