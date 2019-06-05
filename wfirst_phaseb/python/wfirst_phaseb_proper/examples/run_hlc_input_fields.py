import proper
import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle
import wfirst_phaseb_proper
from wfirst_phaseb_proper import trim
import time

def run_hlc_input_fields():

    nlam = 7
    lam0 = 0.575
    bandwidth = 0.1
    minlam = lam0 * (1 - bandwidth/2)
    maxlam = lam0 * (1 + bandwidth/2)
    lam_array = np.linspace( minlam, maxlam, nlam )

    n = 256                 # output image dimension
    final_sampling = 0.1    # output sampling in lam0/D

    # compute coronagraphic field with full prescription

    print( 'Computing field with full prescription' )
    t1 = time.time()
    (fields, sampling) = proper.prop_run_multi('wfirst_phaseb', lam_array, n, QUIET=True, \
        PASSVALUE={'use_errors':1,'polaxis':-2,'zindex':[4],'zval_m':[0.19e-9],'final_sampling_lam0':final_sampling} )
    t2 = time.time()
    print( 'Time (sec) for full prescription = ' + str(t2-t1) )
    images = np.abs(fields)**2
    image_full = np.sum( images, 0 ) / nlam

    # write FPM exit pupil fields with full prescription

    print( 'Writing FPM exit pupil fields using full prescription' )
    (fields, sampling) = proper.prop_run_multi('wfirst_phaseb', lam_array, n, QUIET=True, \
        PASSVALUE={'use_errors':1,'polaxis':-2,'zindex':[4],'zval_m':[0.19e-9],'use_hlc_dm_patterns':0,'use_fpm':0,'end_at_fpm_exit_pupil':1, \
        'output_field_rootname':'hlc_input_field'} )

    # compute coronagraphic field with compact prescription and input fields

    print( 'Computing field with compact prescription and pre-computed input fields' )
    t1 = time.time()
    (fields, sampling) = proper.prop_run_multi('wfirst_phaseb_compact', lam_array, n, QUIET=True, \
        PASSVALUE={'input_field_rootname':'hlc_input_field','polaxis':-2,'final_sampling_lam0':final_sampling} )
    t2 = time.time()
    print( 'Time (sec) for compact prescription with input fields = ' + str(t2-t1) )
    images = np.abs(fields)**2
    image = np.sum( images, 0 ) / nlam

    # move source to 7 lam/D

    print( 'Computing offset source to compute NI' )
    (fields, sampling) = proper.prop_run_multi('wfirst_phaseb_compact', lam_array, n, QUIET=True, \
        PASSVALUE={'source_x_offset':7.0,'final_sampling_lam0':final_sampling} )
    psfs = np.abs(fields)**2
    psf = np.sum( psfs, 0 ) / nlam
    max_psf = np.max(psf)

    ni_full = image_full / max_psf
    ni = image / max_psf

    fig, ax = plt.subplots( nrows=1, ncols=2, figsize=(9,4) )

    im = ax[0].imshow(ni_full, norm=LogNorm(vmin=1e-7,vmax=1e-2), cmap=plt.get_cmap('jet'))
    circ = Circle((n/2,n/2),3/final_sampling,edgecolor='white', facecolor='none')
    ax[0].add_patch(circ)
    circ = Circle((n/2,n/2),9/final_sampling,edgecolor='white', facecolor='none')
    ax[0].add_patch(circ)
    ax[0].set_title('Full prescription')
    fig.colorbar( im, ax=ax[0], shrink=0.5 ) 

    im = ax[1].imshow(ni, norm=LogNorm(vmin=1e-7,vmax=1e-2), cmap=plt.get_cmap('jet'))
    circ = Circle((n/2,n/2),3/final_sampling,edgecolor='white', facecolor='none')
    ax[1].add_patch(circ)
    circ = Circle((n/2,n/2),9/final_sampling,edgecolor='white', facecolor='none')
    ax[1].add_patch(circ)
    ax[1].set_title('Compact prescription')
    fig.colorbar( im, ax=ax[1], shrink=0.5 ) 

    plt.show()
    
if __name__ == '__main__':
    run_hlc_input_fields()
