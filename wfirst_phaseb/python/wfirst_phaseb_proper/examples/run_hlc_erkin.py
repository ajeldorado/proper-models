import proper
import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle
import wfirst_phaseb_proper
from wfirst_phaseb_proper import trim

def run_hlc_erkin():

    nlam = 7
    lam0 = 0.575
    bandwidth = 0.1
    minlam = lam0 * (1 - bandwidth/2)
    maxlam = lam0 * (1 + bandwidth/2)
    lam_array = np.linspace( minlam, maxlam, nlam )

    n = 256                 # output image dimension
    final_sampling = 0.1    # output sampling in lam0/D

    # compute coronagraphic field

    print( "Computing coronagraphic field..." )
    (fields, sampling) = proper.prop_run_multi('wfirst_phaseb', lam_array, n, QUIET=True, \
        PASSVALUE={'cor_type':'hlc_erkin','use_errors':0,'zindex':[4],'zval_m':[0.19e-9],
        'use_hlc_dm_patterns':1, 'final_sampling_lam0':final_sampling} )
    images = np.abs(fields)**2
    image = np.sum( images, 0 ) / nlam

    # move source to 7 lam/D

    print( "Computing offset source to compute NI..." )
    (fields, sampling) = proper.prop_run_multi('wfirst_phaseb_compact', lam_array, n, QUIET=True, \
        PASSVALUE={'cor_type':'hlc_erkin','source_x_offset':7.0,'final_sampling_lam0':final_sampling,
        'use_hlc_dm_patterns':1} )
    psfs = np.abs(fields)**2
    psf = np.sum( psfs, 0 ) / nlam

    ni = image / np.max(psf)

    fig, ax  = plt.subplots( figsize=(8,8) )
    im = ax.imshow(ni, norm=LogNorm(vmin=1e-10,vmax=1e-7), cmap=plt.get_cmap('jet'))
    circ = Circle((n/2,n/2),3/final_sampling,edgecolor='white', facecolor='none')
    ax.add_patch(circ)
    circ = Circle((n/2,n/2),9/final_sampling,edgecolor='white', facecolor='none')
    ax.add_patch(circ)
    fig.colorbar(im, shrink=0.5) 
    plt.show()
    
if __name__ == '__main__':
    run_hlc_erkin()
