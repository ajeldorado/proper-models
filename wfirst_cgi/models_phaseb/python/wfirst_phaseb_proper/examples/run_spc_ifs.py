#   Copyright 2019 California Institute of Technology
# ------------------------------------------------------------------

import proper
import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle
from wfirst_phaseb_proper import trim

def run_spc_ifs():

    nlam = 9
    lam0 = 0.73
    bandwidth = 0.18
    lam_min = lam0 * (1 - bandwidth/2)
    lam_max = lam0 * (1 + bandwidth/2)
    lam = np.linspace( lam_min, lam_max, nlam )

    n = 256                 # output image dimension
    final_sampling = 0.1    # output sampling in lam0/D

    # compute coronagraphic image

    print( "Computing coronagraphic field..." )
    (fields, sampling) = proper.prop_run_multi('wfirst_phaseb', lam, n, QUIET=True, \
        PASSVALUE={'cor_type':'spc-ifs_long','lam0':lam0,'use_errors':0,'zindex':[4],'zval_m':[0.19e-9],'final_sampling_lam0':final_sampling} )
    images = np.abs(fields)**2
    image = np.sum( images, 0 ) / nlam

    # move source to 7 lam/D

    print( "Computing offset source to compute NI..." )
    (fields, sampling) = proper.prop_run_multi('wfirst_phaseb_compact', lam, n, QUIET=True, \
        PASSVALUE={'cor_type':'spc-ifs_long','lam0':lam0,'source_x_offset':7.0,'final_sampling_lam0':final_sampling} )
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
    run_spc_ifs()
