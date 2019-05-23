import proper
import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
from trim import trim

def run_hlc_s_vs_p():

    n = 512                 # output image dimension
    final_sampling = 0.1    # output sampling in lam0/D
    lam = np.array([ 0.54625, 0.5534375, 0.560625, 0.5678125, 0.575, 0.5821875, 0.589375, 0.5965625, 0.60375 ])
    nlam = 9

    # compute coronagraphic field

    (fields, sampling) = proper.prop_run_multi('wfirst_phaseb', lam, n, QUIET=True, \
        PASSVALUE={'fpm_axis':'s','cor_type':'hlc','use_errors':0,'zindex':[4],'zval_m':[0.19e-9],'final_sampling_lam0':final_sampling} )
    images = np.abs(fields)**2
    simage = np.sum( images, 0 ) / nlam

    (fields, sampling) = proper.prop_run_multi('wfirst_phaseb', lam, n, QUIET=True, \
        PASSVALUE={'fpm_axis':'p','cor_type':'hlc','use_errors':0,'zindex':[4],'zval_m':[0.19e-9],'final_sampling_lam0':final_sampling} )
    images = np.abs(fields)**2
    pimage = np.sum( images, 0 ) / nlam

    # move source to 7 lam/D

    (fields, sampling) = proper.prop_run_multi('wfirst_phaseb_compact', lam, n, QUIET=True, \
        PASSVALUE={'cor_type':'hlc','source_x_offset':7.0,'final_sampling_lam0':final_sampling} )
    psfs = np.abs(fields)**2
    psf = np.sum( psfs, 0 ) / nlam

    sni = simage / np.max(psf)
    pni = pimage / np.max(psf)

    fig, ax  = plt.subplots( nrows=1, ncols=2 )

    im = ax[0].imshow(sni, norm=LogNorm(vmin=1e-10,vmax=1e-7))
    ax[0].set_title('S axis')
    fig.colorbar( im, ax=ax[0], shrink=0.5 )
 
    im = ax[1].imshow(pni, norm=LogNorm(vmin=1e-10,vmax=1e-7))
    ax[1].set_title('P axis')
    fig.colorbar( im, ax=ax[1], shrink=0.5 )
 
    plt.show()
    
if __name__ == '__main__':
    run_hlc_s_vs_p()
