import proper
import numpy as np

def ffts( wavefront, direction ):
    n = wavefront.shape[0]  # assumed to be power of 2
    wavefront = np.roll( np.roll(wavefront, -n//2, 0), -n//2, 1 ) # shift to corner
    if proper.use_fftw:
        if direction == -1:
            new_wavefront = proper.prop_fftw( wavefront, directionFFTW='FFTW_FORWARD' ) / np.size(wavefront)
        else:
            new_wavefront = proper.prop_fftw( wavefront, directionFFTW='FFTW_BACKWARD' ) * np.size(wavefront)
    else:
        if direction == -1:
            new_wavefront = np.fft.fft2(wavefront) / np.size(wavefront)
        else:
            new_wavefront = np.fft.ifft2(wavefront) * np.size(wavefront)
    
    new_wavefront = np.roll( np.roll(new_wavefront, n//2, 0), n//2, 1 ) # shift to center 

    return new_wavefront


