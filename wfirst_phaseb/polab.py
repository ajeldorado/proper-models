import numpy as np
from scipy.interpolate import interp1d
import math
import proper
from trim import trim

# wavefront: current wavefront structure
# polfile: rootname of file containing polarization coefficients
# pupil_diam_pix: diameter of pupil in pixels
# condition: polarization circumstance:
#        -2: -45 deg in, Y out
#        -1: -45 deg in, X out
#         1: +45 deg in, X out
#         2: +45 deg in, Y out
#         5: X polarization (mean of +/-45 deg in, X out)
#         6: Y polarization (mean of +/-45 deg in, X out)
#        10: All polarizations (mean of +/-45 deg in, X&Y out)
#    NOTE: the mean conditions (5,6,10) should only be used for sensing;
#    contrast evaluation must be done by computing each in/out condition separately

def polmap( wavefront, polfile, pupil_diam_pix, condition, MUF=1.0 ):
    n = proper.prop_get_gridsize( wavefront )
    lambda_m = proper.prop_get_wavelength(wavefront)

    if condition <= 2:
        (amp, pha) = polab( polfile, lambda_m, pupil_diam_pix, condition )
    elif condition == 5:
        (amp_m45_x, pha_m45_x) = polab( polfile, lambda_m, pupil_diam_pix, -1 )
        (amp_p45_x, pha_p45_x) = polab( polfile, lambda_m, pupil_diam_pix, +1 )
        amp = (amp_m45_x + amp_p45_x) / 2
        pha = (pha_m45_x + pha_p45_x) / 2
    elif condition == 6:
        (amp_m45_y, pha_m45_y) = polab( polfile, lambda_m, pupil_diam_pix, -2 )
        (amp_p45_y, pha_p45_y) = polab( polfile, lambda_m, pupil_diam_pix, +2 )
        amp = (amp_m45_y + amp_p45_y) / 2
        pha = (pha_m45_y + pha_p45_y) / 2
    elif condition == 10:
        (amp_m45_x, pha_m45_x) = polab( polfile, lambda_m, pupil_diam_pix, -1 )
        (amp_p45_x, pha_p45_x) = polab( polfile, lambda_m, pupil_diam_pix, +1 )
        (amp_m45_y, pha_m45_y) = polab( polfile, lambda_m, pupil_diam_pix, -2 )
        (amp_p45_y, pha_p45_y) = polab( polfile, lambda_m, pupil_diam_pix, +2 )
        amp = (amp_m45_x + amp_p45_x + amp_m45_y + amp_p45_y) / 4
        pha = (pha_m45_x + pha_p45_x + pha_m45_y + pha_p45_y) / 4
    else:
        print 'POLMAP: unmatched condition = ', condition
        return

    proper.prop_multiply( wavefront, trim(amp,n) ) 
    proper.prop_add_phase( wavefront, trim(MUF*pha,n) )

    return 



# polfile: rootname of file containing polarization coefficients
# lambda_m: wavelength in meters
# pupil_diam_pix: diameter of pupil in pixels
# condition: polarization circumstance:
#        -2: -45 deg in, Y out
#        -1: -45 deg in, X out
#         1: +45 deg in, X out
#         2: +45 deg in, Y out
# amp, pha: returned aberration maps (pha is WFE in meters)

def polab( polfile, lambda_m, pupil_diam_pix, condition ):
    if abs(condition) == 1:
        dir_out = 0 
    else:
        dir_out = 1         # dir_out: output polarization (1=X, 2=Y)
    if condition < 0:
        dir_in = 0 
    else:
        dir_in = 1          # dir_in: input polarization (negative=-45d, positive=+45d)

    # zernike coefficient files are [nzer, nlam, ndir_in, ndir_out]
    #    nzer = 22 (number of zernikes)
    #    nlam = 6 or 11 (450 - 950 nm in 100 or 50 nm steps)
    #    ndir_in = 2 (input polarization direction, 0=-45 deg, 1=+45 deg)
    #    ndir_out = 2 (output polarization direction, 0=X, 1=Y)

    zamp_array = proper.prop_fits_read( polfile+'_amp.fits' )
    zpha_array = proper.prop_fits_read( polfile+'_pha.fits' )
    nlam = zamp_array.shape[2]
    if nlam == 6:
        lam_array_m = (np.arange(6) * 100 + 450) * 1.0e-9 
    else:
        lam_array_m = (np.arange(11) * 50 + 450) * 1.0e-9

    # interpolate to get zernikes at specified wavelength

    zamp = np.zeros([22])
    zpha = np.zeros([22])

    for iz in range(0, 22):
        famp = interp1d( lam_array_m, zamp_array[dir_out, dir_in, :, iz], kind='cubic' )
        fpha = interp1d( lam_array_m, zpha_array[dir_out, dir_in, :, iz], kind='cubic' )
        lam = lambda_m
        if lam < 0.45e-6: lam = 0.45e-6
        if lam > 0.95e-6: lam = 0.95e-6
        zamp[iz] = famp( lambda_m )
        zpha[iz] = fpha( lambda_m )

    n = int(round(pupil_diam_pix * 1.1))
    n = (n // 2) * 2     # force even 
    x = (np.arange(n) - n//2) / (pupil_diam_pix/2.0)

    amp = np.zeros([n,n])
    pha = np.zeros([n,n])

    for j in range(0, n):
        y = x[j]
        r2 = x**2 + y**2
        r = np.sqrt(r2)
        r3 = r**3
        r4 = r**4
        r5 = r**5
        r6 = r**6
        t = np.arctan2(y,x)

        for itype in range(0,2):        # 0 = amp, 1 = phase
            map = np.zeros([n])

            if itype == 0:
                z = zamp
                map = map + z[0]    # include piston if amplitude map
            else:
                z = zpha

            map = map + z[1] * 2 * x                # x tilt
            map = map + z[2] * 2 * y                # y tilt
            map = map + z[3] * np.sqrt(3) * (2*r2 - 1)            # focus
            map = map + z[4] * np.sqrt(6) * r2 * np.sin(2*t)        # 45 deg astig
            map = map + z[5] * np.sqrt(6) * r2 * np.cos(2*t)        # 0 deg astig
            map = map + z[6] * np.sqrt(8) * (3*r3 - 2*r) * np.sin(t)      # y coma
            map = map + z[7] * np.sqrt(8) * (3*r3 - 2*r) * np.cos(t)    # x coma
            map = map + z[8] * np.sqrt(8) * r3 * np.sin(3*t)        # y trefoil 
            map = map + z[9] * np.sqrt(8) * r3 * np.cos(3*t)        # x trefoil 
            map = map + z[10] * np.sqrt(5) * (6*r4 - 6*r2 + 1)        # spherical
            map = map + z[11] * np.sqrt(10) * (4*r4 - 3*r2) * np.cos(2*t)
            map = map + z[12] * np.sqrt(10) * (4*r4 - 3*r2) * np.sin(2*t)
            map = map + z[13] * np.sqrt(10) * r4 * np.cos(4*t)
            map = map + z[14] * np.sqrt(10) * r4 * np.sin(4*t)
            map = map + z[15] * np.sqrt(12) * (10*r5 - 12*r3 + 3*r) * np.cos(t)
            map = map + z[16] * np.sqrt(12) * (10*r5 - 12*r3 + 3*r) * np.sin(t)
            map = map + z[17] * np.sqrt(12) * (5*r5 - 4*r3) * np.cos(3*t)
            map = map + z[18] * np.sqrt(12) * (5*r5 - 4*r3) * np.sin(3*t)
            map = map + z[19] * np.sqrt(12) * r5 * np.cos(5*t)
            map = map + z[20] * np.sqrt(12) * r5 * np.sin(5*t)
            map = map + z[21] * np.sqrt(7) * (20*r6 - 30*r4 + 12*r2 - 1)

            if itype == 0:
                amp[j,:] = map 
            else:
                pha[j,:] = map

    return amp, pha

