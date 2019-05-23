import numpy as np

def mft2( field_in, dout, D, nout, direction, xoffset=0, yoffset=0, xc=0, yc=0 ):

    nfield_in = field_in.shape[1] 
    nfield_out = int(nout)
 
    x = (np.arange(nfield_in) - nfield_in//2 - xc) 
    y = (np.arange(nfield_in) - nfield_in//2 - yc) 

    u = (np.arange(nfield_out) - nfield_out//2 - xoffset/dout) * (dout/D)
    v = (np.arange(nfield_out) - nfield_out//2 - yoffset/dout) * (dout/D)

    xu = np.outer(x, u)
    yv = np.outer(y, v)

    if direction == -1:
        expxu = dout/D * np.exp(-2.0 * np.pi * -1j * xu)
        expyv = np.exp(-2.0 * np.pi * -1j * yv).T
    else:
        expxu = dout/D * np.exp(-2.0 * np.pi * 1j * xu)
        expyv = np.exp(-2.0 * np.pi * 1j * yv).T

    t1 = np.dot(expyv, field_in)
    t2 = np.dot(t1, expxu)

    return t2
