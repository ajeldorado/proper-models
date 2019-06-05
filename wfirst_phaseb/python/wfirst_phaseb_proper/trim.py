#   Copyright 2019 California Institute of Technology
# ------------------------------------------------------------------

import numpy as np

def trim( input_image, output_dim ):

    input_dim = input_image.shape[1]

    if input_dim == output_dim:
        return input_image
    elif output_dim < input_dim:
        x1 = input_dim // 2 - output_dim // 2
        x2 = x1 + output_dim
        output_image = input_image[x1:x2,x1:x2].copy()
    else:
        output_image = np.zeros((output_dim,output_dim), dtype=input_image.dtype)
        x1 = output_dim // 2 - input_dim // 2
        x2 = x1 + input_dim
        output_image[x1:x2,x1:x2] = input_image

    return output_image

