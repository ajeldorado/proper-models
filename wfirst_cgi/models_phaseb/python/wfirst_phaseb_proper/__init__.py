#   Copyright 2019 California Institute of Technology
# ------------------------------------------------------------------

import os
import os.path as _osp

lib_dir = _osp.abspath(_osp.dirname(__file__))

__version__ = '1.7'

# from .wfirst_phaseb import wfirst_phaseb
# from .wfirst_phaseb_compact import wfirst_phaseb_compact
from .trim import trim
from .polmap import polmap
from .ffts import ffts
from .mft2 import mft2
from .copy_here import copy_here
from .copy_examples_here import copy_examples_here
from .set_data_dir import set_data_dir

map_dir =  '/maps/'
polfile = '/pol/new_toma'

data_dir ="/Users/ajriggs/Repos/proper-models/wfirst_cgi/data_phaseb"
