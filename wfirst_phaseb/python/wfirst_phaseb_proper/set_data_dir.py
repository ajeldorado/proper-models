#   Copyright 2019 California Institute of Technology
# ------------------------------------------------------------------

import os
import wfirst_phaseb_proper

def set_data_dir( ):

    # set the data_dir variable in the package __init__.py file to the current directory

    cwd = os.getcwd()

    with open( wfirst_phaseb_proper.lib_dir+"/__init__.py", "a" ) as f:
        f.write( 'data_dir ="' + cwd + '"\n' )
        
