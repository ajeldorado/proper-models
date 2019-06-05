#   Copyright 2019 California Institute of Technology
# ------------------------------------------------------------------

import shutil
import wfirst_phaseb_proper

def copy_examples_here( ):

    # copy Phase B PROPER examples into local directory

    files = ['run_flatten.py', 'run_hlc_erkin.py', 'run_hlc_input_fields.py', 
            'run_hlc.py', 'run_spc_ifs.py', 'run_spc_wide.py' ]

    for f in files:
        filename = wfirst_phaseb_proper.lib_dir + '/examples/' + f
        try:
            print( "Copying " + f + " to current directory" )
            shutil.copy( filename, './.' )
        except IOError as e:
            raise IOError( "Unable to copy prescription to current directory. %s" % e )

