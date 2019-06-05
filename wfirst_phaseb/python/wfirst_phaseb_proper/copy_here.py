import shutil
import wfirst_phaseb_proper

def copy_here( ):

    # copy Phase B PROPER prescriptions from wfirst_phaseb_proper package into local directory

    prescription_file = wfirst_phaseb_proper.lib_dir + '/wfirst_phaseb.py'
    try:
        shutil.copy( prescription_file, './.' )
    except IOError as e:
        raise IOError( "Unable to copy prescription to current directory. %s" % e )

    prescription_file = wfirst_phaseb_proper.lib_dir + '/wfirst_phaseb_compact.py'
    try:
        shutil.copy( prescription_file, './.' )
    except IOError as e:
        raise IOError( "Unable to copy prescription to current directory. %s" % e )

