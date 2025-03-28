# fix SDS band codes
import os
import re
import shutil
from obspy import read
import matplotlib.pyplot as plt
import sys
import time
from lib.libseisGT_old3 import fix_trace_id, check_read_write

# many functions have been moved to SDS.py in flovopy


def find_soufrierehills_lib(basedir='SoufriereHillsVolcano'):
    """Find the 'lib' directory under 'SoufriereHillsVolcano' inside the 'Developer' folder in HOME."""
    
    home_dir = os.path.expanduser("~")  # Get home directory
    developer_dir = os.path.join(home_dir, "Developer")

    if not os.path.isdir(developer_dir):
        return "Developer directory not found in HOME."

    # Look for 'SoufriereHills' within 1 or 2 levels under Developer
    soufrierehills_path = None
    for root, dirs, _ in os.walk(developer_dir):
        if basedir in dirs:
            soufrierehills_path = os.path.join(root, basedir)
            break  # Stop searching once found

    if not soufrierehills_path:
        return "SoufriereHills directory not found under Developer."

    # Construct the path to 'lib' directory
    lib_path = os.path.join(soufrierehills_path, "lib")

    if os.path.isdir(lib_path):
        return lib_path  # Return the valid path
    else:
        return "lib directory not found under SoufriereHills."


def main(sds_directory, networks=None, write=True, backup=True):

    print('\n', f'Note that networks={networks}')
    if write:
        print('and CHANGES WILL BE MADE')
    else:
        print('changes will not be made - this is a dry-run\n- change value of "write" variable to True if you want changes to be made')

    print('\nFixing filenames like .miniseed or .part')
    fix_sds_filenames(sds_directory, networks=networks, write=write)

    print('\nFixing IDs based on sampling rate and location code')
    fix_sds_ids(sds_directory, networks=networks, write=write)

    print('\nMoving files to correct directory') # e.g. if fix_sds_ids moved a HHZ.D file to DHZ.D, but it does not change the directory
    move_files_to_correct_directory(sds_directory, write=write)

    print('\nRemoving any .DS_Store files and any empty directories')
    remove_ds_store_files_and_empty_dirs(sds_directory)


##################################### end of dealing with SDS archive ##################################

##################################### now deal with loose files not in SDS archive yet #################



if __name__ == "__main__":
    # Example usage
    lib_dir = find_soufrierehills_lib()
    print('lib_dir=',lib_dir)
    sys.path.append(lib_dir)
    import libMVO
    sds_directory = "/data/SDS"  # Change this to your SDS archive path
    BACKUPDIR = '/data/backups'
    NOTMERGEDDIR = '/data/notmerged'
    os.makedirs(BACKUPDIR, exist_ok=True)
    os.makedirs(NOTMERGEDDIR, exist_ok=True)
    networks = ['FL']
    write=True
    backup=write
    main(sds_directory, networks=networks, write=write, backup=backup)

    source_directory = "/data/KSC/eventminiseedfiles"  # Change this to your flat directory
    move_files_to_sds_structure(source_directory, sds_directory, write=write, backup=backup)
