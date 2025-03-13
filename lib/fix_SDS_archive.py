# fix SDS band codes
import os
import re
import shutil
from obspy import read
import matplotlib.pyplot as plt
import sys
import time
from lib.libseisGT_old3 import fix_trace_id, check_read_write

def get_processed_dirs(LOG_FILE):
    """Reads the log file and returns a set of already processed directories."""
    if os.path.exists(LOG_FILE):
        with open(LOG_FILE, "r") as f:
            return set(line.strip() for line in f)
    return set()

def save_processed_dir(directory, LOG_FILE):
    """Logs the processed directory to a file."""
    with open(LOG_FILE, "a") as f:
        f.write(directory + "\n")


def trace2correct_sdsfullpath(base_dir, trace):
    """Determine the correct SDS path for a trace based on its metadata."""
    net = trace.stats.network
    sta = trace.stats.station
    #loc = trace.stats.location if trace.stats.location else "--"
    loc = trace.stats.location
    if len(loc)==1:
        loc = '0' + loc
    chan = trace.stats.channel
    year = str(trace.stats.starttime.year)
    dtype = "D"  # Typically "D" for daily files in SDS
    
    day_of_year = str(trace.stats.starttime.julday).zfill(3)
    filename = f"{net}.{sta}.{loc}.{chan}.{dtype}.{year}.{day_of_year}"
    sds_subdir = os.path.join(year, net, sta, f"{chan}.{dtype}")

    return os.path.join(base_dir, sds_subdir, filename)



def movefile(file_path, newfile_path, write=False):
    print(f'- {newfile_path} does not exist')
    filename = os.path.basename(file_path)
    newfilename = os.path.basename(newfile_path)
    if write:
        print(f'- Will move {filename} to {newfilename}')  
        os.makedirs(os.path.dirname(newfile_path), exist_ok=True)                                    
        shutil.move(file_path, newfile_path)
        if os.path.isfile(newfile_path):
            file_path = newfile_path
        else:
            print(f'- shutil failed to move {filename} to {newfilename}')
    else:
        print(f'- Would move {filename} to {newfilename}')
    with open('movedfiles.log', "a") as f:
        f.write(f'{file_path} to {newfile_path}' + "\n")
                                

def mergefile(root, file_path, newfile_path, write=False, backup=False):
    # merge
    print(f'- Will try to merge {file_path} with {newfile_path}')
    try:
        oldst = read(newfile_path)
        partst = read(file_path)
        for tr in partst:
            try:
                oldst.append(tr)
                oldst.merge(fill_value=0)
            except:
                pass
        # save oldst & delete newst
        if len(oldst)==1 or 'soh' in root.lower():
            if check_write_read(oldst[0]):
                if write:
                    if backup:
                        shutil.copy2(newfile_path, os.path.join(BACKUPDIR, os.path.basename(newfile_path)))
                        shutil.copy2(file_path, os.path.join(BACKUPDIR, os.path.basename(file_path)))
                    os.makedirs(os.path.dirname(newfile_path), exist_ok=True)
                    oldst.write(newfile_path, format='MSEED')
                    os.remove(file_path)
                else:
                    print(f'- would write merged file to {newfile_path} and remove {file_path}')
                with open('mergedfiles.log', "a") as f:
                    f.write(f'{file_path} to {newfile_path} SUCCESS' + "\n")
            else:
                print('- merge/write/read failed')
                if write:
                    shutil.move(file_path, os.path.join(NOTMERGEDDIR, os.path.basename(file_path)))
                with open('mergedfiles.log', "a") as f:
                    f.write(f'{file_path} to {newfile_path} FAILED' + "\n")
        else:
            print(f'- got {len(oldst)} Trace objects from merging. should be 1')
            with open('mergedfiles.log', "a") as f:
                f.write(f'{file_path} to {newfile_path} but got {len(oldst)} traces' + "\n")
    except Exception as e:
        print(f'- merge failed for {newfile_path} and {file_path}?')
        print(e)    
        with open('mergedfiles.log', "a") as f:
            f.write(f'{file_path} to {newfile_path} CRASHED' + "\n")

def fix_sds_filenames(sds_directory, write=False, networks=None, backup=True):
    """Scan the SDS archive, rename files if necessary.
    Walks the directory tree in alphanumeric order, resuming from log if interrupted."""
    LOG_FILE = "fix_sds_filenames.log"
    processed_dirs = get_processed_dirs(LOG_FILE)
    for root, dirs, files in os.walk(sds_directory, topdown=True):
        # Sort directories and files alphanumerically
        dirs.sort()
        files.sort()

        # Skip directories that were already processed
        if root in processed_dirs:
            print(f"Skipping already processed: {root}")
            continue

        # Process files (this is where you would add your processing logic)
        print(f"Processing: {root}")
        for filename in files:
            if (networks and filename[0:2] in networks) or not networks:
                file_path = os.path.join(root, filename)

                # make sure filename is sensible first
                parts=filename.split('.')
                if len(parts)!=7:
                    print('\n', f'Bad filename {filename}')
                    if len(parts)>7:
                        newfilename = '.'.join(parts[0:7])
                        newfile_path = os.path.join(root, newfilename)   
                                             
                        if (('old' in parts[7] or 'seed' in parts[7] or 'ms' in parts[7]) or 'part' in parts[7]) : 
                            #if 'soh' in root.lower() and parts[7]='miniseed' and len(parts)==7: # we don't want to change SOH files .miniseed but .miniseed.part we do
                            #    continue
                            if os.path.isfile(newfile_path):
                                mergefile(root, file_path, newfile_path, write=write, backup=backup)
                            else:
                                movefile(file_path, newfile_path, write=write)

        # Mark this directory as processed
        save_processed_dir(root, LOG_FILE)


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


def fix_sds_ids(sds_directory, write=False, networks=None):
    """Scan the SDS archive, correct band codes, and rename files if necessary.
    Walks the directory tree in alphanumeric order, resuming from log if interrupted."""
    LOG_FILE = "fix_sds_ids.log"
    processed_dirs = get_processed_dirs(LOG_FILE)

    for root, dirs, files in os.walk(sds_directory, topdown=True):
        # Sort directories and files alphanumerically
        dirs.sort()
        files.sort()

        # Skip directories that were already processed
        if root in processed_dirs:
            print(f"Skipping already processed: {root}")
            continue

        # Process files (this is where you would add your processing logic)
        print(f"Processing: {root}")
        for filename in files:
            if (networks and filename[0:2] in networks) or not networks:
                file_path = os.path.join(root, filename)

                # make sure filename is sensible first
                parts=filename.split('.')
                if len(parts)!=7:
                    continue
 
                try:
                    # Read the MiniSEED file
                    stream = read(file_path)

                    for trace in stream:
                        if fix_trace_id(trace):

                            # Get the new SDS-compliant filename and path
                            new_file_path = trace2correct_sdsfullpath(sds_directory, trace)

                            if write:

                                # Ensure the directory exists
                                os.makedirs(os.path.dirname(new_file_path), exist_ok=True)

                                # Write the corrected MiniSEED file
                                stream.write(new_file_path, format="MSEED")

                                # Remove the old file
                                os.remove(file_path)

                                print(f"Renamed to: {new_file_path}")

                            else:
                                print(f'Would write {file_path} to {new_file_path}')

                except Exception as e:
                    print(f"Error processing {file_path}: {e}")

        # Mark this directory as processed
        save_processed_dir(root, LOG_FILE)


def find_channel_directory(path):
    """
    Extracts the channel directory (e.g., 'HHZ.D') from a file path.
    Assumes the directory structure is in SDS format.
    """
    parts = path.split(os.sep)
    for part in parts:
        if part.endswith(".D"):  # Looking for a directory like 'HHZ.D'
            return part
    return None  # Should not happen if the SDS structure is correct

def get_correct_directory_from_trace(base_dir, trace):
    """
    Determines the correct SDS channel directory for a trace.
    """
    '''
    year = str(trace.stats.starttime.year)
    network = trace.stats.network
    station = trace.stats.station
    location = trace.stats.location if trace.stats.location else "--"
    channel = trace.stats.channel
    dtype = "D"  # Always "D" for daily files

    correct_dir = os.path.join(base_dir, year, network, station, f"{channel}.{dtype}")
    return correct_dir
    '''
    return os.path.dirname(trace2correct_sdsfullpath(base_dir, trace))

def get_correct_directory_from_filename(base_dir, filename):
    """ Determines the correct SDS channel directory for a trace."""
    sds_tuple = parse_sds_filename(filename)
    if sds_tuple:
        network, station, location, channel, dtype, year, jday = sds_tuple

    correct_dir = os.path.join(base_dir, year, network, station, f"{channel}.{dtype}")
    return correct_dir

def get_correct_filename_from_trace(base_dir, trace):
    """ Determines the correct SDS file basename for a trace."""
    return os.path.basename(trace2correct_sdsfullpath(base_dir, trace))

def move_files_to_correct_directory(sds_directory, write=False):
    """ Walks through the SDS archive and moves MiniSEED files to the correct channel directory.
        Walks the directory tree in alphanumeric order, resuming from log if interrupted."""
    LOG_FILE = "move_files_to_correct_directory.log"
    processed_dirs = get_processed_dirs(LOG_FILE)
    for root, dirs, files in os.walk(sds_directory, topdown=True):
        # Sort directories and files alphanumerically
        dirs.sort()
        files.sort()

        # Skip directories that were already processed
        if root in processed_dirs:
            print(f"Skipping already processed: {root}")
            continue

        # Process files (this is where you would add your processing logic)
        print(f"Processing: {root}")
        current_channel_dir = find_channel_directory(root)
        if not current_channel_dir:
            continue  # Skip directories that aren't in *.D format

        for filename in files:

            file_path = os.path.join(root, filename)

            #net, sta, loc, chan, d, year, jday = filename.split('.')
            sds_tuple = parse_sds_filename(filename)
            if sds_tuple:
                network, station, location, channel, dtype, year, jday = sds_tuple     
                expected_dir = f'{channel}.{dtype}'
                if current_channel_dir != expected_dir:
                    # move
                    print(f"Moving: {file_path} -> {expected_dir}")

                    if write:

                        # Ensure the target directory exists
                        os.makedirs(expected_dir, exist_ok=True)

                        # Move the file
                        shutil.move(file_path, os.path.join(expected_dir, filename))
 
        # Mark this directory as processed
        save_processed_dir(root, LOG_FILE)                       


def remove_ds_store_files_and_empty_dirs(sds_directory):
    """ Recursively removes all .DS_Store files and deletes any empty directories. """
    for root, dirs, files in os.walk(sds_directory, topdown=False):  # Bottom-up to clean empty dirs
        # Remove .DS_Store files
        for file in files:
            if ".DS_Store" in file:
                file_path = os.path.join(root, file)
                try:
                    os.remove(file_path)
                    print(f"Deleted: {file_path}")
                except Exception as e:
                    print(f"Error deleting {file_path}: {e}")

        # Remove empty directories
        for d in dirs:
            dir_path = os.path.join(root, d)
            if not os.listdir(dir_path):  # Check if the directory is empty
                try:
                    os.rmdir(dir_path)
                    print(f"Removed empty directory: {dir_path}")
                except Exception as e:
                    print(f"Error removing {dir_path}: {e}")

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


def parse_sds_filename(filename):
    """
    Parses an SDS-style MiniSEED filename and extracts its components.
    Assumes filenames follow: NET.STA.LOC.CHAN.TYPE.YEAR.DAY.mseed
    """
    pattern = r"^(\w+)\.(\w+)\.(\w*)\.([A-Z]\w)\.([A-Z])\.(\d{4})\.(\d{3})$"
    pattern = r"^(\w*)\.(\w*)\.(\w*)\.(\w*)\.(\w*)\.(\d{4})\.(\d{3})$"
    match = re.match(pattern, filename)

    if match:
        network, station, location, channel, dtype, year, jday = match.groups()
        location = location if location else "--"  # Default for empty location codes
        if len(location)==1:
            location='0'+location
        return network, station, location, channel, dtype, year, jday
    else:
        return None

def move_files_to_sds_structure(source_dir, sds_base_dir, write=False, backup=False, fix_filename_but_do_not_move=False):
    """ Moves MiniSEED files from a flat directory into the correct SDS archive structure. 
        if fix_filename_but_do_not_move it just corrects the filename in the current directory
    """
    for filename in os.listdir(source_dir):
        # Parse the filename
        if parse_sds_filename(filename) or fix_filename_but_do_not_move:
            print(f'Processing {filename}')

            file_path = os.path.join(source_dir, filename)

            '''
            network, station, location, channel, dtype, year, jday = parse_sds_filename(filename):
            '''

            # read miniseed file
            try:
                this_st = read(file_path, format='MSEED')
            except Exception as e:
                print(e)
                print('Cannot read ',file_path)
                continue

            if len(this_st)==1:
                fix_trace_id(this_st[0])
                if fix_filename_but_do_not_move:
                    target_path = trace2correct_sdsfullpath(source_dir, this_st[0])
                else:
                    target_path = trace2correct_sdsfullpath(sds_base_dir, this_st[0])
                if file_path != target_path:
                    # check if target_path exists. if so, merge (and delete). otherwise move.
                    if os.path.isfile(target_path):
                        mergefile(sds_base_dir, file_path, target_path, write=write, backup=backup)
                    else:
                        movefile(file_path, target_path, write=write)

            else:
                print(f'got {len(this_st)} traces')
                continue
        else:
            print(f"Skipping file (does not match SDS format): {filename}")
            continue           



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
