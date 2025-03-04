# Define the file path where the input list is stored
import os
import sys
DEV_DIR = '/home/thompsong/Developer'
LIB_DIR = os.path.join(DEV_DIR, 'SoufriereHillsVolcano', 'lib')
sys.path.append(LIB_DIR)
os.listdir(LIB_DIR)
from libMVO import fix_trace_id_mvo

def fix_nslc_montserrat(traceid):
# function tr = fix_nslc_montserrat(tr)
    # FIX NET, STA, LOC, CHAN CODES ###
    # fix the network, channel and location
    oldnet, oldsta, oldloc, oldcha = traceid.split('.')

    network = 'MV'
    chan = oldcha
    sta = oldsta
    loc = oldloc
    if chan=='PRS' or chan[0:2]=='AP':
        chan='BDF'
    else:
        if chan[0]=='A':
           if loc == 'J':
               bandcode = 'S'
           else:
               bandcode = chan[0]
        else:
            try:
                if chan[1]=='B':
                    bandcode = 'B'
                else:
                    bandcode = chan[0]
            except:
                bandcode = chan[0]
        instrumentcode = 'H'
        if len(chan)<3:
            orientationcode = oldloc
        else:
            orientationcode = chan[2]

        chan = bandcode + instrumentcode + orientationcode

    if chan[0]=='A':
        print(network, sta, chan, loc)
        sys.exit('bugger!')
    return network + "." + sta + "." + loc + "." + chan


input_file_path = "SoufriereHillsVolcano/mappings.txt"  # Change this to your actual file name

# Read the file into a list
with open(input_file_path, "r") as file:
    input_list = [line.strip() for line in file]

# Display first few lines to verify
print("First few lines of the input list:", input_list[:5])


# Process the list: replace "->" with "," and add quotes
processed_list = [f'"{item.replace("->", ",")}"' for item in input_list]

# Save to file
file_path = "formatted_output.txt"
with open(file_path, "w") as file:
    file.writelines(f"{line.strip()}\n" for line in processed_list)

# Display the formatted file content as a DataFrame for verification
import pandas as pd
df = pd.DataFrame([line.replace('"', '').split(",") for line in processed_list], columns=["Old Trace ID", "New Trace ID"])
# Strip leading and trailing spaces from all string values in the DataFrame
df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)

# Get unique rows of a DataFrame
df = df.drop_duplicates()

#print(df)
from obspy import Trace
df['verified_id'] = None
for i,row in df.iterrows():
    #print('\n', i)
    tr = Trace()
    origid = row['Old Trace ID'].strip()
    tr.id = origid
    tr.stats.sampling_rate = 75.19
    fix_trace_id_mvo(tr)
    #altid = fix_nslc_montserrat(origid)
    #print(f'Propose to convert {origid} to {tr.id} or {altid}')
    #choice = input('Agree?')
    df.at[i,'verified_id'] = tr.id

print(df)

# Create a dictionary where verified_id is the key, and the value is a list of corresponding 'Old Trace ID'
mapping_dict = df.groupby('verified_id')['Old Trace ID'].apply(list).to_dict()

reverse_mapping = {}
forward_mapping = {}

# Display the resulting dictionary
for k,v in mapping_dict.items():
    net,sta,loc,chan = k.split('.')
    if len(v)>0:
        for count, this_v in enumerate(v):
            loc_count = count * 10
            if len(loc)>0 and loc.isnumeric():
                loc_count += int(loc)

            new_k = '.'.join([net, sta, str(loc_count).zfill(2), chan])
            reverse_mapping[new_k] = this_v
            forward_mapping[this_v] = new_k
    else:
        reverse_mapping[k] = v[0]
        forward_mapping[v[0]] = k

import json

# Define the file path to save the dictionary
forward_dict_file_path = "trace_ids_forward_mapping.json"

# Save the dictionary to a JSON file
with open(forward_dict_file_path, "w") as file:
    json.dump(forward_mapping, file, indent=4)

# Define the file path to save the dictionary
reverse_dict_file_path = "trace_ids_forward_mapping.json"

# Save the dictionary to a JSON file
with open(reverse_dict_file_path, "w") as file:
    json.dump(reverse_mapping, file, indent=4)

print(forward_mapping)