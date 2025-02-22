import os
import json
import tkinter as tk
from tkinter import filedialog, messagebox
from obspy.clients.filesystem.sds import Client as SDSClient
from obspy import UTCDateTime, Stream
import numpy as np

CONFIG_FILE = "sds_config.json"

def load_last_directory():
    """Load the last used SDS archive directory from a config file."""
    if os.path.exists(CONFIG_FILE):
        try:
            with open(CONFIG_FILE, "r") as f:
                config = json.load(f)
            return config.get("last_sds_directory", "/")  # Fallback to root if missing
        except json.JSONDecodeError:
            return "/"
    return "/"

def save_last_directory(directory):
    """Save the last used SDS archive directory to a config file."""
    with open(CONFIG_FILE, "w") as f:
        json.dump({"last_sds_directory": directory}, f)

def browse_sds_directory():
    """Open a dialog to select the SDS archive directory, using the last used path."""
    initial_dir = load_last_directory()
    sds_dir = filedialog.askdirectory(title="Select SDS Archive Directory", initialdir=initial_dir)
    
    if sds_dir:  # If the user selected a directory, update the entry and save
        sds_entry.delete(0, tk.END)
        sds_entry.insert(0, sds_dir)
        save_last_directory(sds_dir)

def toggle_date_input():
    """Enables/disables input fields based on selected date mode."""
    if date_mode.get() == "ymd":
        year_entry.config(state="normal")
        month_entry.config(state="normal")
        day_entry.config(state="normal")
        julian_entry.config(state="disabled")
    else:
        year_entry.config(state="normal")
        month_entry.config(state="disabled")
        day_entry.config(state="disabled")
        julian_entry.config(state="normal")



def count_nonzero_samples(trace):
    """Returns the count of non-zero samples in a trace."""
    return np.count_nonzero(trace.data)

def plot_seismic_data():
    """Loads, downsamples, and plots the seismic and infrasound channels with the highest number of non-zero samples."""
    sds_directory = sds_entry.get().strip()
    network = network_entry.get().strip()
    year_str = year_entry.get().strip()
    print('\nProcessing inputs')

    if date_mode.get() == "ymd":
        month_str = month_entry.get().strip()
        day_str = day_entry.get().strip()
    else:
        julian_str = julian_entry.get().strip()

    if not sds_directory or not network or not year_str:
        messagebox.showerror("Input Error", "Please fill in all fields.")
        return

    try:
        year = int(year_str)
        
        if date_mode.get() == "ymd":
            month = int(month_str)
            day = int(day_str)
            date = UTCDateTime(year, month, day)
        else:
            julian_day = int(julian_str)
            if julian_day < 1 or julian_day > 366:
                messagebox.showerror("Input Error", "Julian day must be between 1 and 366.")
                return
            date = UTCDateTime(f"{year}-{julian_day:03d}")

        print(f'Reading waveform data for {date} and network={network} from {sds_directory}')
        client = SDSClient(sds_directory)

        # Load all available data for the given network and date
        stream = client.get_waveforms(network=network, station="*", location="*", channel="*", starttime=date, endtime=date + 86400)

        if len(stream) == 0:
            messagebox.showwarning("No Data", f"No waveform data found for {date.date}.")
            return
       
        print(stream.__str__(extended=True))

        # Downsample to 1 Hz if necessary
        target_sampling_rate = 1.0
        print(f'Downsampling traces', end=' ')
        for trace in stream:
            print(trace.id, end=' ')
            if trace.stats.sampling_rate > target_sampling_rate:
                decimation_factor = int(trace.stats.sampling_rate / target_sampling_rate)
                if decimation_factor > 1:
                    trace.decimate(decimation_factor, no_filter=True)
        print('\n')

        # Group traces by station-location and select the ones with the most non-zero samples
        station_location_data = {}
        print('Subsetting channels to plot')
        for trace in stream:
            key = (trace.stats.station, trace.stats.location)
            channel = trace.stats.channel
            nonzero_count = count_nonzero_samples(trace)

            if key not in station_location_data:
                station_location_data[key] = {"seismic": None, "infrasound": None}

            # Seismic selection: Highest non-zero sample count among ?HZ, ?HN, ?HE
            if channel[1]=='H':
                if (
                    station_location_data[key]["seismic"] is None or
                    nonzero_count > count_nonzero_samples(station_location_data[key]["seismic"])
                ):
                    station_location_data[key]["seismic"] = trace  # Pick the trace with the most non-zero samples

            # Infrasound selection: Highest non-zero sample count among ?DF, ?DG, ?D
            if channel[1]=='D':
                if (
                    station_location_data[key]["infrasound"] is None or
                    nonzero_count > count_nonzero_samples(station_location_data[key]["infrasound"])
                ):
                    station_location_data[key]["infrasound"] = trace  # Pick the trace with the most non-zero samples

        # Create a single stream for all selected traces
        selected_stream = Stream()
        for channels in station_location_data.values():
            if channels["seismic"]:
                selected_stream += channels["seismic"]
            if channels["infrasound"]:
                selected_stream += channels["infrasound"]

        # Plot all selected traces in one figure
        print('Plotting')
        if len(selected_stream) > 0:
            print(f"Plotting {len(selected_stream)} traces in one figure.")
            selected_stream.plot(equal_scale=False)

    except Exception as e:
        messagebox.showerror("Error", f"An error occurred: {e}")



# GUI Setup
root = tk.Tk()
root.title("SDS Archive Data Loader")

# SDS Directory Selection
tk.Label(root, text="SDS Archive Path:").grid(row=0, column=0, padx=5, pady=5, sticky="w")
sds_entry = tk.Entry(root, width=50)
sds_entry.grid(row=0, column=1, padx=5, pady=5)
tk.Button(root, text="Browse", command=lambda: sds_entry.insert(0, filedialog.askdirectory())).grid(row=0, column=2, padx=5, pady=5)

# Network Name
tk.Label(root, text="Network Name:").grid(row=1, column=0, padx=5, pady=5, sticky="w")
network_entry = tk.Entry(root, width=20)
network_entry.grid(row=1, column=1, padx=5, pady=5)

# Date Mode Selection (Radio Buttons)
date_mode = tk.StringVar(value="ymd")  # Default to Year-Month-Day
tk.Label(root, text="Date Entry Mode:").grid(row=2, column=0, padx=5, pady=5, sticky="w")
tk.Radiobutton(root, text="Year-Month-Day", variable=date_mode, value="ymd", command=toggle_date_input).grid(row=2, column=1, sticky="w")
tk.Radiobutton(root, text="Year-Julian Day", variable=date_mode, value="jday", command=toggle_date_input).grid(row=2, column=2, sticky="w")

# Year Input
tk.Label(root, text="Year:").grid(row=3, column=0, padx=5, pady=5, sticky="w")
year_entry = tk.Entry(root, width=10)
year_entry.grid(row=3, column=1, padx=5, pady=5)

# Month & Day Inputs
tk.Label(root, text="Month:").grid(row=4, column=0, padx=5, pady=5, sticky="w")
month_entry = tk.Entry(root, width=5)
month_entry.grid(row=4, column=1, padx=5, pady=5)

tk.Label(root, text="Day:").grid(row=4, column=2, padx=5, pady=5, sticky="w")
day_entry = tk.Entry(root, width=5)
day_entry.grid(row=4, column=3, padx=5, pady=5)

# Julian Day Input (Initially Disabled)
tk.Label(root, text="Julian Day (1-366):").grid(row=5, column=0, padx=5, pady=5, sticky="w")
julian_entry = tk.Entry(root, width=10, state="disabled")
julian_entry.grid(row=5, column=1, padx=5, pady=5)

# Load & Plot Button
plot_button = tk.Button(root, text="Load & Plot Data", command=plot_seismic_data)
plot_button.grid(row=6, column=0, columnspan=3, pady=10)

# Initialize Date Input Mode
toggle_date_input()

root.mainloop()
