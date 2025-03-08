import numpy as np
import tkinter as tk
from tkinter import messagebox
from obspy import read, UTCDateTime, Stream
from obspy.signal.trigger import coincidence_trigger
import random
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

### --- Step 1: GUI for Selecting Traces --- ###
class TraceSelectorApp:
    def __init__(self, master, stream):
        self.master = master
        self.master.title("Select Traces")

        self.stream = stream
        self.selected_traces = []

        self.fig, self.axs = plt.subplots(len(stream), 1, figsize=(8, len(stream) * 1.5), sharex=True)
        plt.subplots_adjust(hspace=0.5)

        self.checkboxes = []
        self.vars = []

        # Create checkboxes and plot traces
        for i, tr in enumerate(stream):
            var = tk.BooleanVar(value=False)
            chk = tk.Checkbutton(master, text=tr.id, variable=var)
            chk.pack(anchor="w")
            self.checkboxes.append(chk)
            self.vars.append(var)

            # Plot the traces
            data_norm = tr.data / max(abs(tr.data)) if max(abs(tr.data)) != 0 else tr.data
            self.axs[i].plot(tr.times("matplotlib"), data_norm, label=tr.id)
            self.axs[i].legend(loc="upper right")

        self.canvas = FigureCanvasTkAgg(self.fig, master)
        self.canvas.get_tk_widget().pack()

        # Selection button
        self.select_button = tk.Button(master, text="Select Traces", command=self.select_traces)
        self.select_button.pack()

    def select_traces(self):
        self.selected_traces = [self.stream[i] for i, var in enumerate(self.vars) if var.get()]
        if len(self.selected_traces) == 0:
            messagebox.showerror("Error", "Please select at least one trace!")
        else:
            self.master.quit()  # Close GUI when selection is done

### --- Step 2: Event Picking (Start & End) --- ###
def pick_event_times(selected_stream):
    """Pick event start and end times in the **same figure** as the traces, with correctly scaled x-axis."""
    
    fig, ax = plt.subplots(figsize=(10, 4))

    # Get start and end time of the stream
    start_time = selected_stream[0].stats.starttime.matplotlib_date
    end_time = selected_stream[0].stats.endtime.matplotlib_date

    # Plot all selected traces
    for tr in selected_stream:
        times = tr.times("matplotlib") + start_time  # Convert relative times to absolute matplotlib times
        data_norm = tr.data / max(abs(tr.data)) if max(abs(tr.data)) != 0 else tr.data
        ax.plot(times, data_norm, label=tr.id)

    ax.set_xlim(start_time, end_time)  # Ensure x-axis matches actual time range
    ax.set_title("🖱 Click to select event start time")
    ax.legend(loc="upper right")

    cursor_line = ax.axvline(x=start_time, color='r', linestyle='dotted', lw=1)  # Dotted vertical cursor
    picked_lines = []

    def on_mouse_move(event):
        """Move dotted vertical cursor line."""
        if event.xdata:
            cursor_line.set_xdata(event.xdata)
            fig.canvas.draw_idle()

    def on_mouse_click(event):
        """Handle mouse clicks to pick event start and end times."""
        if len(picked_lines) < 2 and event.xdata:
            ax.axvline(x=event.xdata, color='r', linestyle='solid', lw=2)  # Draw permanent vertical line
            picked_lines.append(UTCDateTime(event.xdata))  # Convert to UTCDateTime

            if len(picked_lines) == 1:
                ax.set_title("🖱 Click to select event end time")
            elif len(picked_lines) == 2:
                plt.close()

            fig.canvas.draw_idle()

    fig.canvas.mpl_connect('motion_notify_event', on_mouse_move)
    fig.canvas.mpl_connect('button_press_event', on_mouse_click)
    
    plt.show()

    if len(picked_lines) != 2:
        print("⚠️ Error: You must pick exactly two points!")
        return None

    event_start, event_end = picked_lines
    print(f"✅ Picked event start: {event_start}, end: {event_end}")
    return [event_start, event_end]

### --- Step 3: Monte Carlo Simulation for Best Parameters --- ###
def run_monte_carlo(stream, event_start, event_end, n_trials=100):
    """Runs Monte Carlo simulations with different trigger parameters."""
    algorithms = ["recstalta", "classicstalta", "zdetect", "carlstatrig"]
    best_params = None
    best_misfit = float("inf")
    min_triggers = 4

    for _ in range(n_trials):
        algorithm = random.choice(algorithms)
        sta = random.uniform(1, 5)
        lta = random.uniform(sta * 4, sta * 25)
        thr_on = random.uniform(1.5, 5)
        thr_off = random.uniform(thr_on / 10, thr_on / 2)

        print(f"Trying {algorithm} with STA={sta}, LTA={lta}, THR_ON={thr_on}, THR_OFF={thr_off}")

        # Run coincidence_trigger
        triggers = coincidence_trigger(algorithm, thr_on, thr_off, stream, min_triggers, sta=sta, lta=lta)

        detected_start_times = [UTCDateTime(trig["time"]) for trig in triggers]
        detected_end_times = [UTCDateTime(trig["time"]) + trig["duration"] for trig in triggers]

        if detected_start_times and detected_end_times:
            detected_start = min(detected_start_times)
            detected_end = max(detected_end_times)
            misfit = abs(detected_start - event_start) + abs(detected_end - event_end)

            print(f"Detected event: {detected_start} to {detected_end}, Misfit: {misfit}")

            if misfit < best_misfit:
                best_misfit = misfit
                best_params = {
                    "algorithm": algorithm,
                    "sta": sta,
                    "lta": lta,
                    "thr_on": thr_on,
                    "thr_off": thr_off,
                    "misfit": misfit,
                }

    return best_params

### --- Step 4: Main Workflow --- ###
def run_event_detection(stream, n_trials=50):
    """
    Full workflow for trace selection, event picking, and parameter tuning.
    """

    # Step 1: Select Traces using Tkinter GUI
    root = tk.Tk()
    app = TraceSelectorApp(root, stream)
    root.mainloop()
    selected_stream = Stream(app.selected_traces)

    if len(selected_stream) == 0:
        print("⚠️ Error: No traces selected. Exiting.")
        return None, None, None

    # Step 2: Pick Event Start & End
    picked_times = pick_event_times(selected_stream)
    if picked_times is None:
        print("⚠️ Error: You must pick exactly two points. Exiting.")
        return None, None, None
    event_start, event_end = picked_times

    # Step 3: Optimize Parameters
    best_params = run_monte_carlo(selected_stream, event_start, event_end, n_trials)

    # Step 4: Display Best Parameters
    if best_params:
        print("\n✅ Best Parameters Found:")
        print(best_params)
        return selected_stream, picked_times, best_params
    else:
        print("⚠️ No optimal parameters found.")
        return selected_stream, picked_times, None  

### --- Step 5: Run the Standalone Application --- ###
if __name__ == "__main__":
    # Load waveform data (modify with actual file)
    stream = read('test.mseed', format='MSEED')  # Modify to load actual data
    selected_stream, picked_times, best_params = run_event_detection(stream, n_trials=50)