import numpy as np
import tkinter as tk
from tkinter import messagebox
from obspy import read, UTCDateTime, Stream
from obspy.signal.trigger import coincidence_trigger
import random
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from pprint import pprint

class SeismicGUI:
    def __init__(self, master, stream):
        self.master = master
        self.master.title("Seismic Event Detection")

        self.stream = stream
        self.selected_traces = stream.copy()  # Default: All traces selected
        self.picked_times = []
        self.mode = "select_traces"

        self.frame_left = tk.Frame(master)
        self.frame_left.pack(side=tk.LEFT, fill=tk.Y, padx=10, pady=10)

        self.vars = []
        self.checkboxes = []
        for i, tr in enumerate(stream):
            var = tk.BooleanVar(value=True)  # Default to checked
            chk = tk.Checkbutton(self.frame_left, text=tr.id, variable=var)
            chk.pack(anchor="w")
            self.checkboxes.append(chk)
            self.vars.append(var)

        self.fig, self.axs = plt.subplots(len(stream), 1, figsize=(10, 6), sharex=True)
        plt.subplots_adjust(hspace=0.3)

        self.canvas = FigureCanvasTkAgg(self.fig, master)
        self.canvas.get_tk_widget().pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        self.plot_traces()

        self.button = tk.Button(self.frame_left, text="Select Traces", command=self.select_traces)
        self.button.pack(pady=10)

        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.canvas.mpl_connect('button_press_event', self.on_mouse_click)

    def plot_traces(self):
        """Plots traces separately in linked subplots, ensuring correct time scaling and linking."""
        for ax in self.axs:
            ax.clear()

        start_time = min(tr.stats.starttime for tr in self.stream)
        end_time = max(tr.stats.endtime for tr in self.stream)

        start_time_num = mdates.date2num(start_time.datetime)
        end_time_num = mdates.date2num(end_time.datetime)

        if self.mode == "select_traces":
            self.fig.suptitle("Select traces using checkboxes, then click 'Select Traces'")

            for i, tr in enumerate(self.stream):
                absolute_times = mdates.date2num([start_time + t for t in tr.times()])
                data_norm = tr.data / max(abs(tr.data)) if max(abs(tr.data)) != 0 else tr.data
                self.axs[i].plot(absolute_times, data_norm, label=tr.id)
                self.axs[i].legend(loc="upper right")
                self.axs[i].set_ylabel(tr.id)

            self.axs[-1].set_xlabel("Time (UTC)")

        elif self.mode == "pick_times":
            self.fig.suptitle("Click to select event start and end times")

            for i, tr in enumerate(self.selected_traces):
                absolute_times = mdates.date2num([tr.stats.starttime + t for t in tr.times()])
                data_norm = tr.data / max(abs(tr.data)) if max(abs(tr.data)) != 0 else tr.data
                self.axs[i].plot(absolute_times, data_norm, label=tr.id)
                self.axs[i].legend(loc="upper right")
                self.axs[i].set_ylabel(tr.id)

            self.axs[-1].set_xlabel("Time (UTC)")

            self.axs[0].set_xlim(start_time_num, end_time_num)

            self.cursor_lines = [ax.axvline(x=start_time_num, color='r', linestyle='dotted', lw=1) for ax in self.axs]

        for ax in self.axs:
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d %H:%M:%S"))

        self.canvas.draw()

    def select_traces(self):
        """Handles trace selection and switches to event picking mode."""
        self.selected_traces = [self.stream[i] for i, var in enumerate(self.vars) if var.get()]

        if not self.selected_traces:
            messagebox.showerror("Error", "Please select at least one trace!")
            return

        self.mode = "pick_times"
        self.button.config(text="Pick Start and End Times", state=tk.DISABLED)
        self.plot_traces()

    def on_mouse_move(self, event):
        """Move dotted cursor in pick_times mode across all subplots."""
        if self.mode == "pick_times" and event.xdata:
            for cursor_line in self.cursor_lines:
                cursor_line.set_xdata(event.xdata)
            self.canvas.draw()

    def on_mouse_click(self, event):
        """Handles event start and end time selection."""
        if self.mode == "pick_times" and event.xdata:
            picked_time = UTCDateTime(mdates.num2date(event.xdata))
            for ax in self.axs:
                ax.axvline(x=event.xdata, color='r', linestyle='solid', lw=2)

            self.picked_times.append(picked_time)
            print(f"✅ Picked event time: {picked_time}")

            if len(self.picked_times) == 1:
                self.fig.suptitle("Click to select event end time")
            elif len(self.picked_times) == 2:
                plt.close(self.fig)  # ✅ Closes only the current figure
                plt.close('all')  # ✅ Ensures all figures are closed                
                self.master.quit()

            self.canvas.draw()

def run_monte_carlo(stream, event_start, event_end, n_trials=100):
    """Runs Monte Carlo simulations with different trigger parameters."""
    algorithms = ["recstalta", "classicstalta", "zdetect", "carlstatrig", "delayedstalta"]
    # classicstalta: Classic STA/LTA algorithm
        # classic_sta_lta(tr.data, nsta, nlta) nsta is number of samples
    # recstalta: Recursive STA/LTA algorithm
        # recursive_sta_lta(tr.data, nsta, nlta)
    # zdetect: Z-detect algorithm
        # z_detect(tr.data, nsta) https://docs.obspy.org/packages/autogen/obspy.signal.trigger.z_detect.html#obspy.signal.trigger.z_detect
    # carlstatrig: Carl-Sta-Trig algorithm
        # carl_sta_trig(tr.data, nsta, nlta, ratio, quiet) https://docs.obspy.org/packages/autogen/obspy.signal.trigger.carl_sta_trig.html#obspy.signal.trigger.carl_sta_trig
    # delayed_sta_lta: Delayed STA/LTA algorithm
        # delayed_sta_lta(tr.data, nsta, nlta)
    best_params = None
    best_misfit = float("inf")
    min_triggers = 4

    for _ in range(n_trials):
        algorithm = random.choice(algorithms)
        sta = random.uniform(1, 5)
        lta = random.uniform(sta * 4, sta * 25)
        thr_on = random.uniform(1.5, 5)
        thr_off = random.uniform(thr_on / 10, thr_on / 2)
        ratio = random.uniform(1, 2)

        print(f"Trying {algorithm} with STA={sta}, LTA={lta}, THR_ON={thr_on}, THR_OFF={thr_off}")

        # Run coincidence_trigger
        # coincidence_trigger(trigger_type, thr_on, thr_off, stream, thr_coincidence_sum, \
        # trace_ids=None, max_trigger_length=1000000.0, delete_long_trigger=False, trigger_off_extension=0, details=False, event_templates={}, similarity_threshold=0.7, **options)
        if algorithm == "zdetect":
            triggers = coincidence_trigger(algorithm, thr_on, thr_off, stream, min_triggers, sta=sta)
        elif algorithm == "carlstatrig":
            triggers = coincidence_trigger(algorithm, thr_on, thr_off, stream, min_triggers, sta=sta, lta=lta, ratio=1, quiet=True)
        else:
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

def run_event_detection(stream, n_trials=50):
    """Runs GUI and returns selected traces, event times, and best parameters."""
    root = tk.Tk()
    app = SeismicGUI(root, stream)
    root.mainloop()

    if len(app.picked_times) != 2:
        print("⚠️ Error: No valid event times selected.")
        return None, None, None, None

    event_start, event_end = app.picked_times

    best_params = run_monte_carlo(app.selected_traces, event_start, event_end, n_trials)
    out_stream = Stream(traces=app.selected_traces)

    return out_stream, [event_start, event_end], best_params

if __name__ == "__main__":
    stream = read('test.mseed', format='MSEED')  # Modify to load actual data
    selected_stream, picked_times, best_params = run_event_detection(stream, n_trials=50)

    print(f"\n Selected {len(selected_stream)} traces")
    print(f" Event Start: {picked_times[0]}, Event End: {picked_times[1]}")
    print(f" Best Parameters: {best_params}")

