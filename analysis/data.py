import pickle
import os
import numpy


def load_snapshots(snapshot_dir):
    snapshots = []
    for filename in os.listdir(snapshot_dir):
        if filename.endswith('.pkl') and "snapshot" in filename:
            with open(snapshot_dir+filename, 'rb') as inputfile:
                snapshots.append(pickle.load(inputfile))
    return snapshots


def change_to_time_key(metrics):
    keys = metrics.keys()
    metrics_list = [dict(zip(keys, vals))
                    for vals in zip(*(metrics[k] for k in keys))]
    return {metrics['times'][i].number: metrics_list[i]
            for i in range(len(metrics_list))}


class Data:
    metrics_by_time = None
    dts = None

    def __init__(self, config):
        self.snapshots = load_snapshots(config.input)
        print(f"Loaded snapshots of {len(self.snapshots)} timesteps.")

        self.consts = pickle.load(open(config.input+"constants.pkl", "rb"))
        print("Loaded constants:", list(self.consts.keys()))

        self.params = pickle.load(open(config.input+"parameters.pkl", "rb"))
        print("Loaded parameters:", list(vars(self.params).keys()))

        self.metrics = pickle.load(
            open(config.input+"cluster_metrics.pkl", "rb"))
        print("Loaded metrics:", list(self.metrics.keys()))

        if config.input == ("/home/s1478621/job_outputs/"
                            "s65561_detailed_continue/"):
            for key in self.metrics:
                self.metrics[key] = self.metrics[key][9515:-3838]

            print(f"Metrics start at t={self.metrics['times'][0]}"
                  f" and end at t={self.metrics['times'][-1]}")

        if config.original_run is not None:
            self.og_metrics = pickle.load(open(config.original_run +
                                          "cluster_metrics.pkl", "rb"))
            print("Loaded original run metrics:", list(self.og_metrics.keys()))
        else:
            self.og_metrics = None

        for key in self.metrics:
            assert (len(self.metrics[key]) == len(self.snapshots)),\
                (f"len({key})={len(self.metrics[key])},"
                 f"there are {len(self.snapshots)} snaps")

    def by_time(self):
        if self.metrics_by_time is None:
            self.metrics_by_time = change_to_time_key(self.metrics)
            if self.og_metrics is not None:
                og_metrics_by_time = change_to_time_key(self.og_metrics)
                self.metrics_by_time[0.0] = og_metrics_by_time[0.0]
            if 0.0 not in self.metrics_by_time:
                raise IndexError(("t=0.0 not found, please provide original"
                                  "run data with -og"))
            return self.metrics_by_time
        else:
            return self.metrics_by_time

    def delta_ts(self):
        if self.dts is None:
            delta_ts = []
            for t_crc in self.metrics['t_crc']:
                delta_ts.append(0.01*t_crc.number)
            self.dts = numpy.array(delta_ts)
            return self.dts
        else:
            return self.dts
