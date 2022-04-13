import pickle
import os


def load_snapshots(snapshot_dir):
    snapshots = []
    for filename in os.listdir(snapshot_dir):
        if filename.endswith('.pkl') and "snapshot" in filename:
            with open(snapshot_dir+filename, 'rb') as inputfile:
                snapshots.append(pickle.load(inputfile))
    return snapshots


class Data:
    def __init__(self, arguments):
        self.snapshots = load_snapshots(arguments.input)
        print(f"Loaded snapshots of {len(self.snapshots)} timesteps.")

        self.consts = pickle.load(open(arguments.input+"constants.pkl", "rb"))
        print("Loaded constants:", list(self.consts.keys()))

        self.params = pickle.load(open(arguments.input+"parameters.pkl", "rb"))
        print("Loaded parameters:", list(vars(self.params).keys()))

        self.metrics = pickle.load(
            open(arguments.input+"cluster_metrics.pkl", "rb"))
        print("Loaded metrics:", list(self.metrics.keys()))

        if arguments.original_run is not None:
            self.og_metrics = pickle.load(open(arguments.original_run +
                                          "cluster_metrics.pkl", "rb"))
            print("Loaded original run metrics:", list(self.og_metrics.keys()))
        else:
            self.og_metrics = None

        for key in self.metrics:
            assert (len(self.metrics[key]) == len(self.snapshots)),\
                f"len({key})={len(key)}, there are {len(self.snapshots)} snaps"
