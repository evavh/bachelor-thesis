import os
import pickle
import argparse

from amuse.units import nbody_system


def create_directory(directory_name):
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Data directory to be analysed",
                        default="analysis/data/")
    parser.add_argument("-o", "--output", help="Directory to put output",
                        default="analysis/plots/")
    parser.add_argument("--scatter", help="generate scatter plots.",
                        action='store_true')
    parser.add_argument("--scatter_core", help="generate scatter of core.",
                        action='store_true')
    parser.add_argument("-og", "--original_run", help="run from t0",
                        type=str)
    arguments = parser.parse_args()
    if not arguments.input.endswith("/"):
        arguments.input += "/"
    if not arguments.output.endswith("/"):
        arguments.output += "/"
    if arguments.original_run is not None \
            and not arguments.original_run.endswith("/"):
        arguments.original_run += "/"

    return arguments


def load_snapshots(snapshot_dir):
    snapshots = []
    for filename in os.listdir(snapshot_dir):
        if filename.endswith('.pkl') and "snapshot" in filename:
            with open(snapshot_dir+filename, 'rb') as inputfile:
                snapshots.append(pickle.load(inputfile))
    return snapshots


def load_data(arguments):
    snapshots = load_snapshots(arguments.input)
    print(f"Loaded snapshots of {len(snapshots)} timesteps.")

    consts = pickle.load(open(arguments.input+"constants.pkl", "rb"))
    print("Loaded constants:", list(consts.keys()))

    params = pickle.load(open(arguments.input+"parameters.pkl", "rb"))
    print("Loaded parameters:", list(vars(params).keys()))

    metrics = pickle.load(open(arguments.input+"cluster_metrics.pkl", "rb"))
    print("Loaded metrics:", list(metrics.keys()))

    if arguments.original_run is not None:
        og_metrics = pickle.load(open(arguments.original_run +
                                      "cluster_metrics.pkl", "rb"))
        print("Loaded original run metrics:", list(og_metrics.keys()))
    else:
        og_metrics = None

    for key in metrics:
        assert (len(metrics[key]) == len(snapshots)),\
            f"len({key})={len(key)}, there are {len(snapshots)} snaps."

    return snapshots, consts, params, metrics, og_metrics


def metrics_to_time_key(metrics):
    keys = metrics.keys()
    metrics_list = [dict(zip(keys, vals))
                    for vals in zip(*(metrics[k] for k in keys))]
    return {metrics['times'][i].value_in(nbody_system.time): metrics_list[i]
            for i in range(len(metrics_list))}
