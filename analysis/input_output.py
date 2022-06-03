import os
import argparse
import re

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
    parser.add_argument("--binary_ids", help="override binary ids.",
                        default=None, type=int, nargs=2)
    parser.add_argument("--interesting", help="mark stars with this id in plot.",
                        default=None, type=int, nargs=1)
    parser.add_argument("--energy_threshold", help="Eb threshold in kT.",
                        default=10.0, type=float)
    parser.add_argument("--fast_plot", help="skip loading snapshots.",
                        action='store_true')
    parser.add_argument("--scatter", help="generate scatter plots.",
                        action='store_true')
    parser.add_argument("--scatter_core", help="generate scatter of core.",
                        action='store_true')
    parser.add_argument("--load_work", help="load work data from pickle.",
                        action='store_true')
    parser.add_argument("--calc_work", help="calculate work (pricey).",
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


def extract_times(input_dir):
    files = os.listdir(input_dir)

    times = []
    for file in files:
        matches = re.match(r"snapshot_(.*) time.pkl", file)
        if matches is not None:
            times.append(float(matches[1]))

    return times | nbody_system.time
