import os
import pickle
import pandas
import argparse

from amuse.units import nbody_system


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--delta_t", help="time between snapshots",
                        default=1.0, type=float)
    parser.add_argument("-v", "--variable_delta", help="delta_t = 0.01t_crc",
                        action='store_true')
    parser.add_argument("-n", help="number of stars to simulate",
                        default=100, type=int)
    parser.add_argument("-s", "--random_seed", help="random number seed",
                        default=None, type=int)
    parser.add_argument("-e", "--bs_tolerance", help="tolerance for Brutus",
                        default=None, type=float)
    parser.add_argument("-t", "--t_end", help="time to force end the sim",
                        default=None, type=float)
    parser.add_argument("-T", "--start_time", help="snapshot time to load",
                        default=None, type=float)
    parser.add_argument("-o", "--output_folder", help="where to put output",
                        default="simulation/output/", type=str)
    parser.add_argument("-i", "--snapshot_input", help="snaps to start from",
                        default=None, type=str)
    parser.add_argument("-b", "--minimum_Eb_kT", help="minimum binding E / kT",
                        default=10, type=float)
    parser.add_argument("-r", "--reverse", help="run the sim in reverse",
                        action='store_true')

    params = parser.parse_args()

    if params.reverse:
        print("Running a reversed simulation.")
        assert params.bs_tolerance is not None, "Please specify bs tolerance"
        assert params.t_end is not None, \
            "When reversed there must be an end time."
        assert params.start_time is not None, \
            "When reversed there must be a start time."
        assert params.start_time > params.t_end, \
            "When reversed start time must be later than end time."

        params.t_end = params.start_time - params.t_end \
            + params.start_time
        print(f"Set end time to {params.t_end} for the integrator's sake.")
    elif params.start_time is not None and params.t_end is not None:
        assert params.start_time < params.t_end, \
            "When not reversed start time must be before end time."

    params.delta_t = params.delta_t | nbody_system.time
    if params.t_end is not None:
        params.t_end = params.t_end | nbody_system.time
    if params.start_time is not None:
        assert params.snapshot_input is not None, \
            "When there is a start time there must be snapshots."
        params.start_time = params.start_time | nbody_system.time

    if not params.output_folder.endswith('/'):
        params.output_folder += '/'

    return params


def create_directory(directory_name):
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)


def remove_file(filename):
    if os.path.exists(filename):
        os.remove(filename)


def pickle_object(object, filename, params):
    object_path = params.output_folder+filename
    pickle.dump(object, open(object_path, "wb"))


def unpickle_object(filename, params):
    object_path = params.snapshot_input+'/'+filename
    return pickle.load(open(object_path, 'rb'))


def round_csv(filename, decimals):
    data = pandas.read_csv(filename)
    data = data.iloc[3:]
    data = data.astype(float).round(3)
    data.to_csv(filename[:-4]+"_rounded.csv", index=False)
