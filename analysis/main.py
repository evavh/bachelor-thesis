from amuse.io import read_set_from_file
from amuse.units import nbody_system
from matplotlib import pyplot
import pickle
import os
import argparse


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
    return parser.parse_args()


def scatterplot(stars, time, output_dir):
    pyplot.scatter(stars.x.value_in(nbody_system.length),
                   stars.y.value_in(nbody_system.length),
                   s=10*stars.mass/stars.mass.max())
    pyplot.xlabel("x")
    pyplot.ylabel("y")
    pyplot.savefig(output_dir+"scatter/"+str(time)+".png")
    pyplot.clf()


def radiiplot(metrics, output_dir):
    pyplot.plot(metrics["times"].value_in(nbody_system.time),
                metrics["rcore"].value_in(nbody_system.length), c='b')
    pyplot.plot(metrics["times"].value_in(nbody_system.time),
                metrics["rvir"].value_in(nbody_system.length), c='r')
    pyplot.plot(metrics["times"].value_in(nbody_system.time),
                metrics["r10pc"].value_in(nbody_system.length), c='k', ls="--",
                lw=1)
    pyplot.plot(metrics["times"].value_in(nbody_system.time),
                metrics["r50pc"].value_in(nbody_system.length), c='k', ls="--",
                lw=3)
    pyplot.xlabel("time")
    pyplot.ylabel("radius")
    pyplot.semilogy()
    pyplot.savefig(output_dir+"radii.png")


if __name__ == '__main__':
    arguments = parse_arguments()
    input_dir = arguments.input
    output_dir = arguments.output
    scatter = arguments.scatter

    create_directory(output_dir)
    create_directory(output_dir+"scatter")
    create_directory(input_dir)

    snapshots = read_set_from_file(input_dir+"snapshots.hdf5", "hdf5")
    tmax = len(list(snapshots.history))
    print(f"Loaded snapshots of {tmax} timesteps.")
    metrics = pickle.load(open(input_dir+"cluster_metrics.pkl", "rb"))
    print("Loaded metrics:", list(metrics.keys()))

    print(f"The first binaries are: {metrics['first_binaries']}.")
    print(f"Their energies are: {metrics['first_binary_energies_kT']}.")
    print(f"They were found at t={metrics['first_binary_time']}.")

    if scatter:
        for stars in snapshots.history:
            time = stars.get_timestamp().number
            print(f"Plotting t={time} out of {tmax}", end="\r")
            scatterplot(stars, time, output_dir)

    radiiplot(metrics, output_dir)
