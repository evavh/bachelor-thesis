from amuse.units import nbody_system
from matplotlib import pyplot
import numpy
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

    return snapshots, consts, params, metrics


def scatterplot(stars, time, arguments, first_binary):
    x_of_stars = stars.x.value_in(nbody_system.length)
    y_of_stars = stars.y.value_in(nbody_system.length)

    colours = []
    sizes = []
    for star in stars:
        if first_binary is not None:
            if star.id == first_binary[0].id or star.id == first_binary[1].id:
                colours.append('red')
                sizes.append(10)
            else:
                colours.append('blue')
                sizes.append(0.5)
        else:
            colours.append('blue')
            sizes.append(0.5)

    pyplot.scatter(x_of_stars, y_of_stars, s=sizes,
                   c=colours)
    pyplot.xlabel("x")
    pyplot.ylabel("y")
    pyplot.savefig(arguments.output+"scatter/"+str(time)+".svg", format='svg')
    pyplot.clf()


def radiiplot(metrics, arguments):
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
    pyplot.savefig(arguments.output+"radii.png")


if __name__ == '__main__':
    arguments = parse_arguments()

    create_directory(arguments.output)
    create_directory(arguments.output+"scatter")
    create_directory(arguments.input)

    snapshots, consts, params, metrics = load_data(arguments)

    t_max = snapshots[0].get_timestamp()

    times = metrics['times']
    binaries = metrics['binaries']
    binding_Es_kT = metrics['binding_energies_kT']

    binaries_found = False

    for time, binaries, binding_E_kT in zip(times, binaries, binding_Es_kT):
        if binaries == []:
            print("No binaries at", time)
        else:
            print((f"Binaries {binaries} with binding energy {binding_E_kT} "
                   f"found at {time}"))
            if not binaries_found:
                binaries_found = True
                first_binaries = binaries
                first_binding_energies = binding_E_kT
                first_binaries_time = time

    if binaries_found:
        print("The first binaries are:")
        for binary in first_binaries:
            print(binary[0].id, ",", binary[1].id)
        print(f"Their energies are: {first_binding_energies}")
        print(f"They were found at t={first_binaries_time}.")
    else:
        print("No binaries found.")

    if arguments.scatter:
        for stars in snapshots:
            time = stars.get_timestamp().number
            print(f"Plotting t={time} out of {t_max}", end="\r")
            if binaries_found:
                scatterplot(stars, time, arguments.output, first_binaries[0])
            else:
                scatterplot(stars, time, arguments.output, None)

    radiiplot(metrics, arguments)
