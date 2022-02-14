from amuse.units import nbody_system
from matplotlib import pyplot
import math
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


def t_rh(N, r_h, G, M):
    return 0.138*N*r_h**(3/2)/(G**(1/2)*M**(1/2)*math.log(0.4*N))


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
    print('')

    t_max = metrics['times'][-1]
    t_rhi = t_rh(params.n, metrics['r50pc'][0], nbody_system.G,
                 1 | nbody_system.mass)
    times_crc = metrics['t_crc']

    if params.variable_delta:
        tau = numpy.cumsum(numpy.full(len(times_crc), 0.01))
    else:
        tau = numpy.cumsum(params.delta_t/times_crc)

    print(f"t_crc: {times_crc}, delta_t: {params.delta_t}, tau: {tau}")

    binaries_found = False

    for time, binaries, binding_E_kT in zip(metrics['times'],
                                            metrics['binaries'],
                                            metrics['binding_Es_kT']):
        if binaries == []:
            print("No binaries at", time)
        else:
            print((f"Binaries with binding energy {binding_E_kT} "
                   f"found at {time}"))
            if not binaries_found:
                binaries_found = True
                first_binaries = binaries
                first_binding_energies = binding_E_kT
                t_bin = time

    if binaries_found:
        print("The first binaries are:")
        for binary in first_binaries:
            print(binary[0].id, ",", binary[1].id)
        print(f"Their energies are: {first_binding_energies}")
        print(f"They were found at t = {round(t_bin/t_rhi, 1)} t_rhi.")
    else:
        print("No binaries found.")

    print(f"t_max = {round(t_max/t_rhi, 1)} t_rhi")

    if arguments.scatter:
        for stars in snapshots:
            time = stars.get_timestamp().number
            print(f"Plotting t={time} out of {t_max}", end="\r")
            if binaries_found:
                scatterplot(stars, time, arguments.output, first_binaries[0])
            else:
                scatterplot(stars, time, arguments.output, None)

    radiiplot(metrics, arguments)
