from amuse.units import nbody_system
from amuse.datamodel import Particles
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
    arguments = parser.parse_args()
    if not arguments.input.endswith("/"):
        arguments.input += "/"
    if not arguments.output.endswith("/"):
        arguments.output += "/"

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

    return snapshots, consts, params, metrics


def t_rh(N, r_h, G, M):
    return 0.138*N*r_h**(3/2)/(G**(1/2)*M**(1/2)*math.log(0.4*N))


def norm(vector):
    return numpy.sqrt(vector.dot(vector))


def f_ik_j(component_k, star_j):
    G = nbody_system.G
    m_j = star_j.mass
    m_ik = component_k.mass
    r_j = numpy.array([star_j.x.number, star_j.y.number, star_j.z.number])\
        | star_j.x.unit
    r_ik = numpy.array([component_k.x.number, component_k.y.number,
                        component_k.z.number]) | star_j.x.unit

    force_vector = -G*m_j*m_ik*(r_ik-r_j)/(norm(r_ik-r_j)**3)
    assert (force_vector.unit == nbody_system.mass * nbody_system.acceleration)

    return force_vector


def power_function(tuple, star):
    dEdt = 0 | nbody_system.energy / nbody_system.time
    for k in tuple:
        v_k = numpy.array([k.vx.number, k.vy.number, k.vz.number]) | k.vx.unit
        v_cm = tuple.center_of_mass_velocity()
        dEdt -= f_ik_j(k, star).dot(v_k - v_cm)

    assert (dEdt.unit == nbody_system.energy / nbody_system.time)

    return dEdt


def work_function(snapshots, tuple_ids, star_id, start_time, end_time):
    start_index = None
    for index, stars, time in zip(range(len(snapshots)), snapshots,
                                  metrics['times']):
        if time >= start_time and start_index is None:
            start_index = index

        if time >= end_time:
            end_index = index
            break

    power_functions = []
    for stars in snapshots[start_index:end_index]:
        tuple = Particles(0)
        for particle in stars:
            if particle.id in tuple_ids:
                tuple.add_particle(particle)
            if particle.id == star_id:
                star = particle

        power = power_function(tuple, star) * (1.0 | nbody_system.time)
        power_functions.append(power)

    work = numpy.cumsum(power_functions)
    assert (work[0].unit == nbody_system.energy),\
        f"[work]: {work[0].unit}, [energy]: {nbody_system.energy}"

    return work


def scatterplot(stars, time, arguments, first_binary=None, core_radius=None,
                density_centre=None):
    x_of_stars = stars.x.value_in(nbody_system.length)
    y_of_stars = stars.y.value_in(nbody_system.length)
    core_radius = core_radius.value_in(nbody_system.length)
    density_centre = density_centre.value_in(nbody_system.length)
    x_centre = density_centre[0]
    y_centre = density_centre[1]

    colours = []
    sizes = []
    for star in stars:
        if first_binary is not None:
            if star.id == first_binary[0].id or star.id == first_binary[1].id:
                colours.append('red')
                sizes.append(10)
            else:
                colours.append('blue')
                sizes.append(0.1)
        else:
            colours.append('blue')
            sizes.append(0.1)

    pyplot.scatter(x_of_stars, y_of_stars, s=sizes,
                   c=colours)
    pyplot.xlabel("x")
    pyplot.ylabel("y")
    if core_radius is not None:
        pyplot.xlim(x_centre-core_radius/2, x_centre+core_radius/2)
        pyplot.ylim(y_centre-core_radius/2, y_centre+core_radius/2)
    else:
        pyplot.xlim(-8, 8)
        pyplot.ylim(-8, 8)
    axes = pyplot.gca()
    axes.set_aspect('equal')
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
    for key in metrics:
        assert (len(metrics[key]) == len(snapshots)),\
            f"len({key})={len(key)}, there are {len(snapshots)} snaps."
    print('')

    t_max = metrics['times'][-1]
    t_rhi = t_rh(params.n, metrics['r50pc'][0], nbody_system.G,
                 1 | nbody_system.mass)
    times_crc = metrics['t_crc']

    if params.variable_delta:
        tau = numpy.cumsum(numpy.full(len(times_crc), 0.01))
    else:
        tau = numpy.cumsum(params.delta_t/times_crc)

    binaries_found = False
    for time, binaries, binding_E_kT in zip(metrics['times'],
                                            metrics['binaries'],
                                            metrics['binding_energies_kT']):
        if binaries != []:
            if not binaries_found:
                print((f"{len(binaries)} binaries found at {time}"))
                binaries_found = True
                first_binaries = binaries
                first_binding_energies = binding_E_kT
                t_bin = time

    if binaries_found:
        print("The first binaries are:")
        for binary in first_binaries:
            print(binary[0].id.number, ",", binary[1].id.number)
        print(f"Their energies are: {first_binding_energies}")
        print(f"They were found at t = {round(t_bin/t_rhi, 1)} t_rhi.")
    else:
        print("No binaries found.")

    print(f"t_max = {round(t_max/t_rhi, 1)} t_rhi")

    if arguments.scatter:
        for stars, time in zip(snapshots, metrics['times']):
            print(f"Plotting t={time} out of {t_max}", end="\r")
            if binaries_found:
                scatterplot(stars, time, arguments, first_binaries[0])
            else:
                scatterplot(stars, time, arguments, None)

    radiiplot(metrics, arguments)
