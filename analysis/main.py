from amuse.units import nbody_system
from amuse.datamodel import Particles
from matplotlib import pyplot

import math
import numpy

import input_output


def find_first_binary(metrics, t_rhi):
    binaries_found = False
    first_binary = None
    t_bin = None

    for time, binaries in zip(metrics['times'], metrics['binaries']):
        if binaries != []:
            if not binaries_found:
                binaries_found = True
                first_binary = binaries[0]
                first_binary_ids = (first_binary[0].id.number,
                                    first_binary[1].id.number)
                t_bin = time

                print(f"The first binary is {first_binary_ids}")
                print((f"It was found at t = {t_bin.number} = "
                       f"{round(t_bin/t_rhi, 1)} t_rhi."))
                if len(binaries) > 1:
                    print((f"NOTE: {len(binaries)-1} more were found at this"
                           " time!"))

    if not binaries_found:
        print("No binaries found.")

    return first_binary, t_bin


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


def scatterplot(stars, time, output_folder, first_binary=None,
                metrics_by_time=None):
    if not output_folder.endswith("/"):
        output_folder += "/"

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
                sizes.append(0.1)
        else:
            colours.append('blue')
            sizes.append(0.1)

    pyplot.scatter(x_of_stars, y_of_stars, s=sizes, c=colours)
    pyplot.xlabel("x")
    pyplot.ylabel("y")
    if metrics_by_time is not None:
        time_number = time.value_in(nbody_system.time)
        core_radius = metrics_by_time[time_number]['rcore']
        core_radius = core_radius.value_in(nbody_system.length)
        density_centre = metrics_by_time[time_number]['density_centre']
        density_centre = density_centre.value_in(nbody_system.length)
        x_centre = density_centre[0]
        y_centre = density_centre[1]

        pyplot.xlim(x_centre-core_radius/2, x_centre+core_radius/2)
        pyplot.ylim(y_centre-core_radius/2, y_centre+core_radius/2)
    else:
        pyplot.xlim(-8, 8)
        pyplot.ylim(-8, 8)
    axes = pyplot.gca()
    axes.set_aspect('equal')
    pyplot.savefig(output_folder+str(time)+".svg", format='svg')
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
    arguments = input_output.parse_arguments()

    input_output.create_directory(arguments.output)
    input_output.create_directory(arguments.input)

    snapshots, consts, params, metrics = input_output.load_data(arguments)
    print('')

    metrics_by_time = input_output.metrics_to_time_key(metrics)

    t_max = metrics['times'][-1]
    t_rhi = t_rh(params.n, metrics['r50pc'][0], nbody_system.G,
                 1 | nbody_system.mass)
    times_crc = metrics['t_crc']

    if params.variable_delta:
        tau = numpy.cumsum(numpy.full(len(times_crc), 0.01))
    else:
        tau = numpy.cumsum(params.delta_t/times_crc)

    first_binary, t_bin = find_first_binary(metrics, t_rhi)
    binaries_found = (first_binary is not None)

    print(f"t_max = {round(t_max/t_rhi, 1)} t_rhi")

    if arguments.scatter:
        output_folder = arguments.output+"scatter"
        input_output.create_directory(output_folder)

        for stars, time in zip(snapshots, metrics['times']):
            print(f"Plotting t={time} out of {t_max}", end="\r")
            if binaries_found:
                scatterplot(stars, time, output_folder, first_binary)
            else:
                scatterplot(stars, time, output_folder)

    if arguments.scatter_core:
        output_folder = arguments.output+"core_scatter"
        input_output.create_directory(output_folder)

        for stars, time in zip(snapshots, metrics['times']):
            print(f"Plotting t={time} out of {t_max}", end="\r")
            if binaries_found:
                scatterplot(stars, time, output_folder,
                            first_binary, metrics_by_time)
            else:
                scatterplot(stars, time, output_folder,
                            None, metrics_by_time)

    radiiplot(metrics, arguments)
