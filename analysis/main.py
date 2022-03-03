from amuse.units import nbody_system
from matplotlib import pyplot

import numpy

import input_output
import formulas


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


def calculate_xylim(metrics_by_time, time, radius_key):
    time_number = time.value_in(nbody_system.time)

    radius = metrics_by_time[time_number][radius_key]
    radius = radius.value_in(nbody_system.length)

    density_centre = metrics_by_time[time_number]['density_centre']
    density_centre = density_centre.value_in(nbody_system.length)

    x_centre = density_centre[0]
    y_centre = density_centre[1]

    neg_xlim = x_centre - radius/2
    neg_ylim = y_centre - radius/2
    pos_xlim = x_centre + radius/2
    pos_ylim = y_centre + radius/2

    return (neg_xlim, neg_ylim, pos_xlim, pos_ylim)


def scatterplot(stars, time, output_folder, xylims, first_binary=None):
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

    neg_xlim, neg_ylim, pos_xlim, pos_ylim = xylims
    pyplot.xlim(neg_xlim, pos_xlim)
    pyplot.ylim(neg_ylim, pos_ylim)

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
    t_rhi = formulas.t_rh(params.n, metrics['r50pc'][0], nbody_system.G,
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
            xylims = calculate_xylim(metrics_by_time, time, 'rvir')

            if binaries_found:
                scatterplot(stars, time, output_folder, xylims, first_binary)
            else:
                scatterplot(stars, time, output_folder, xylims)

    if arguments.scatter_core:
        output_folder = arguments.output+"core_scatter"
        input_output.create_directory(output_folder)

        for stars, time in zip(snapshots, metrics['times']):
            print(f"Plotting t={time} out of {t_max}", end="\r")
            xylims = calculate_xylim(metrics_by_time, time, 'rcore')

            if binaries_found:
                scatterplot(stars, time, output_folder, xylims, first_binary)
            else:
                scatterplot(stars, time, output_folder, xylims)

    radiiplot(metrics, arguments)
