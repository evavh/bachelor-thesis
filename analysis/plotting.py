from amuse.units import nbody_system
from matplotlib import pyplot

import formulas


def get_xylim(metrics_by_time, time, radius_key):
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


def scatter(stars, time, output_folder, xylims, first_binary=None):
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

    N = len(stars)
    kT = 1/(6*N)
    Eb = formulas.binding_energy(first_binary[0], first_binary[1])
    pyplot.suptitle((f"t = {time.number} (nbody time),"
                     f"Eb = {round(Eb.number/kT, 1)} kT"))
    pyplot.savefig(output_folder+str(time)+".svg", format='svg')
    pyplot.clf()


def radii(metrics, arguments):
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
