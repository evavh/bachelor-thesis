from amuse.units import nbody_system
from matplotlib import pyplot
import numpy

import formulas
import main


def get_xylim(metrics_by_time, time, radius_key):
    radius = metrics_by_time[time.number][radius_key]
    radius = radius.value_in(nbody_system.length)

    density_centre = metrics_by_time[time.number]['density_centre']
    density_centre = density_centre.value_in(nbody_system.length)

    x_centre = density_centre[0]
    y_centre = density_centre[1]

    neg_xlim = x_centre - radius
    neg_ylim = y_centre - radius
    pos_xlim = x_centre + radius
    pos_ylim = y_centre + radius

    return (neg_xlim, neg_ylim, pos_xlim, pos_ylim)


def scatter(snapshot, time, output_folder, xylims, rvir,
            first_binary_ids=None):
    if not output_folder.endswith("/"):
        output_folder += "/"

    first_binary = None
    if first_binary_ids is not None:
        first_binary = main.ids_to_stars(snapshot, first_binary_ids)

    x_of_stars = snapshot.x/rvir
    y_of_stars = snapshot.y/rvir

    colours = []
    sizes = []
    for star in snapshot:
        if first_binary is not None:
            if star == first_binary[0] or star == first_binary[1]:
                colours.append('red')
                sizes.append(10)
            else:
                colours.append('blue')
                sizes.append(0.1)
        else:
            colours.append('blue')
            sizes.append(0.1)

    pyplot.scatter(x_of_stars, y_of_stars, s=sizes, c=colours)
    pyplot.xlabel("$x/r_v$")
    pyplot.ylabel("$y/r_v$")

    neg_xlim, neg_ylim, pos_xlim, pos_ylim = xylims
    pyplot.xlim(neg_xlim, pos_xlim)
    pyplot.ylim(neg_ylim, pos_ylim)

    axes = pyplot.gca()
    axes.set_aspect('equal')

    N = len(snapshot)
    kT = 1/(6*N)
    if first_binary is None:
        pyplot.suptitle(f"t = {time.number} (nbody time)")
    else:
        Eb = formulas.binding_energy(first_binary[0], first_binary[1])
        pyplot.suptitle((f"t = {time.number} (nbody time), "
                         f"Eb = {round(Eb.number/kT, 1)} kT"))
    pyplot.savefig(output_folder+str(time)+".svg", format='svg')
    pyplot.clf()


def radii(metrics, config):
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
    pyplot.savefig(config.output+"radii.svg", format='svg')
    pyplot.clf()


def number_of_binaries(metrics, config, t_rhi):
    number_of_binaries = []

    for binaries in metrics['binaries']:
        number_of_binaries.append(len(binaries))

    times = metrics['times']/t_rhi

    pyplot.plot(times, number_of_binaries)
    pyplot.xlabel("$t/t_rh,i$")
    pyplot.ylabel("$n_b$")

    pyplot.gca().set_aspect('equal')
    pyplot.savefig(config.output+"number_of_binaries.svg", format='svg')
    pyplot.clf()


def integration_time(metrics, config):
    times = metrics['times'].value_in(nbody_system.time)

    pyplot.plot(times[1:], metrics['integration_time'][1:])
    pyplot.xlabel("t")
    pyplot.ylabel("Integration time [s]")

    pyplot.savefig(config.output+"integration_time.svg", format='svg')
    pyplot.clf()


def N_core(snapshots, metrics, config):
    density_centres = metrics['density_centre']
    core_radii = metrics['rcore']
    times = metrics['times']

    N_core = []

    for snapshot, density_centre, core_radius in zip(snapshots,
                                                     density_centres,
                                                     core_radii):
        core_stars = main.stars_in_area(snapshot, density_centre, core_radius)
        N_core.append(len(core_stars))

    pyplot.plot(times.number, N_core)
    pyplot.xlabel("t")
    pyplot.ylabel("Number of stars in the core")

    pyplot.savefig(config.output+"N_core.svg", format='svg')
    pyplot.clf()


def work_function(work_for_star, metrics, config, Eb, start, stop):
    times = metrics['times'][start:stop]

    for star in work_for_star:
        pyplot.plot(times.value_in(nbody_system.time), work_for_star[star],
                    label=str(star))

    pyplot.plot(times.value_in(nbody_system.time),
                Eb[start:stop], label="The binary")
    pyplot.plot(times.value_in(nbody_system.time),
                numpy.sum(list(work_for_star.values()), axis=0),
                label="Total work by top stars")

    pyplot.xlabel("t")
    pyplot.ylabel("Work on / Eb of binary [kT]")
    pyplot.legend()

    pyplot.savefig(config.output+"work_function.svg", format='svg')
    pyplot.clf()
