from amuse.units import nbody_system
from matplotlib import pyplot
import numpy

import formulas
import helpers


def get_xylim(data, time, radius_key):
    metrics_by_time = data.by_time()
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


def scatter(config, snapshot, time, output_folder, xylims, rvir,
            first_binary_ids=None):
    if not output_folder.endswith("/"):
        output_folder += "/"

    x_of_stars = snapshot.x/rvir
    y_of_stars = snapshot.y/rvir

    pyplot.scatter(x_of_stars, y_of_stars, s=0.1, c='k')

    first_binary = None
    if first_binary_ids is not None:
        first_binary = helpers.ids_to_stars(snapshot, first_binary_ids)

        x_of_binary = first_binary.x/rvir
        y_of_binary = first_binary.y/rvir

        pyplot.scatter(x_of_binary, y_of_binary, s=30.0,
                       fc='w', ec='k', marker='o')

        for x, y, id in zip(x_of_binary, y_of_binary, first_binary_ids):
            pyplot.annotate(round(id), (x + 0.01, y + 0.01), fontsize=5)

    interesting_ids = config.interesting
    if interesting_ids is not None:
        interesting_stars = helpers.ids_to_stars(snapshot, interesting_ids)

        x_of_interesting = interesting_stars.x/rvir
        y_of_interesting = interesting_stars.y/rvir

        for x, y, id in zip(x_of_interesting, y_of_interesting, interesting_ids):
            pyplot.annotate(round(id), (x, y), fontsize=5)

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
        pyplot.suptitle(f"t = {round(time.number, 5)} (n-body time)")
    else:
        Eb = formulas.binding_energy(first_binary[0], first_binary[1])
        pyplot.suptitle((f"t = {round(time.number, 5)} (n-body time)\n"
                         f"Binary ({round(first_binary_ids[0])},"
                         f"{round(first_binary_ids[1])}) Eb ="
                         f" {round(Eb.number/kT, 1)} kT"))
    pyplot.savefig(output_folder+str(time)+".svg", format='svg')
    pyplot.clf()


def radii(data, config):
    metrics = data.metrics
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


def number_of_binaries(data, config, t_rhi):
    metrics = data.metrics
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


def integration_time(data, config):
    metrics = data.metrics
    times = metrics['times'].value_in(nbody_system.time)

    pyplot.plot(times[1:], metrics['integration_time'][1:])
    pyplot.xlabel("t")
    pyplot.ylabel("Integration time [s]")

    pyplot.savefig(config.output+"integration_time.svg", format='svg')
    pyplot.clf()


def N_core(data, config):
    snapshots = data.snapshots
    metrics = data.metrics

    density_centres = metrics['density_centre']
    core_radii = metrics['rcore']
    times = metrics['times']

    N_core = []

    for snapshot, density_centre, core_radius in zip(snapshots,
                                                     density_centres,
                                                     core_radii):
        core_stars = helpers.stars_in_area(snapshot, density_centre,
                                           core_radius)
        N_core.append(len(core_stars))

    pyplot.plot(times.number, N_core)
    pyplot.xlabel("t")
    pyplot.ylabel("Number of stars in the core")

    pyplot.savefig(config.output+"N_core.svg", format='svg')
    pyplot.clf()


def Eb(data, config, Eb):
    print("Plotting Eb")
    if data is None:
        times = range(len(Eb))
    else:
        times = data.metrics['times'].number
    pyplot.plot(times, Eb, color='k', linewidth=2.0)
    pyplot.axhline(y=config.energy_threshold, color='k', linestyle='--')
    pyplot.xlabel("t")
    pyplot.ylabel("Eb of the binary")
    pyplot.savefig(config.output+"Eb.svg", format='svg')
    pyplot.clf()
    print("Finished plotting Eb")


def work_function(work_for_star, data, config, Eb, start, stop, taus=None):
    metrics = data.metrics
    if taus is None:
        times = metrics['times'][start:stop].number
    else:
        times = taus[start:stop]

    for star in work_for_star:
        pyplot.plot(times, work_for_star[star],
                    label=str(star), color='k', linewidth=0.7)
    pyplot.plot(times,
                Eb[start:stop], label="The binary", color='k', linewidth=2.0)
    # pyplot.plot(times,
    #             numpy.sum(list(work_for_star.values()), axis=0),
    #             label="Total work by top stars", color='k', linewidth=0.7)

    pyplot.xlabel(r"$\tau$")
    pyplot.ylabel("Work on / Eb of binary [kT]")
    pyplot.legend()

    pyplot.savefig(config.output+"work_function.svg", format='svg')
    pyplot.clf()
