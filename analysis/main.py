from amuse.units import nbody_system

import numpy
import itertools
import datetime

import input_output
import formulas
import plotting


def ids_to_stars(snapshot, ids):
    def criterium(id): return id in ids
    return snapshot.select(criterium, ['id'])


def stars_to_ids(stars):
    ids = []
    for star in stars:
        ids.append(star.id)
    return ids


def time_to_index(time, metrics):
    index = None
    for i, t in zip(range(len(metrics['times'])), metrics['times']):
        if t >= time and index is None:
            index = i

    return index


def limits_from_radius(radius, density_centre):
    x_min = density_centre[0] - radius
    x_max = density_centre[0] + radius

    y_min = density_centre[1] - radius
    y_max = density_centre[1] + radius

    z_min = density_centre[2] - radius
    z_max = density_centre[2] + radius

    return ((x_min, x_max), (y_min, y_max), (z_min, z_max))


def stars_in_area(snapshot, density_centre, radius):
    limits = limits_from_radius(radius, density_centre)

    def within_limits(x, y, z):
        x_lims, y_lims, z_lims = limits
        x_min, x_max = x_lims
        y_min, y_max = y_lims
        z_min, z_max = z_lims

        return (x > x_min and x < x_max and
                y > y_min and y < y_max and
                z > z_min and z < z_max)

    return snapshot.select(within_limits, ['x', 'y', 'z'])


def find_first_binary(snapshots, metrics, t_rhi):
    binaries_found = False
    first_binary = None
    first_binary_ids = None
    t_bin = None

    for time, binaries in zip(metrics['times'], metrics['binaries']):
        if binaries != [] and not binaries_found:
            binaries_found = True
            first_binary = binaries[0]
            first_binary_ids = (first_binary[0].id.number,
                                first_binary[1].id.number)
            t_bin = time

            print(f"The first binary is {first_binary_ids}")
            print((f"It has reached 10kT by t = {t_bin.number} = "
                   f"{round(t_bin/t_rhi, 1)} t_rhi."))
            if len(binaries) > 1:
                print((f"NOTE: {len(binaries)-1} more were found at this"
                       " time!"))

            return first_binary_ids, t_bin

    print("No binaries found.")
    return None, None


if __name__ == '__main__':
    arguments = input_output.parse_arguments()

    input_output.create_directory(arguments.output)
    input_output.create_directory(arguments.input)

    snapshots, consts, params, metrics, og_metrics = \
        input_output.load_data(arguments)
    print('')

    metrics_by_time = input_output.metrics_to_time_key(metrics)

    if og_metrics is not None:
        t0_metrics = input_output.metrics_to_time_key(og_metrics)[0.0]
        metrics_by_time[0.0] = t0_metrics

    if 0.0 not in metrics_by_time:
        print("t=0.0 not found, please provide original run data with -og")

    kT = 1/(6*params.n)
    t_min = metrics['times'][0]
    t_max = metrics['times'][-1]
    t_rhi = formulas.t_rh(params.n, metrics_by_time[0.0]['r50pc'],
                          nbody_system.G, 1 | nbody_system.mass)
    times_crc = metrics['t_crc']

    if params.variable_delta:
        taus = numpy.cumsum(numpy.full(len(times_crc), 0.01))
    else:
        taus = numpy.cumsum(params.delta_t/times_crc)

    first_binary_ids, t_bin_10 = find_first_binary(snapshots, metrics, t_rhi)
    binaries_found = (first_binary_ids is not None)

    t_bin_0 = None
    Eb = []
    for snapshot, time in zip(snapshots, metrics['times']):
        first_binary = ids_to_stars(snapshot, first_binary_ids)
        if formulas.binding_energy(*first_binary) > 0 | nbody_system.energy:
            if t_bin_0 is None:
                t_bin_0 = time
                print((f"It has formed by t = {t_bin_0.number} = "
                       f"{round(t_bin_0/t_rhi, 1)} t_rhi."))

        Eb.append(formulas.binding_energy(*first_binary)
                  .value_in(nbody_system.energy))

    Eb = numpy.array(Eb)/kT

    print(f"t_max = {round(t_max/t_rhi, 1)} t_rhi")

    density_centre = metrics_by_time[t_bin_0.number]['density_centre']
    core_radius = metrics_by_time[t_bin_0.number]['rcore']

    t_bin_0_index = time_to_index(t_bin_0, metrics)
    core_stars = stars_in_area(snapshots[t_bin_0_index], density_centre,
                               core_radius)
    print(f"There are {len(core_stars)} stars in the core at t_bin.")

    core_star_ids = stars_to_ids(core_stars)
    star_works = {}
    total_star_works = {}
    print("Starting work function calculation.")
    start_of_calc = datetime.datetime.now()

    work_start = -14
    t_min = metrics['times'][work_start-1]

    work_start_i = time_to_index(t_min, metrics)
    work_end_i = time_to_index(t_bin_10, metrics)

    for star_id in core_star_ids:
        key = star_id.number
        star_works[key], total_star_works[key] = \
            formulas.work_function(snapshots, metrics, first_binary_ids,
                                   star_id, work_start_i, work_end_i)
        star_works[key] /= kT
        total_star_works[key] /= kT
    calc_time_s = datetime.datetime.now() - start_of_calc
    calc_time_s = calc_time_s.total_seconds()
    print(f"Work function calculation finished after {calc_time_s} s.")

    star_works = dict(sorted(star_works.items(),
                             key=lambda x: abs(total_star_works[x[0]]),
                             reverse=True))
    print(total_star_works)
    top_stars = dict(itertools.islice(star_works.items(), 10))

    if arguments.scatter:
        folder_name = "scatter"
        radius_key = 'rvir'

    if arguments.scatter_core:
        folder_name = "scatter_core"
        radius_key = 'rcore'

    if arguments.scatter or arguments.scatter_core:
        output_folder = arguments.output+folder_name
        input_output.create_directory(output_folder)

        for snapshot, time in zip(snapshots, metrics['times']):
            print(f"Plotting t={time} out of {t_max}", end="\r")
            xylims = plotting.get_xylim(metrics_by_time, time, radius_key)
            rvir = metrics_by_time[time.number]['rvir']

            if binaries_found:
                plotting.scatter(snapshot, time, output_folder, xylims, rvir,
                                 first_binary_ids)
            else:
                plotting.scatter(snapshot, time, output_folder, xylims, rvir)

    plotting.radii(metrics, arguments)
    plotting.number_of_binaries(metrics, arguments, t_rhi)
    plotting.integration_time(metrics, arguments)
    plotting.N_core(snapshots, metrics, arguments)
    plotting.work_function(top_stars, metrics, arguments, Eb, work_start_i,
                           work_end_i)
