from amuse.units import nbody_system

import numpy
import itertools
import datetime
import pickle

import input_output
import formulas
import plotting
from data import Data


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


def find_core_stars(snapshots, metrics, metrics_by_time, t_bin_0):
    density_centre = metrics_by_time[t_bin_0.number]['density_centre']
    core_radius = metrics_by_time[t_bin_0.number]['rcore']

    t_bin_0_index = time_to_index(t_bin_0, metrics)
    core_stars = stars_in_area(snapshots[t_bin_0_index], density_centre,
                               core_radius)
    print(f"There are {len(core_stars)} stars in the core at t_bin.")

    core_star_ids = stars_to_ids(core_stars)
    return core_star_ids


def calculate_work(snapshots, metrics, first_binary_ids, work_start_i,
                   work_end_i):
    star_works = {}
    total_star_works = {}
    print("Starting work function calculation.")
    start_of_calc = datetime.datetime.now()

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

    return star_works, total_star_works


if __name__ == '__main__':
    config = input_output.parse_arguments()

    input_output.create_directory(config.output)
    input_output.create_directory(config.input)

    data = Data(config)
    print('')

    kT = 1/(6*data.params.n)
    t_min = data.metrics['times'][0]
    t_max = data.metrics['times'][-1]
    t_rhi = formulas.t_rh(data.params.n, data.by_time()[0.0]['r50pc'])
    times_crc = data.metrics['t_crc']

    if data.params.variable_delta:
        taus = numpy.cumsum(numpy.full(len(times_crc), 0.01))
    else:
        taus = numpy.cumsum(data.params.delta_t/times_crc)

    print(f"t_max = {t_max.number} = {round(t_max/t_rhi, 1)} t_rhi")

    first_binary_ids, t_bin_10 = find_first_binary(data.snapshots,
                                                   data.metrics,
                                                   t_rhi)
    binaries_found = (first_binary_ids is not None)

    if binaries_found:
        t_bin_0 = None
        Eb = []
        for snapshot, time in zip(data.snapshots, data.metrics['times']):
            first_binary = ids_to_stars(snapshot, first_binary_ids)
            if formulas.binding_energy(*first_binary) \
                    > 0 | nbody_system.energy:
                if t_bin_0 is None:
                    t_bin_0 = time
                    print((f"It has formed by t = {t_bin_0.number} = "
                           f"{round(t_bin_0/t_rhi, 1)} t_rhi."))

            Eb.append(formulas.binding_energy(*first_binary)
                      .value_in(nbody_system.energy))

        Eb = numpy.array(Eb)/kT

        if config.load_work:
            star_works = pickle.load(open(config.output+"star_works.pkl",
                                          'rb'))
            total_star_works = pickle.load(open(config.output +
                                                "total_star_works.pkl", 'rb'))
            work_start_i, work_end_i = pickle.load(open(config.output +
                                                        "work_indexes.pkl",
                                                        'rb'))
        else:
            core_star_ids = find_core_stars(data.snapshots, data.metrics,
                                            data.by_time(), t_bin_0)

            work_start_i = time_to_index(t_min, data.metrics)
            work_end_i = time_to_index(t_bin_10, data.metrics)

            star_works, total_star_works = calculate_work(data.snapshots,
                                                          data.metrics,
                                                          first_binary_ids,
                                                          work_start_i,
                                                          work_end_i)
            pickle.dump(star_works, open(config.output+"star_works.pkl",
                                         "wb"))
            pickle.dump(total_star_works,
                        open(config.output+"total_star_works.pkl", "wb"))
            pickle.dump((work_start_i, work_end_i),
                        open(config.output+"work_indexes.pkl", 'wb'))

        top_stars = dict(itertools.islice(star_works.items(), 10))

    if config.scatter:
        folder_name = "scatter"
        radius_key = 'rvir'

    if config.scatter_core:
        folder_name = "scatter_core"
        radius_key = 'rcore'

    if config.scatter or config.scatter_core:
        output_folder = config.output+folder_name
        input_output.create_directory(output_folder)

        for snapshot, time in zip(data.snapshots, data.metrics['times']):
            print(f"Plotting t={time} out of {t_max}", end="\r")
            xylims = plotting.get_xylim(data.by_time(), time, radius_key)
            rvir = data.by_time()[time.number]['rvir']

            if binaries_found:
                plotting.scatter(snapshot, time, output_folder, xylims, rvir,
                                 first_binary_ids)
            else:
                plotting.scatter(snapshot, time, output_folder, xylims, rvir)

    plotting.radii(data.metrics, config)
    plotting.number_of_binaries(data.metrics, config, t_rhi)
    plotting.integration_time(data.metrics, config)
    plotting.N_core(data.snapshots, data.metrics, config)
    plotting.work_function(top_stars, data.metrics, config, Eb,
                           work_start_i, work_end_i)
