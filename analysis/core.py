import datetime
import numpy

import helpers
import formulas
import plotting
import input_output

from amuse.units import nbody_system


def find_first_binary(data, t_rhi):
    metrics = data.metrics

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
                   f"{round(t_bin/t_rhi, 2)} t_rhi."))
            if len(binaries) > 1:
                print((f"NOTE: {len(binaries)-1} more were found at this"
                       " time!"))

            return first_binary_ids, t_bin

    print("No binaries found.")
    return None, None


def scan_binary_metrics(data, first_binary_ids, t_rhi):
    t_bin_0 = None
    Eb = []

    for snapshot, time in zip(data.snapshots, data.metrics['times']):
        first_binary = helpers.ids_to_stars(snapshot, first_binary_ids)
        if formulas.binding_energy(*first_binary) \
                > 0 | nbody_system.energy:
            if t_bin_0 is None:
                t_bin_0 = time
                print((f"It has formed by t = {t_bin_0} = "
                       f"{round(t_bin_0/t_rhi, 2)} t_rhi."))

        Eb.append(formulas.binding_energy(*first_binary)
                  .value_in(nbody_system.energy))

    kT = 1/(6*data.params.n)
    Eb = numpy.array(Eb)/kT

    return Eb, t_bin_0


def find_core_stars(data, t_bin_0):
    snapshots = data.snapshots
    metrics_by_time = data.by_time()

    density_centre = metrics_by_time[t_bin_0.number]['density_centre']
    core_radius = metrics_by_time[t_bin_0.number]['rcore']

    t_bin_0_index = helpers.time_to_index(t_bin_0.number, data)
    core_stars = helpers.stars_in_area(snapshots[t_bin_0_index],
                                       density_centre,
                                       core_radius)
    print(f"There are {len(core_stars)} stars in the core at t_bin.")

    core_star_ids = helpers.stars_to_ids(core_stars)
    return core_star_ids


# output is in kT
def calculate_work(data, first_binary_ids, work_start_i,
                   work_end_i, core_star_ids):
    kT = 1/(6*data.params.n)

    star_works = {}
    total_star_works = {}
    print("Starting work function calculation.")
    start_of_calc = datetime.datetime.now()

    for star_id in core_star_ids:
        key = star_id.number
        star_works[key], total_star_works[key] = \
            formulas.work_function(data, first_binary_ids,
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


def produce_scatterplots(data, config, first_binary_ids, t_max):
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
            xylims = plotting.get_xylim(data, time, radius_key)
            rvir = data.by_time()[time.number]['rvir']

            if first_binary_ids is not None:
                plotting.scatter(snapshot, time, output_folder, xylims, rvir,
                                 first_binary_ids)
            else:
                plotting.scatter(snapshot, time, output_folder, xylims, rvir)
