from amuse.units import nbody_system

import numpy

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


def binary_wont_dissolve(snapshots, starting_snap, candidate_ids):
    for snapshot in snapshots[starting_snap:]:
        stars = ids_to_stars(snapshot, candidate_ids)
        if formulas.binding_energy(stars[0], stars[1])\
                < 0 | nbody_system.energy:
            return False
    return True


def find_first_binary(snapshots, metrics, t_rhi):
    for index, time, binaries in zip(range(len(metrics['times'])),
                                     metrics['times'], metrics['binaries']):
        for candidate in binaries:
            candidate_ids = stars_to_ids(candidate)
            if binary_wont_dissolve(snapshots, index, candidate_ids):
                first_binary_ids = candidate_ids
                t_bin = time

                print(f"The first hard enough binary is {first_binary_ids}")
                print((f"It was found at t = {t_bin.number} = "
                       f"{round(t_bin/t_rhi, 1)} t_rhi."))

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

    t_max = metrics['times'][-1]
    t_rhi = formulas.t_rh(params.n, metrics_by_time[0.0]['r50pc'],
                          nbody_system.G, 1 | nbody_system.mass)
    times_crc = metrics['t_crc']

    if params.variable_delta:
        taus = numpy.cumsum(numpy.full(len(times_crc), 0.01))
    else:
        taus = numpy.cumsum(params.delta_t/times_crc)

    E_tot_unit = numpy.array(metrics['potential_energy']
                             + metrics['kinetic_energy']
                             + metrics['total_binary_energy'])
    E_tot = []
    for E in E_tot_unit:
        E_tot.append(E.value_in(nbody_system.energy))
    E_tot = numpy.array(E_tot)

    dEs = (E_tot[0] - E_tot)/E_tot

    for t, dE in zip(metrics['times'], dEs):
        print(t.number, dE)

    first_binary_ids, t_bin = find_first_binary(snapshots, metrics, t_rhi)
    binaries_found = (first_binary_ids is not None)
    if binaries_found:
        first_binary = ids_to_stars(snapshots[0], first_binary_ids)

    print(f"t_max = {round(t_max/t_rhi, 1)} t_rhi")

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
