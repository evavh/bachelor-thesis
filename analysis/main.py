from amuse.units import nbody_system

import numpy

import input_output
import formulas
import plotting


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
        folder_name = "scatter"
        radius_key = 'rvir'

    if arguments.scatter_core:
        folder_name = "scatter_core"
        radius_key = 'rcore'

    if arguments.scatter or arguments.scatter_core:
        output_folder = arguments.output+folder_name
        input_output.create_directory(output_folder)

        for stars, time in zip(snapshots, metrics['times']):
            print(f"Plotting t={time} out of {t_max}", end="\r")
            xylims = plotting.get_xylim(metrics_by_time, time, radius_key)

            if binaries_found:
                plotting.scatter(stars, time, output_folder, xylims,
                                 first_binary)
            else:
                plotting.scatter(stars, time, output_folder, xylims)

    plotting.radii(metrics, arguments)
