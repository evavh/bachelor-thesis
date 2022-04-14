from amuse.units import nbody_system

import numpy
import itertools
import pickle

import helpers
import core
import input_output
import formulas
import plotting
from data import Data


if __name__ == '__main__':
    config = input_output.parse_arguments()

    input_output.create_directory(config.output)
    input_output.create_directory(config.input)

    data = Data(config)
    print('')

    kT = 1/(6*data.params.n)
    t_min = data.metrics['times'][0]
    t_max = data.metrics['times'][-1]
    t_rhi = formulas.t_rh(data)
    times_crc = data.metrics['t_crc']

    if data.params.variable_delta:
        taus = numpy.cumsum(numpy.full(len(times_crc), 0.01))
    else:
        taus = numpy.cumsum(data.params.delta_t/times_crc)

    print(f"t_max = {t_max.number} = {round(t_max/t_rhi, 1)} t_rhi")

    first_binary_ids, t_bin_10 = core.find_first_binary(data, t_rhi)

    if first_binary_ids is not None:
        t_bin_0 = None
        Eb = []
        for snapshot, time in zip(data.snapshots, data.metrics['times']):
            first_binary = helpers.ids_to_stars(snapshot, first_binary_ids)
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
            core_star_ids = core.find_core_stars(data, t_bin_0)

            work_start_i = helpers.time_to_index(t_min, data)
            work_end_i = helpers.time_to_index(t_bin_10, data)

            star_works, total_star_works = core.calculate_work(data,
                                                               first_binary_ids,
                                                               work_start_i,
                                                               work_end_i,
                                                               core_star_ids,
                                                               kT)
            pickle.dump(star_works, open(config.output+"star_works.pkl",
                                         "wb"))
            pickle.dump(total_star_works,
                        open(config.output+"total_star_works.pkl", "wb"))
            pickle.dump((work_start_i, work_end_i),
                        open(config.output+"work_indexes.pkl", 'wb'))

        top_stars = dict(itertools.islice(star_works.items(), 10))

    core.produce_scatterplots(data, config, first_binary_ids, t_max)

    plotting.radii(data, config)
    plotting.number_of_binaries(data, config, t_rhi)
    plotting.integration_time(data, config)
    plotting.N_core(data, config)
    plotting.work_function(top_stars, data, config, Eb,
                           work_start_i, work_end_i)
