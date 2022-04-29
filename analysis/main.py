from amuse.units import nbody_system

import numpy
import itertools
import pickle
import sys

import helpers
import core
import input_output
import formulas
import plotting
from data import Data


def flushed_print(string):
    print(string)
    sys.stdout.flush()


if __name__ == '__main__':
    flushed_print("Analysis script started.")
    config = input_output.parse_arguments()

    input_output.create_directory(config.output)
    input_output.create_directory(config.input)

    flushed_print("Starting data reading.")
    data = Data(config)
    print('')

    if not data.params.variable_delta and config.calc_work:
        print("Work calc uses variable delta, invalid otherwise!")

    kT = 1/(6*data.params.n)
    t_min = data.metrics['times'][0]
    t_max = data.metrics['times'][-1]
    t_rhi = formulas.t_rh(data)
    times_crc = data.metrics['t_crc']

    if data.params.variable_delta:
        taus = numpy.cumsum(numpy.full(len(times_crc), 0.01))
    else:
        taus = numpy.cumsum(data.params.delta_t/times_crc)

    print(f"t_max = {t_max.number} = {round(t_max/t_rhi, 2)} t_rhi")

    first_binary_ids, t_bin_10 = core.find_first_binary(data, t_rhi)

    core.produce_scatterplots(data, config, first_binary_ids, t_max)
    plotting.radii(data, config)
    plotting.number_of_binaries(data, config, t_rhi)
    plotting.integration_time(data, config)

    if not config.fast_plot:
        plotting.N_core(data, config)

    if first_binary_ids is None:
        quit()

    if config.fast_plot:
        Eb, t_bin_0 = pickle.load(open(config.output+"Eb-t_bin_0.pkl", 'rb'))
        print(f"The binary has formed by {t_bin_0}")
    else:
        Eb, t_bin_0 = core.scan_binary_metrics(data, first_binary_ids, t_rhi)
        pickle.dump((Eb, t_bin_0), open(config.output+"Eb-t_bin_0.pkl", "wb"))

    if config.load_work:
        star_works = pickle.load(open(config.output+"star_works.pkl",
                                      'rb'))
        total_star_works = pickle.load(open(config.output +
                                            "total_star_works.pkl", 'rb'))
        work_start_i, work_end_i = pickle.load(open(config.output +
                                                    "work_indexes.pkl",
                                                    'rb'))
    elif config.calc_work:
        core_star_ids = core.find_core_stars(data, t_bin_0)

        work_start_i = helpers.time_to_index(t_bin_0.number - 1,
                                             data)
        work_end_i = helpers.time_to_index(t_bin_10.number, data)
        print(f"Calculating work from {t_bin_0.number-1} to {t_bin_10}.")

        star_works, total_star_works = core.calculate_work(data,
                                                           first_binary_ids,
                                                           work_start_i,
                                                           work_end_i,
                                                           core_star_ids)
        pickle.dump(star_works, open(config.output+"star_works.pkl",
                                     "wb"))
        pickle.dump(total_star_works,
                    open(config.output+"total_star_works.pkl", "wb"))
        pickle.dump((work_start_i, work_end_i),
                    open(config.output+"work_indexes.pkl", 'wb'))
    else:
        quit()

    top_stars = dict(itertools.islice(star_works.items(), 10))
    plotting.work_function(top_stars, data, config, Eb,
                           work_start_i, work_end_i)
