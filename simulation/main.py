import input_output
import setting_up
from gravity import Gravity
import core

import sys
import datetime

from amuse.units import nbody_system
from amuse.community.brutus.interface import Brutus
from amuse.community.ph4.interface import ph4

from amuse.io import write_set_to_file
from amuse.rfi.core import is_mpd_running


def flushed_print(string):
    print(string)
    sys.stdout.flush()


if __name__ == '__main__':
    flushed_print("Simulation script started.")

    assert is_mpd_running()

    metrics = setting_up.initialize_metrics()
    CONSTS = {'accuracy': 0.01,
              'initial_accuracy': 0.0025,
              'epsilon_squared': 0 | nbody_system.length**2}
    params = input_output.parse_arguments()
    params.random_seed = setting_up.set_random_seed(params.random_seed)

    input_output.create_directory(params.output_folder)

    input_output.pickle_object(params, "parameters.pkl", params)
    input_output.pickle_object(CONSTS, "constants.pkl", params)

    kT = 1/(6*params.n)
    minimum_Eb = params.minimum_Eb_kT * kT

    flushed_print("\nStarting simulation setup.")
    stars, time = setting_up.initialize_stars(params, CONSTS)
    if params.reverse:
        stars = setting_up.reverse_velocities(stars)
        gravity = Gravity(Brutus, CONSTS, time, stars, params)
    else:
        gravity = Gravity(ph4, CONSTS, time, stars, params)

    binaries = []
    binding_energies = []
    integration_time = 0

    E0 = gravity.total_energy()
    max_dE = 0

    print("")

    while True:

        if params.reverse:
            reverse_time = params.start_time - (time - params.start_time)

            print(f"Saving metrics and snapshot at t={reverse_time.number}")
            flushed_print((f"This is {time.number} in integrator time, "
                           f"with start time {params.start_time.number}"))
            metrics = core.update_metrics(metrics, reverse_time, stars,
                                          gravity.integrator,
                                          binaries, binding_energies, kT,
                                          integration_time)
            input_output.pickle_object(stars, f"snapshot_{reverse_time}.pkl",
                                       params)
        else:
            reverse_time = None

            flushed_print(f"Saving metrics and snapshot at t={time.number}")
            metrics = core.update_metrics(metrics, time, stars,
                                          gravity.integrator,
                                          binaries,
                                          binding_energies, kT,
                                          integration_time)
            input_output.pickle_object(stars, f"snapshot_{time}.pkl", params)

        input_output.pickle_object(metrics, "cluster_metrics.pkl", params)

        E_tot = gravity.total_energy()

        dE = (E_tot - E0)/E0
        if abs(dE) > max_dE:
            max_dE = dE
        print(f"dE = {dE}")

        if params.t_end is not None and time >= params.t_end:
            break

        flushed_print(("\nDone saving, starting integration at "
                       f"t={time.number} (integrator time)"))
        integration_start_time = datetime.datetime.now()

        t_crc = metrics['t_crc'][-1]
        if params.variable_delta:
            time += 0.01*t_crc
            print(f"delta_t is 0.01*t_crc = {0.01*t_crc} for this iteration")
        else:
            time += params.delta_t
            print(f"delta_t is {params.delta_t} for this iteration")

        while gravity.time() < time:
            gravity.evolve_model(time)

        gravity.copy_from_worker()

        print(f"Finished integrating until t={time.number} (integrator time)")
        integration_time = datetime.datetime.now() - integration_start_time
        integration_time = integration_time.total_seconds()
        print(f"{integration_time}s elapsed.")

        flushed_print("Starting binary finding.")

        binaries, binding_energies = core.find_binaries(stars, minimum_Eb)
        if len(binaries) > 0 and params.t_end is None:
            params.t_end = time + (5 | nbody_system.time)
            print(f"Set t_end to {params.t_end}")

    gravity.stop()
    print(f"\nMaximum dE = {max_dE}")

    write_set_to_file(stars, params.output_folder+"final_state.csv", "csv")
    input_output.round_csv(params.output_folder+"final_state.csv", 3)
