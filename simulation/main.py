import file_io
import setting_up
from gravity import Gravity

import math
import numpy
import argparse
import sys
import datetime

from amuse.units import nbody_system
from amuse.units import units
from amuse.ext.LagrangianRadii import LagrangianRadii
from amuse.community.brutus.interface import Brutus
from amuse.community.ph4.interface import ph4

from amuse.io import write_set_to_file
from amuse.rfi.core import is_mpd_running

from matplotlib import pyplot


def flushed_print(string):
    print(string)
    sys.stdout.flush()


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--delta_t", help="time between snapshots",
                        default=1.0, type=float)
    parser.add_argument("-v", "--variable_delta", help="delta_t = 0.01t_crc",
                        action='store_true')
    parser.add_argument("-n", help="number of stars to simulate",
                        default=100, type=int)
    parser.add_argument("-s", "--random_seed", help="random number seed",
                        default=None, type=int)
    parser.add_argument("-t", "--t_end", help="time to force end the sim",
                        default=None, type=float)
    parser.add_argument("-T", "--start_time", help="snapshot time to load",
                        default=None, type=float)
    parser.add_argument("-o", "--output_folder", help="where to put output",
                        default="simulation/output/", type=str)
    parser.add_argument("-i", "--snapshot_input", help="snaps to start from",
                        default=None, type=str)
    parser.add_argument("-b", "--minimum_Eb_kT", help="minimum binding E / kT",
                        default=10, type=float)
    parser.add_argument("-r", "--reverse", help="run the sim in reverse",
                        action='store_true')

    params = parser.parse_args()

    if params.reverse:
        print("Running a reversed simulation.")
        assert params.t_end is not None, \
            "When reversed there must be an end time."
        assert params.start_time is not None, \
            "When reversed there must be a start time."
        assert params.start_time > params.t_end, \
            "When reversed start time must be later than end time."

        params.t_end = params.start_time - params.t_end \
            + params.start_time
        print(f"Set end time to {params.t_end} for the integrator's sake.")
    elif params.start_time is not None and params.t_end is not None:
        assert params.start_time < params.t_end, \
            "When not reversed start time must be before end time."

    params.delta_t = params.delta_t | nbody_system.time
    if params.t_end is not None:
        params.t_end = params.t_end | nbody_system.time
    if params.start_time is not None:
        assert params.snapshot_input is not None, \
            "When there is a start time there must be snapshots."
        params.start_time = params.start_time | nbody_system.time

    if not params.output_folder.endswith('/'):
        params.output_folder += '/'

    return params


def set_random_seed(random_seed):
    if random_seed is None:
        numpy.random.seed()
        random_seed = numpy.random.randint(1, pow(2, 31)-1)
    numpy.random.seed(random_seed)
    return random_seed


def scatterplot(stars, filename):
    pyplot.scatter(stars.x.value_in(nbody_system.length),
                   stars.y.value_in(nbody_system.length),
                   s=10*stars.mass/stars.mass.max())
    pyplot.xlabel("x")
    pyplot.ylabel("y")
    pyplot.savefig(filename)
    pyplot.clf()


def calculate_t_crc(core_density):
    G = nbody_system.G
    under_sqrt = 4*math.pi*G*core_density/3
    return 1/under_sqrt.sqrt()


def binding_energy(star1, star2):
    mu = star1.mass*star2.mass/(star1.mass+star2.mass)
    dr = (star2.position - star1.position).lengths()
    dv = (star2.velocity - star1.velocity).lengths()
    Eb = nbody_system.G*star1.mass*star2.mass/dr - 0.5*mu*dv*dv

    return Eb


def find_binaries(stars, minimum_Eb):
    G = nbody_system.G
    binding_energies = []
    binaries = []

    for star in stars:
        mu = star.mass*stars.mass/(star.mass+stars.mass)
        dr = (stars.position - star.position).lengths()
        dv = (stars.velocity - star.velocity).lengths()
        Eb = G*star.mass*stars.mass/dr - 0.5*mu*dv*dv

        # find index of second largest binding energy
        # (largest is binding energy to self, which is infinite)
        # .number removes units
        maxEb_index = numpy.argpartition(-Eb.number, 2)[1]
        maxEb = Eb.number[maxEb_index]
        partner = stars[maxEb_index]
        if maxEb > minimum_Eb and star.id < partner.id:
            binding_energies.append(maxEb)
            binaries.append((star, partner))

    return binaries, binding_energies


def update_metrics(metrics, time, stars, gravity, binaries,
                   binding_energies, kT, integration_time):
    metrics["times"].append(time)

    metrics["binaries"].append(binaries)

    binary_ids = []
    for binary in binaries:
        star1 = binary[0]
        star2 = binary[1]
        binary_ids.append((star1.id, star2.id))
    metrics["binary_ids"].append(binary_ids)

    metrics["binding_energies_kT"].append([x/kT for x in binding_energies])

    metrics["rvir"].append(stars.virial_radius())

    Rl = LagrangianRadii(stars, massf=[0.1, 0.5, 1.0] | units.none)
    metrics["r10pc"].append(Rl[0])
    metrics["r50pc"].append(Rl[1])

    density_centre, core_radius, core_density = \
        gravity.particles.densitycentre_coreradius_coredens()
    metrics["density_centre"].append(density_centre)
    metrics["rcore"].append(core_radius)
    metrics["core_density"].append(core_density)

    metrics["total_mass"].append(gravity.total_mass)
    metrics["potential_energy"].append(gravity.potential_energy)
    metrics["kinetic_energy"].append(gravity.kinetic_energy)

    metrics['t_crc'].append(calculate_t_crc(core_density))

    metrics['integration_time'].append(integration_time)

    return metrics


if __name__ == '__main__':
    flushed_print("Simulation script started.")

    assert is_mpd_running()

    metrics = setting_up.initialize_metrics()
    CONSTS = {'accuracy': 0.01,
              'initial_accuracy': 0.0025,
              'epsilon_squared': 0 | nbody_system.length**2,
              'bs_tolerance': 1e-10}
    params = parse_arguments()
    params.random_seed = set_random_seed(params.random_seed)

    file_io.create_directory(params.output_folder)

    file_io.pickle_object(params, "parameters.pkl", params)
    file_io.pickle_object(CONSTS, "constants.pkl", params)

    kT = 1/(6*params.n)
    minimum_Eb = params.minimum_Eb_kT * kT

    flushed_print("\nStarting simulation setup.")
    stars, time = setting_up.initialize_stars(params, CONSTS)
    if params.reverse:
        stars = setting_up.reverse_velocities(stars)
        gravity = Gravity(Brutus, CONSTS, time, stars)
    else:
        gravity = Gravity(ph4, CONSTS, time, stars)

    binaries = []
    binding_energies = []
    integration_time = 0

    E0 = gravity.total_energy()
    max_dE = 0

    print("")

    while True:

        if params.reverse:
            reverse_time = params.start_time - (time - params.start_time)
        else:
            reverse_time = None

        if params.reverse:
            print(f"Saving metrics and snapshot at t={reverse_time.number}")
            flushed_print((f"This is {time.number} in integrator time, "
                           f"with start time {params.start_time.number}"))
            metrics = update_metrics(metrics, reverse_time, stars,
                                     gravity.integrator,
                                     binaries, binding_energies, kT,
                                     integration_time)
            file_io.pickle_object(stars, f"snapshot_{reverse_time}.pkl",
                                  params)
        else:
            flushed_print(f"Saving metrics and snapshot at t={time.number}")
            metrics = update_metrics(metrics, time, stars, gravity.integrator,
                                     binaries,
                                     binding_energies, kT, integration_time)
            file_io.pickle_object(stars, f"snapshot_{time}.pkl", params)

        file_io.pickle_object(metrics, "cluster_metrics.pkl", params)

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

        binaries, binding_energies = find_binaries(stars, minimum_Eb)
        if len(binaries) > 0 and params.t_end is None:
            params.t_end = time + (5 | nbody_system.time)
            print(f"Set t_end to {params.t_end}")

    gravity.stop()
    print(f"\nMaximum dE = {max_dE}")

    scatterplot(stars, params.output_folder+"final_state.png")
    write_set_to_file(stars, params.output_folder+"final_state.csv", "csv")
    file_io.round_csv(params.output_folder+"final_state.csv", 3)
