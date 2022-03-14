import file_io
import setting_up

import math
import numpy
import argparse
import sys
import datetime

from amuse.units import nbody_system
from amuse.units import units
from amuse.ext.LagrangianRadii import LagrangianRadii

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

    parameters = parser.parse_args()

    parameters.delta_t = parameters.delta_t | nbody_system.time
    if parameters.t_end is not None:
        parameters.t_end = parameters.t_end | nbody_system.time
    if parameters.start_time is not None:
        parameters.start_time = parameters.start_time | nbody_system.time

    if not parameters.output_folder.endswith('/'):
        parameters.output_folder += '/'

    return parameters


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
    metrics["total_binary_energy"].append(gravity.get_binary_energy())

    metrics['t_crc'].append(calculate_t_crc(core_density))

    metrics['integration_time'].append(integration_time)

    return metrics


if __name__ == '__main__':
    flushed_print("Simulation script started.")

    assert is_mpd_running()

    metrics = setting_up.initialize_metrics()
    CONSTS = {'accuracy': 0.01,
              'initial_accuracy': 0.0025,
              'epsilon_squared': 0 | nbody_system.length**2}
    params = parse_arguments()
    params.random_seed = set_random_seed(params.random_seed)

    file_io.pickle_object(params, "parameters.pkl", params)
    file_io.pickle_object(CONSTS, "constants.pkl", params)

    kT = 1/(6*params.n)
    minimum_Eb = params.minimum_Eb_kT * kT
    file_io.create_directory(params.output_folder)

    flushed_print("Starting simulation setup.")
    stars, time = setting_up.initialize_stars(params, CONSTS)
    gravity = setting_up.setup_integrator(stars, CONSTS)
    gravity.parameters.begin_time = time

    channel = gravity.particles.new_channel_to(stars)
    stopping_condition = gravity.stopping_conditions.collision_detection
    stopping_condition.enable()

    binaries = []
    binding_energies = []
    binary_queue = []
    first_binary = None
    new_first_needed = False

    integration_time = 0

    while True:
        flushed_print(
            f"First star vx at t={time.number}: {stars[0].vx.number}")
        print(f"Saving metrics and snapshot at t={time.number}")
        metrics = update_metrics(metrics, time, stars, gravity, binaries,
                                 binding_energies, kT, integration_time)
        file_io.pickle_object(metrics, "cluster_metrics.pkl", params)

        print(f"First star at t={time.number}: {stars[0]}")
        file_io.pickle_object(stars, f"snapshot_{time}.pkl", params)

        if params.t_end is not None and time >= params.t_end:
            break

        print(f"Starting integration at t={time.number}")
        integration_start_time = datetime.datetime.now()

        t_crc = metrics['t_crc'][-1]
        if params.variable_delta:
            time += 0.01*t_crc
            print(f"delta_t is 0.01*t_crc = {0.01*t_crc} for this iteration")
        else:
            time += params.delta_t
            print(f"delta_t is {params.delta_t} for this iteration")

        while gravity.get_time() < time:
            gravity.evolve_model(time)
            if stopping_condition.is_set():
                print(f"Collision detected at t={gravity.get_time().number}")
                break
        channel.copy()

        print(f"Finished integrating until t={time.number}")
        integration_time = datetime.datetime.now() - integration_start_time
        integration_time = integration_time.total_seconds()
        print(f"{integration_time}s elapsed.")

        print("Starting binary finding.")

        binaries, binding_energies = find_binaries(stars, minimum_Eb)
        if len(binaries) > 0:
            for binary in binaries:
                binary_queue.append((binary, time))
                print((f"Added {binary[0].id.number, binary[1].id.number} "
                       f"to queue, now {len(binary_queue)} long."))

            if params.t_end is None:
                first_binary = binaries[0]
                params.t_end = time + (20 | nbody_system.time)
                print("Setting t_end")

        for binary, time in binary_queue:
            if binding_energy(binary[0], binary[1]) < 0 | nbody_system.energy:
                binary_queue.remove((binary, time))
                print((f"Removed {binary[0].id.number, binary[1].id.number} "
                       f"from queue, now {len(binary_queue)} long."))

                if binary == first_binary:
                    new_first_needed = True
                    print("This was the first binary!")

        if new_first_needed:
            first_binary = binary_queue[0][0]
            print((f"{first_binary[0].id.number, first_binary[1].id.number} "
                   "is now the new first binary."))
            print(f"End time was {params.t_end}")
            params.t_end = binary_queue[0][1] + (20 | nbody_system.time)
            print(f"New end time: {params.t_end}")

    gravity.stop()

    scatterplot(stars, params.output_folder+"final_state.png")
    write_set_to_file(stars, params.output_folder+"final_state.csv", "csv")
    file_io.round_csv(params.output_folder+"final_state.csv", 3)
