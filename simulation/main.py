import file_io
import setting_up

import numpy
import argparse

from amuse.units import nbody_system
from amuse.units import units
from amuse.ext.LagrangianRadii import LagrangianRadii

from amuse.io import write_set_to_file
from amuse.rfi.core import is_mpd_running

from matplotlib import pyplot


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--delta_t", help="time between snapshots",
                        default=1.0, type=float)
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
                   binding_energies, kT):
    metrics["times"].append(time)

    metrics["binaries"].append(binaries)
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

    return metrics


if __name__ == '__main__':
    print("Simulation script started.")

    assert is_mpd_running()

    metrics = setting_up.initialize_metrics()
    CONSTS = {'accuracy': 0.1,
              'epsilon_squared': 0 | nbody_system.length**2}
    params = parse_arguments()
    params.random_seed = set_random_seed(params.random_seed)

    file_io.pickle_object(params, "parameters.pkl", params)
    file_io.pickle_object(CONSTS, "constants.pkl", params)

    kT = 1/(6*params.n)
    minimum_Eb = params.minimum_Eb_kT * kT
    file_io.create_directory(params.output_folder)

    print("Starting simulation setup.")
    stars, time = setting_up.initialize_stars(params, CONSTS)
    gravity = setting_up.setup_integrator(stars, CONSTS)
    channel = gravity.particles.new_channel_to(stars)
    stopping_condition = gravity.stopping_conditions.collision_detection
    stopping_condition.enable()

    binaries = []
    binding_energies = []

    while True:
        print(f"First star vx at t={time.number}: {stars[0].vx.number}")
        print(f"Saving metrics and snapshot at t={time.number}")
        metrics = update_metrics(metrics, time, stars, gravity, binaries,
                                 binding_energies, kT)
        file_io.pickle_object(metrics, "cluster_metrics.pkl", params)

        print(f"First star at t={time.number}: {stars[0]}")
        file_io.pickle_object(stars, f"snapshot_{time}.pkl", params)

        if params.t_end is not None and time >= params.t_end:
            break

        print(f"Starting integration at t={time.number}")
        time += params.delta_t

        while gravity.get_time() < time:
            gravity.evolve_model(time)
            if stopping_condition.is_set():
                print(f"Collision detected at t={gravity.get_time().number}")
                break
        channel.copy()

        print(f"Finished integrating until t={time.number}")
        print("Starting binary finding.")

        binaries, binding_energies = find_binaries(stars, minimum_Eb)
        if len(binaries) > 0 and params.t_end is None:
            params.t_end = time + (20 | nbody_system.time)

    gravity.stop()

    scatterplot(stars, params.output_folder+"final_state.png")
    write_set_to_file(stars, params.output_folder+"final_state.csv", "csv")
