import sys
import numpy
import getopt
import os
import pickle

from amuse.ic.plummer import new_plummer_model
from amuse.ic.salpeter import new_salpeter_mass_distribution_nbody
from amuse.units import nbody_system
from amuse.units import units
from amuse.community.ph4.interface import ph4 as grav
from amuse.ext.LagrangianRadii import LagrangianRadii

from amuse.io import write_set_to_file, read_set_from_file
from amuse.rfi.core import is_mpd_running

from matplotlib import pyplot

ADD_MASS_FUNCTION = False


def create_directory(directory_name):
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)


def remove_file(filename):
    if os.path.exists(filename):
        os.remove(filename)


def scatterplot(stars, filename):
    pyplot.scatter(stars.x.value_in(nbody_system.length),
                   stars.y.value_in(nbody_system.length),
                   s=10*stars.mass/stars.mass.max())
    pyplot.xlabel("x")
    pyplot.ylabel("y")
    pyplot.savefig(filename)
    pyplot.clf()


def new_cluster_model(N, eps2):
    print("making a Plummer model")
    stars = new_plummer_model(N)

    id = numpy.arange(N)
    stars.id = id+1 | units.none

    print("setting particle masses and radii")
    if ADD_MASS_FUNCTION:
        scaled_mass = new_salpeter_mass_distribution_nbody(N)
        stars.mass = scaled_mass
    stars.radius = 0.0 | nbody_system.length

    print("centering stars")
    stars.move_to_center()
    print("scaling stars to virial equilibrium")
    stars.scale_to_standard(smoothing_length_squared=eps2)

    return stars


def find_snapshot(snapshots, start_time):
    for snapshot in snapshots.history:
        time = snapshot.get_timestamp()
        if time >= start_time:
            return snapshot, time
    raise Exception("Start time not found in snapshot file.")


if __name__ == '__main__':
    N = 100
    t_end = -1 | nbody_system.time
    start_time = -1 | nbody_system.time
    delta_t = 1.0 | nbody_system.time
    accuracy_parameter = 0.1
    softening_length = 0 | nbody_system.length
    random_seed = None
    output_folder = "simulation/output"
    snapshot_input = output_folder
    minimum_Eb_kT = 10

    try:
        opts, args = getopt.getopt(sys.argv[1:],
                                   "a:d:e:n:s:t:T:o:i:b:")
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(1)

    for o, a in opts:
        if o == "-a":
            accuracy_parameter = float(a)
        elif o == "-d":
            delta_t = float(a) | nbody_system.time
        elif o == "-e":
            softening_length = float(a) | nbody_system.length
        elif o == "-n":
            N = int(a)
        elif o == "-s":
            random_seed = int(a)
        elif o == "-t":
            t_end = float(a) | nbody_system.time
        elif o == "-T":
            start_time = float(a) | nbody_system.time
        elif o == "-o":
            output_folder = a
        elif o == "-i":
            snapshot_input = a
        elif o == "-b":
            minimum_Eb_kT = float(a)
        else:
            print("unexpected argument", o)

    create_directory(output_folder)
    remove_file(output_folder+"/snapshots.hdf5")

    def set_random_seed(random_seed):
        if random_seed is None:
            numpy.random.seed()
            random_seed = numpy.random.randint(1, pow(2, 31)-1)
        numpy.random.seed(random_seed)
        return random_seed

    random_seed = set_random_seed(random_seed)
    print("random seed =", random_seed)

    assert is_mpd_running()

    metrics = {"N": N,
               "times": [] | nbody_system.time,
               "rvir": [] | nbody_system.length,
               "rcore": [] | nbody_system.length,
               "r10pc": [] | nbody_system.length,
               "r50pc": [] | nbody_system.length,
               "density_centre": [],
               "core_density": [],
               "total_mass": [],
               "potential_energy": [],
               "kinetic_energy": [],
               "total_binary_energy": []}

    binaries_found = False

    if softening_length == -1 | nbody_system.length:
        eps2 = 0.25*(float(N))**(-0.666667) \
            | nbody_system.length**2
    else:
        eps2 = softening_length*softening_length

    if start_time == -1 | nbody_system.time:
        stars = new_cluster_model(N, eps2)
        time = 0.0 | nbody_system.time
    else:
        snapshots = read_set_from_file(snapshot_input, 'hdf5')
        stars, time = find_snapshot(snapshots, start_time)

    write_set_to_file(stars.savepoint(time),
                      output_folder+"/snapshots.hdf5", "hdf5",
                      append_to_file=True)

    gravity = grav(number_of_workers=1, redirection="none")

    gravity.initialize_code()
    print("Initialized code")
    gravity.parameters.set_defaults()
    print("Set defaults")

    gravity.parameters.timestep_parameter = accuracy_parameter
    gravity.parameters.epsilon_squared = eps2
    gravity.parameters.use_gpu = 0

    print("adding particles")
    gravity.particles.add_particles(stars)
    gravity.commit_particles()

    print('')
    print("number_of_stars =", N)
    print("evolving to time =", t_end,
          "in steps of", delta_t)

    # Channel to copy values from the code to the set in memory.
    channel = gravity.particles.new_channel_to(stars)

    stopping_condition = gravity.stopping_conditions.collision_detection
    stopping_condition.enable()

    kT = 1/(6*N)
    minimum_Eb = minimum_Eb_kT * kT
    print("kT =", kT)
    print("minimum Eb =", minimum_Eb, "=", minimum_Eb_kT, "kT")

    zero = 0 | nbody_system.time
    while True:
        print("Starting integration at time", time)
        time += delta_t

        if t_end > zero and time >= t_end:
            break

        while gravity.get_time() < time:
            gravity.evolve_model(time)
            if stopping_condition.is_set():
                print('Collision detected at time',
                      gravity.get_time())

                break

        # Copy values from the module to the set in memory.
        channel.copy()
        channel.copy_attribute("index_in_code", "id")

        print("Finished integration, starting binary finding.")

        # Search for and print out bound pairs using a numpy-accelerated
        # N^2 search.
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

        if len(binaries) > 0:
            print("Binding energies:", binding_energies)
            binding_energies_kT = [x/kT for x in binding_energies]
            print("Binding energies in kT:", binding_energies_kT)
            if not binaries_found:
                binaries_found = True
                metrics["first_binaries"] = binaries
                metrics["first_binary_energies_kT"] = binding_energies_kT
                metrics["first_binary_time"] = time
            if t_end < zero:
                t_end = time + (20 | nbody_system.time)

        print("Finished binary finding, starting filling metrics struct.")

        metrics["times"].append(time)

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

        print("Finished filling metrics struct, starting writing snapshot.")

        write_set_to_file(stars.savepoint(time),
                          output_folder+"/snapshots.hdf5", "hdf5",
                          append_to_file=True)

        print("Finished writing snapshot, starting pickling metrics.")

        metrics_filename = output_folder+"/cluster_metrics.pkl"
        pickle.dump(metrics, open(metrics_filename, "wb"))

        print("Finished writing data to files, going to next loop")

    print('')
    gravity.stop()

    scatterplot(stars, output_folder+"/final_state.png")

    write_set_to_file(stars, output_folder+"/final_state.csv", "csv")
