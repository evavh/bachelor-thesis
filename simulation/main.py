import sys
import numpy
import getopt
import os

from amuse.ic.plummer import new_plummer_model
from amuse.ic.salpeter import new_salpeter_mass_distribution_nbody
from amuse.units import nbody_system
from amuse.units import units
from amuse.community.ph4.interface import ph4 as grav
from amuse.community.kepler.interface import Kepler
from amuse.ext.LagrangianRadii import LagrangianRadii

from amuse.io import write_set_to_file
from amuse import datamodel
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


def print_log(s, gravity, E0=0.0 | nbody_system.energy):
    M = gravity.total_mass
    U = gravity.potential_energy
    T = gravity.kinetic_energy
    try:
        Ebin = gravity.get_binary_energy()
    except:
        Ebin = 0 | nbody_system.energy
    Etop = T + U
    E = Etop + Ebin
    if E0 == 0 | nbody_system.energy:
        E0 = E
    print("")
    print("time = ", gravity.get_time(), "m=", M, "de=", (E-E0)/E0)
    print("Energies = ", E, U, T)

    return E


def get_component_binary_elements(comp1, comp2):
    kep = Kepler(redirection="none")
    kep.initialize_code()

    mass = comp1.mass + comp2.mass
    pos = comp2.position - comp1.position
    vel = comp2.velocity - comp1.velocity
    kep.initialize_from_dyn(mass, pos[0], pos[1], pos[2],
                            vel[0], vel[1], vel[2])
    a, e = kep.get_elements()
    r = kep.get_separation()
    E, J = kep.get_integrals()  # per unit reduced mass, note
    kep.stop()

    return mass, a, e, r, E


def find_binaries(particles,
                  minimal_binding_energy=-1.e-4,
                  G=nbody_system.G):

    # Search for and print out bound pairs using a numpy-accelerated
    # N^2 search.

    binding_energies = []
    binaries = datamodel.Particles(0)
    for p in particles:
        mu = p.mass*particles.mass/(p.mass+particles.mass)
        dr = (particles.position - p.position).lengths()
        dv = (particles.velocity - p.velocity).lengths()
        E = 0.5*mu*dv*dv - G*p.mass*particles.mass/dr
        indices = numpy.argsort(E.number)
        sorted_E = E[indices]
        Emin = sorted_E[1].number
        if Emin < minimal_binding_energy and p.id < particles[indices[1]].id:
            print('bound', p.id, particles[indices[1]].id, Emin)
            binding_energies.append(Emin)
            binaries.add_particle(p)
            binaries.add_particle(particles[indices[1]])
    return binaries, binding_energies


def test_multiples(infile=None, number_of_stars=40,
                   end_time=10 | nbody_system.time,
                   delta_t=1 | nbody_system.time,
                   n_workers=1, use_gpu=1, use_gpu_code=1,
                   accuracy_parameter=0.1,
                   softening_length=-1 | nbody_system.length):

    if infile is not None:
        print("input file =", infile)
    print("end_time =", end_time)
    print("delta_t =", delta_t)
    print("n_workers =", n_workers)
    print("use_gpu =", use_gpu)
    print("\ninitializing the gravity module")

    metrics = {"times": [] | nbody_system.time,
               "rvir": [] | nbody_system.length,
               "rcore": [] | nbody_system.length,
               "r10pc": [] | nbody_system.length,
               "r50pc": [] | nbody_system.length,
               "density_centre": [],
               "core_density": []}

    # -----------------------------------------------------------------

    if softening_length == -1 | nbody_system.length:
        eps2 = 0.25*(float(number_of_stars))**(-0.666667) \
            | nbody_system.length**2
    else:
        eps2 = softening_length*softening_length

    # -----------------------------------------------------------------

    print("making a Plummer model")
    stars = new_plummer_model(number_of_stars)

    id = numpy.arange(number_of_stars)
    stars.id = id+1 | units.none

    print("setting particle masses and radii")
    if ADD_MASS_FUNCTION:
        scaled_mass = new_salpeter_mass_distribution_nbody(number_of_stars)
        stars.mass = scaled_mass
    stars.radius = 0.0 | nbody_system.length

    print("centering stars")
    stars.move_to_center()
    print("scaling stars to virial equilibrium")
    stars.scale_to_standard(smoothing_length_squared=eps2)

    time = 0.0 | nbody_system.time

    # -----------------------------------------------------------------

    # Note that there are actually three GPU options to test:
    #
    # 1. use the GPU code and allow GPU use (default)
    # 2. use the GPU code but disable GPU use (-g)
    # 3. use the non-GPU code (-G)

    if use_gpu_code == 1:
        try:
            gravity = grav(number_of_workers=n_workers,
                           redirection="none", mode="gpu")
        except:
            gravity = grav(number_of_workers=n_workers, redirection="none")
    else:
        gravity = grav(number_of_workers=n_workers, redirection="none")

    gravity.initialize_code()
    print("Initialized code")
    gravity.parameters.set_defaults()
    print("Set defaults")

    gravity.parameters.timestep_parameter = accuracy_parameter
    gravity.parameters.epsilon_squared = eps2
    gravity.parameters.use_gpu = use_gpu

    print("adding particles")
    gravity.particles.add_particles(stars)
    gravity.commit_particles()

    print('')
    print("number_of_stars =", number_of_stars)
    print("evolving to time =", end_time,
          "in steps of", delta_t)

    E0 = print_log('ph4', gravity)
    dEmult = 0.0

    # Channel to copy values from the code to the set in memory.
    channel = gravity.particles.new_channel_to(stars)

    stopping_condition = gravity.stopping_conditions.collision_detection
    stopping_condition.enable()

    zero = 0 | nbody_system.time
    while True:
        time += delta_t

        if end_time > zero and time >= end_time:
            break

        while gravity.get_time() < time:
            gravity.evolve_model(time)
            if stopping_condition.is_set():
                print('Collision detected at time',
                      gravity.get_time())

                E = print_log('ph4', gravity, E0)
                print('dEmult =', dEmult, 'dE =', (E-E0)-dEmult)
                break

        # Copy values from the module to the set in memory.
        channel.copy()
        channel.copy_attribute("index_in_code", "id")

        binaries, Eb = find_binaries(stars, minimal_binding_energy=-0.0001)
        if len(binaries) > 0:
            print("Binding energies:", Eb)
            mass, a, e, r, E = get_component_binary_elements(
                binaries[0], binaries[1])
            print("Binary:", time, mass, a, e, r, E)
            if numpy.min(Eb) < -0.1 and end_time < zero:
                end_time = time + (20 | nbody_system.time)

        Rl = LagrangianRadii(stars, massf=[0.1, 0.5, 1.0] | units.none)
        metrics["r10pc"].append(Rl[0])
        metrics["r50pc"].append(Rl[1])

        density_centre, core_radius, core_density = \
            gravity.particles.densitycentre_coreradius_coredens()
        metrics["rvir"].append(stars.virial_radius())
        metrics["rcore"].append(core_radius)
        metrics["times"].append(time)
        metrics["density_centre"].append(density_centre)
        metrics["core_density"].append(core_density)

        write_set_to_file(stars.savepoint(time),
                          "simulation/output/snapshots.hdf5", "hdf5",
                          append_to_file=True)

        E = print_log('ph4', gravity, E0)
        print('dEmult =', dEmult, 'dE =', (E-E0))

    print('')
    gravity.stop()

    scatterplot(stars, "simulation/output/final_state.png")

    write_set_to_file(stars, "simulation/output/final_state.csv", "csv")

    pyplot.plot(metrics["times"].value_in(nbody_system.time),
                metrics["rcore"].value_in(nbody_system.length), c='b')
    pyplot.plot(metrics["times"].value_in(nbody_system.time),
                metrics["rvir"].value_in(nbody_system.length), c='r')
    pyplot.plot(metrics["times"].value_in(nbody_system.time),
                metrics["r10pc"].value_in(nbody_system.length), c='k', ls="--", lw=1)
    pyplot.plot(metrics["times"].value_in(nbody_system.time),
                metrics["r50pc"].value_in(nbody_system.length), c='k', ls="--", lw=3)
    pyplot.xlabel("time")
    pyplot.ylabel("radius")
    pyplot.semilogy()
    pyplot.savefig("simulation/output/radii.png")


if __name__ == '__main__':

    infile = None
    N = 100
    t_end = -1 | nbody_system.time
    delta_t = 1.0 | nbody_system.time
    n_workers = 1
    use_gpu = 0
    use_gpu_code = 0
    accuracy_parameter = 0.1
    softening_length = 0 | nbody_system.length
    random_seed = -1

    try:
        opts, args = getopt.getopt(sys.argv[1:], "a:c:d:e:f:gGn:s:t:w:")
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
        elif o == "-f":
            infile = a
        elif o == "-g":
            use_gpu = 1
        elif o == "-G":
            use_gpu = 1
            use_gpu_code = 1
        elif o == "-n":
            N = int(a)
        elif o == "-s":
            random_seed = int(a)
        elif o == "-t":
            t_end = float(a) | nbody_system.time
        elif o == "-w":
            n_workers = int(a)
        else:
            print("unexpected argument", o)

    create_directory("simulation/output")
    remove_file("simulation/output/snapshots.hdf5")

    if random_seed <= 0:
        numpy.random.seed()
        random_seed = numpy.random.randint(1, pow(2, 31)-1)
    numpy.random.seed(random_seed)
    print("random seed =", random_seed)

    assert is_mpd_running()
    test_multiples(infile, N, t_end, delta_t, n_workers,
                   use_gpu, use_gpu_code,
                   accuracy_parameter, softening_length)
