import numpy
import collections

from amuse.ic.plummer import new_plummer_model
from amuse.units import nbody_system
from amuse.units import units
from amuse.community.ph4.interface import ph4 as grav

from amuse.io import read_set_from_file


def initialize_metrics():
    metrics = collections.defaultdict(list)

    metrics['times'] = [] | nbody_system.time
    metrics['rvir'] = [] | nbody_system.length
    metrics['rcore'] = [] | nbody_system.length
    metrics['r10pc'] = [] | nbody_system.length
    metrics['r50pc'] = [] | nbody_system.length

    return metrics


def new_cluster_model(N, eps2):
    stars = new_plummer_model(N)

    stars.id = numpy.arange(N)+1 | units.none

    stars.radius = 0.0 | nbody_system.length

    stars.move_to_center()
    stars.scale_to_standard(smoothing_length_squared=eps2)

    return stars


def find_snapshot(snapshots, start_time):
    for snapshot in snapshots.history:
        time = snapshot.get_timestamp()
        if time >= start_time:
            return snapshot, time
    raise Exception("Start time not found in snapshot file.")


def initialize_stars(params, CONSTS):
    if params.start_time is None:
        stars = new_cluster_model(params.n, CONSTS['epsilon_squared'])
        time = 0.0 | nbody_system.time
    else:
        snapshots = read_set_from_file(params.snapshot_input, 'hdf5')
        stars, time = find_snapshot(snapshots, params.start_time)

    return stars, time


def setup_integrator(stars, CONSTS):
    gravity = grav(number_of_workers=1, redirection="none")

    gravity.initialize_code()
    gravity.parameters.set_defaults()

    gravity.parameters.timestep_parameter = CONSTS['accuracy']
    gravity.parameters.epsilon_squared = CONSTS['epsilon_squared']
    gravity.parameters.use_gpu = 0

    gravity.particles.add_particles(stars)
    gravity.commit_particles()

    return gravity
