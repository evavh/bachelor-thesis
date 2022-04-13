import os
import numpy
import collections

from amuse.ic.plummer import new_plummer_model
from amuse.units import nbody_system
from amuse.units import units
from amuse.community.ph4.interface import ph4 as grav

import file_io


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


def find_snapshot(params):
    start_time = params.start_time
    for filename in os.listdir(params.snapshot_input):
        if filename.endswith('.pkl') and str(start_time) in filename:
            stars = file_io.unpickle_object(filename, params)
            return stars, start_time
    raise Exception("Start time not found in snapshot folder.")


def initialize_stars(params, CONSTS):
    if params.start_time is None:
        stars = new_cluster_model(params.n, CONSTS['epsilon_squared'])
        time = 0.0 | nbody_system.time
    else:
        stars, time = find_snapshot(params)

    return stars, time


def reverse_velocities(stars):
    stars.vx *= -1
    stars.vy *= -1
    stars.vz *= -1
    return stars
