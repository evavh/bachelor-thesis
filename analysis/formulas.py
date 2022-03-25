import math
import numpy

from amuse.units import nbody_system
from amuse.datamodel import Particles


def t_rh(N, r_h, G, M):
    return 0.138*N*r_h**(3/2)/(G**(1/2)*M**(1/2)*math.log(0.4*N))


def binding_energy(star1, star2):
    mu = star1.mass*star2.mass/(star1.mass+star2.mass)
    dr = (star2.position - star1.position).lengths()
    dv = (star2.velocity - star1.velocity).lengths()
    Eb = nbody_system.G*star1.mass*star2.mass/dr - 0.5*mu*dv*dv

    return Eb


def norm(vector):
    return numpy.sqrt(vector.dot(vector))


def f_ik_j(component_k, star_j):
    G = nbody_system.G
    m_j = star_j.mass
    m_ik = component_k.mass
    r_j = numpy.array([star_j.x.number, star_j.y.number, star_j.z.number])\
        | star_j.x.unit
    r_ik = numpy.array([component_k.x.number, component_k.y.number,
                        component_k.z.number]) | star_j.x.unit

    force_vector = -G*m_j*m_ik*(r_ik-r_j)/(norm(r_ik-r_j)**3)
    assert (force_vector.unit == nbody_system.mass * nbody_system.acceleration)

    return force_vector


def power_function(tuple, star):
    dEdt = 0 | nbody_system.energy / nbody_system.time
    for k in tuple:
        v_k = numpy.array([k.vx.number, k.vy.number, k.vz.number]) | k.vx.unit
        v_cm = tuple.center_of_mass_velocity()
        dEdt -= f_ik_j(k, star).dot(v_k - v_cm)

    assert (dEdt.unit == nbody_system.energy / nbody_system.time)

    return dEdt


def work_function(snapshots, metrics, tuple_ids, star_id,
                  start_time, end_time):
    start_index = None
    for index, snapshot, time in zip(range(len(snapshots)), snapshots,
                                     metrics['times']):
        if time >= start_time and start_index is None:
            start_index = index

        if time >= end_time:
            end_index = index
            break

    power_functions = []
    for snapshot in snapshots[start_index:end_index]:
        tuple = Particles(0)
        for particle in snapshot:
            if particle.id in tuple_ids:
                tuple.add_particle(particle)
            if particle.id == star_id:
                star = particle

        power = power_function(tuple, star) * (1.0 | nbody_system.time)
        power_functions.append(power.value_in(nbody_system.energy))

    work = numpy.cumsum(power_functions)
    total_work = numpy.sum(work)

    return work, total_work
