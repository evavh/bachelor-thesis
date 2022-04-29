import math
import numpy

from amuse.units import nbody_system

import helpers


def t_rh(data):
    N = data.params.n
    r_h = data.by_time()[0.0]['r50pc']
    G = nbody_system.G
    M = 1 | nbody_system.mass
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
    assert force_vector.unit == nbody_system.mass * nbody_system.acceleration,\
        (f"force vector has unit {force_vector.unit}, expected "
         f"{nbody_system.mass * nbody_system.acceleration}")

    return force_vector


def power_function(tuple, star):
    dEdt = 0 | nbody_system.energy / nbody_system.time
    for k in tuple:
        v_k = numpy.array([k.vx.number, k.vy.number, k.vz.number]) | k.vx.unit
        v_cm = tuple.center_of_mass_velocity()
        dEdt -= f_ik_j(k, star).dot(v_k - v_cm)

    assert (dEdt.unit == nbody_system.energy / nbody_system.time)

    return dEdt


def work_function(data, tuple_ids, star_id,
                  start_index, end_index):
    snapshots = data.snapshots
    power_functions = []

    for snapshot in snapshots[start_index:end_index]:
        tuple = helpers.ids_to_stars(snapshot, tuple_ids)
        star = helpers.ids_to_stars(snapshot, [star_id])[0]

        power = power_function(tuple, star) * (1.0 | nbody_system.time)
        power_functions.append(power.value_in(nbody_system.energy))

    power_functions = numpy.array(power_functions)

    dt = data.delta_ts()[start_index:end_index]
    work = numpy.cumsum(power_functions*dt)
    total_abs_work = numpy.sum(numpy.abs(power_functions*dt))

    return work, total_abs_work
