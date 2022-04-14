import formulas

import numpy

from amuse.units import units
from amuse.units import nbody_system
from amuse.ext.LagrangianRadii import LagrangianRadii


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

    metrics['t_crc'].append(formulas.crc(core_density))

    metrics['integration_time'].append(integration_time)

    return metrics
