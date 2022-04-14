import math
from amuse.units import nbody_system


def crc(core_density):
    G = nbody_system.G
    under_sqrt = 4*math.pi*G*core_density/3
    return 1/under_sqrt.sqrt()
