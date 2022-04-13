import math

from amuse.community.brutus.interface import Brutus
from amuse.community.ph4.interface import ph4
from amuse.units import nbody_system

import setting_up


class Gravity:
    def __init__(self, integrator, CONSTS, time, stars):
        if integrator not in [Brutus, ph4]:
            print(f"Invalid or unimplemented integrator: {integrator}")

        self.integrator = integrator(number_of_workers=1, redirection="none")
        self.integrator.initialize_code()

        self.integrator.parameters.set_defaults()
        self.integrator.parameters.begin_time = time

        if integrator == Brutus:
            print("Using Brutus")

            self.integrator.parameters.bs_tolerance = CONSTS['bs_tolerance']
            word_length = 4*abs(math.log(CONSTS['bs_tolerance'], 10)) + 32
            self.integrator.parameters.word_length = word_length

        elif integrator == ph4:
            print("Using ph4")

            self.integrator.parameters.timestep_parameter = CONSTS['accuracy']
            self.integrator.parameters.initial_timestep_fac =\
                CONSTS['initial_accuracy']
            self.integrator.parameters.epsilon_squared = \
                CONSTS['epsilon_squared']
            self.integrator.parameters.use_gpu = 0

        self.integrator.commit_parameters()

        self.integrator.particles.add_particles(stars)
        self.integrator.commit_particles()
        self.channel = self.integrator.particles.new_channel_to(stars)


if __name__ == '__main__':
    CONSTS = {'accuracy': 0.01,
              'initial_accuracy': 0.0025,
              'epsilon_squared': 0 | nbody_system.length**2,
              'bs_tolerance': 1e-20}
    stars = setting_up.new_cluster_model(100, CONSTS['epsilon_squared'])
    time = 0.0 | nbody_system.time

    grav = Gravity(Brutus, CONSTS, time, stars)
