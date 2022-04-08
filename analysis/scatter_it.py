from amuse.units import nbody_system
from matplotlib import pyplot
import argparse
import pickle

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Snapshot to be plotted",
                    type=str)
parser.add_argument("-o", "--output", help="Filename for plot",
                    type=str)
input = parser.parse_args().input
output = parser.parse_args().output

with open(input, 'rb') as inputfile:
    snapshot = pickle.load(inputfile)

x_of_stars = snapshot.x.value_in(nbody_system.length)
y_of_stars = snapshot.y.value_in(nbody_system.length)

pyplot.scatter(x_of_stars, y_of_stars)
pyplot.xlabel("x")
pyplot.ylabel("y")

axes = pyplot.gca()
axes.set_aspect('equal')
pyplot.savefig(output, format='svg')
pyplot.clf()
