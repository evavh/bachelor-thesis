from amuse.io import read_set_from_file
from amuse.units import nbody_system
from matplotlib import pyplot
import pickle
import os


def create_directory(directory_name):
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)


def scatterplot(stars, time):
    pyplot.scatter(stars.x.value_in(nbody_system.length),
                   stars.y.value_in(nbody_system.length),
                   s=10*stars.mass/stars.mass.max())
    pyplot.xlabel("x")
    pyplot.ylabel("y")
    pyplot.savefig("analysis/plots/scatter"+str(time)+".png")
    pyplot.clf()


if __name__ == '__main__':
    create_directory("analysis/plots")
    create_directory("analysis/plots/scatter")
    create_directory("analysis/data")

    snapshots = read_set_from_file("analysis/data/snapshots.hdf5", "hdf5")
    metrics = pickle.load(open("analysis/data/cluster_metrics.pkl", "rb"))

    for stars in snapshots.history:
        time = stars.get_timestamp()
        print("t=", time, "length=", len(stars))
        scatterplot(stars, time)
