from amuse.io import read_set_from_file
from amuse.units import nbody_system
from matplotlib import pyplot


def scatterplot(stars, time):
    pyplot.scatter(stars.x.value_in(nbody_system.length),
                   stars.y.value_in(nbody_system.length),
                   s=10*stars.mass/stars.mass.max())
    pyplot.xlabel("x")
    pyplot.ylabel("y")
    pyplot.savefig("analysis/plots/"+str(time)+".png")
    pyplot.clf()


if __name__ == '__main__':
    snapshots = read_set_from_file("analysis/snapshots.hdf5", "hdf5")
    for stars in snapshots.history:
        time = stars.get_timestamp()
        print("t=", time, "length=", len(stars))
        scatterplot(stars, time)
