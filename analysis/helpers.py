import numpy


def ids_to_stars(snapshot, ids):
    def criterium(id): return id in ids
    return snapshot.select(criterium, ['id'])


def stars_to_ids(stars):
    ids = []
    for star in stars:
        ids.append(star.id)
    return ids


def time_to_index(time, data):
    metrics = data.metrics

    index = None
    for i, t in zip(range(len(metrics['times'])), metrics['times']):
        if t.number >= time and index is None:
            index = i

    return index


def limits_from_radius(radius, density_centre):
    x_min = density_centre[0] - radius
    x_max = density_centre[0] + radius

    y_min = density_centre[1] - radius
    y_max = density_centre[1] + radius

    z_min = density_centre[2] - radius
    z_max = density_centre[2] + radius

    return ((x_min, x_max), (y_min, y_max), (z_min, z_max))


def stars_in_area(snapshot, density_centre, radius):
    limits = limits_from_radius(radius, density_centre)

    def within_limits(x, y, z):
        x_lims, y_lims, z_lims = limits
        x_min, x_max = x_lims
        y_min, y_max = y_lims
        z_min, z_max = z_lims

        return (x > x_min and x < x_max and
                y > y_min and y < y_max and
                z > z_min and z < z_max)

    return snapshot.select(within_limits, ['x', 'y', 'z'])


def work_slice(work, start, end):
    assert start <= end

    if start == 0:
        return work[start:end]

    return work[start:end] - work[start-1]


def slice_works(star_works, start_index, end_index):
    cropped_star_works = {}
    cropped_work_totals = {}

    for star in star_works:
        cropped_star_works[star] = work_slice(star_works[star],
                                              start_index,
                                              end_index)
        cropped_work_totals[star] = numpy.sum(cropped_star_works[star])

    return cropped_star_works, cropped_work_totals
