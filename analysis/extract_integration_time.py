import argparse
import re
from matplotlib import pyplot

import input_output
import plotting


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="stdout file to be analysed",
                        type=str)
    parser.add_argument("-o", "--output", help="where to put the plot",
                        type=str)
    parser.add_argument("-d", "--dir", help="directory where snapshots are",
                        type=str)

    arguments = parser.parse_args()
    if not arguments.dir.endswith('/'):
        arguments.dir += '/'
    if not arguments.output.endswith('/'):
        arguments.output += '/'

    return arguments


if __name__ == '__main__':
    config = parse_arguments()
    with open(config.input, 'r') as file:
        stdout = file.read()

    integration_times = []
    matches = re.findall(r'(.*)s elapsed.', stdout)
    for match in matches:
        integration_times.append(float(match))

    times = input_output.extract_times(config.dir).number

    pyplot.plot(times[1:], integration_times[:225], color='k')
    pyplot.xlabel("t [n-body time]")
    pyplot.ylabel("Integration time [s]")

    pyplot.savefig(config.output+"integration_time.svg", format='svg')
    pyplot.clf()
