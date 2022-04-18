import os
import argparse


def create_directory(directory_name):
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Data directory to be analysed",
                        default="analysis/data/")
    parser.add_argument("-o", "--output", help="Directory to put output",
                        default="analysis/plots/")
    parser.add_argument("--scatter", help="generate scatter plots.",
                        action='store_true')
    parser.add_argument("--scatter_core", help="generate scatter of core.",
                        action='store_true')
    parser.add_argument("--load_work", help="load work data from pickle.",
                        action='store_true')
    parser.add_argument("--calc_work", help="calculate work (pricey).",
                        action='store_true')
    parser.add_argument("-og", "--original_run", help="run from t0",
                        type=str)
    arguments = parser.parse_args()
    if not arguments.input.endswith("/"):
        arguments.input += "/"
    if not arguments.output.endswith("/"):
        arguments.output += "/"
    if arguments.original_run is not None \
            and not arguments.original_run.endswith("/"):
        arguments.original_run += "/"

    return arguments
