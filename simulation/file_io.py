import os
import pickle
import pandas


def create_directory(directory_name):
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)


def remove_file(filename):
    if os.path.exists(filename):
        os.remove(filename)


def pickle_object(object, filename, params):
    object_path = params.output_folder+filename
    pickle.dump(object, open(object_path, "wb"))


def unpickle_object(filename, params):
    object_path = params.snapshot_input+'/'+filename
    return pickle.load(open(object_path, 'rb'))


def round_csv(filename, decimals):
    data = pandas.read_csv(filename)
    data = data.iloc[3:]
    data = data.astype(float).round(3)
    data.to_csv(filename[:-4]+"_rounded.csv", index=False)
