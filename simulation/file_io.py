import os
import pickle


def create_directory(directory_name):
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)


def remove_file(filename):
    if os.path.exists(filename):
        os.remove(filename)


def pickle_object(object, filename, params):
    object_path = params.output_folder+filename
    pickle.dump(object, open(object_path, "wb"))
