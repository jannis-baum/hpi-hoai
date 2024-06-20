import os

from definitions import data_dir

_filename_to_path = dict[str, str]()
for (root, dirs, files) in os.walk(data_dir):
    for f in files:
        _filename_to_path[f] = os.path.join(root, f)

def get_path(filename):
    return _filename_to_path[filename]
