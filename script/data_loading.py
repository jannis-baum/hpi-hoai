import os

import pandas as pd

from definitions import data_dir

_filename_to_path = dict[str, str]()
for (root, dirs, files) in os.walk(data_dir):
    for f in files:
        _filename_to_path[f] = os.path.join(root, f)

def find_path(filename):
    return _filename_to_path[filename]

def make_path(*components: str) -> str:
    if len(components) == 0:
        raise ValueError('No path components provided')
    if components[-1] in _filename_to_path:
        raise FileExistsError(f'File {components[-1]} already exists')

    path = os.path.join(data_dir, *components)
    directory = os.path.dirname(path)
    if not os.path.exists(directory):
        os.makedirs(directory)

    _filename_to_path[components[-1]] = path
    return path

study_info = pd.read_csv(find_path('info.csv'), index_col='case_id')
