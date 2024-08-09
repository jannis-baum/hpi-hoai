from hashlib import md5
import os
import pickle
from typing import Callable, TypeVar

import pandas as pd

from definitions import data_dir

_filename_to_path = dict[str, str]()

def index():
    for (root, _, files) in os.walk(data_dir):
        for f in files:
            _filename_to_path[f] = os.path.join(root, f)
index()

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

def _hash(args):
    return str(int.from_bytes(
        md5(str(args).encode()).digest(),
        'big'
    ))

class Hashable():
    identity = []
    def hash(self):
        return _hash(self.identity)

T = TypeVar('T')
def retrieve_or_compute(compute_fn: Callable[[], T], *identifiers) -> T:
    filename = _hash([i.hash() if isinstance(i, Hashable) else i for i in identifiers])
    try:
        with open(find_path(filename), 'rb') as fp:
            result = pickle.load(fp)
            return result
    except:
        path = make_path('generic-cache', filename)
        result = compute_fn()
        with open(path, 'wb') as fp:
            pickle.dump(result, fp)
        return result
