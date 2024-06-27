from keras.models import load_model
import numpy as np
from pkg_resources import resource_filename
from spliceai.utils import one_hot_encode

_paths = ('models/spliceai{}.h5'.format(x) for x in range(1, 6))
models = [load_model(resource_filename('spliceai', x)) for x in _paths]

def predict(seq: str) -> tuple[np.ndarray, np.ndarray]:
    encoded = one_hot_encode(seq)[None, :]
    y = np.mean([models[m].predict(encoded) for m in range(5)], axis=0)
    acceptor_prob = y[0, :, 1]
    donor_prob = y[0, :, 2]
    return (acceptor_prob, donor_prob)
