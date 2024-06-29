from keras.models import load_model
import numpy as np
from pkg_resources import resource_filename
from spliceai.utils import one_hot_encode

from script.gene import Gene
from script.sequences import Annotator

class SpliceSitePredictor:
    def __init__(self):
        raise NotImplementedError('Abstract class cannot be instantiated')

    # probabilities donor, acceptor, index mapping to reference genome
    def predict(self, annotator: Annotator, gene: Gene) -> tuple[np.ndarray, np.ndarray, list[list[int]]]:
        raise NotImplementedError('Method needs to be implemented by subclass')

class _SpliceAI(SpliceSitePredictor):
    _lost_padding = 5000

    def __init__(self):
        self._models = None
        self.predict = self._predict

    def _get_models(self) -> list:
        if self._models is None:
            _paths = ('models/spliceai{}.h5'.format(x) for x in range(1, 6))
            self._models = [load_model(resource_filename('spliceai', x)) for x in _paths]
        return self._models

    def _predict(self, annotator: Annotator, gene: Gene) -> tuple[np.ndarray, np.ndarray, list[list[int]]]:
        seq, var_idx = annotator.get_gene_seq(gene, padding=self._lost_padding)
        encoded = one_hot_encode(seq)[None, :]

        y = np.mean([model.predict(encoded) for model in self._get_models()], axis=0)
        donor_prob = y[0, :, 2]
        acceptor_prob = y[0, :, 1]

        return (donor_prob, acceptor_prob, var_idx[self._lost_padding:-self._lost_padding])

spliceai = _SpliceAI()
