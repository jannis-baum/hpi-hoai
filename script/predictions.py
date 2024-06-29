from keras.models import load_model
import numpy as np
from pkg_resources import resource_filename
from spliceai.utils import one_hot_encode
import torch
from pangolin import model as _pangolin_model

from script.gene import Gene
from script.sequences import Annotator

class SpliceSitePredictor:
    def __init__(self):
        raise NotImplementedError('Abstract class cannot be instantiated')

    # probabilities donor, acceptor, index mapping to reference genome
    def predict(self, annotator: Annotator, gene: Gene) -> tuple[np.ndarray, list[list[int]]]:
        raise NotImplementedError('Method needs to be implemented by subclass')

# ------------------------------------------------------------------------------
# MARK: SpliceAI

# adapted from https://github.com/Illumina/SpliceAI/tree/master#examples
class _SpliceAI(SpliceSitePredictor):
    _lost_padding = 5000
    # SpliceAI has different interpretation of which base pair exactly is the
    # "donor" site
    _d_offset = 2

    def __init__(self):
        self._models = None
        self.predict = self._predict

    def _get_models(self) -> list:
        if self._models is None:
            _paths = ('models/spliceai{}.h5'.format(x) for x in range(1, 6))
            self._models = [load_model(resource_filename('spliceai', x)) for x in _paths]
        return self._models

    def _predict(self, annotator: Annotator, gene: Gene) -> tuple[np.ndarray, list[list[int]]]:
        seq, var_idx = annotator.get_gene_seq(gene, padding=self._lost_padding)
        encoded = one_hot_encode(seq)[None, :]

        y = np.mean([model.predict(encoded) for model in self._get_models()], axis=0)
        donor_prob = y[0, :, 2]
        acceptor_prob = y[0, :, 1]
        # site can only be donor OR acceptor -> we can take the maximum
        splice_prob = np.maximum.reduce([
            np.concatenate([np.zeros(self._d_offset), donor_prob]),
            np.concatenate([acceptor_prob, np.zeros(self._d_offset)])
       ])[:-self._d_offset]

        return (splice_prob, var_idx[self._lost_padding:-self._lost_padding])

spliceai = _SpliceAI()

# ------------------------------------------------------------------------------
# MARK: Pangolin

# adapted from https://github.com/tkzeng/Pangolin/blob/main/scripts/custom_usage.py
class _Pangolin(SpliceSitePredictor):
    # 0 = Heart, P(splice)
    # 1 = Heart, usage
    # 2 = Liver, P(splice)
    # 3 = Liver, usage
    # 4 = Brain, P(splice)
    # 5 = Brain, usage
    # 6 = Testis, P(splice)
    # 7 = Testis, usage
    _model_num = 6
    _index_map = {0:1, 1:2, 2:4, 3:5, 4:7, 5:8, 6:10, 7:11}

    _lost_padding = 5000

    def __init__(self):
        # use gpu > apple silicon > cpu
        self._device = torch.device('cuda:0' if torch.cuda.is_available() else (torch.device('mps') if torch.backends.mps.is_available() else 'cpu')) 
        self._models = None

        self.predict = self._predict

    def _get_models(self) -> list:
        if self._models is None:
            self._models = list()
            for j in range(1, 6):
                model = _pangolin_model.Pangolin(_pangolin_model.L, _pangolin_model.W, _pangolin_model.AR)
                model.to(self._device)
                weights = torch.load(resource_filename('pangolin', f'models/final.{j}.{self._model_num}.3'), map_location=self._device)
                model.load_state_dict(weights)
                model.eval()
                self._models.append(model)

        return self._models

    def _predict(self, annotator: Annotator, gene: Gene) -> tuple[np.ndarray, list[list[int]]]:
        seq, var_idx = annotator.get_gene_seq(gene)
        encoded = one_hot_encode(seq).T

        x = torch.from_numpy(np.expand_dims(encoded, axis=0)).float().to(self._device)
        score = []
        # Average across 5 models
        for model in self._get_models():
            with torch.no_grad():
                score.append(model(x)[0][self._index_map[self._model_num],:].cpu().numpy())
        y = np.mean(score, axis=0)

        return (y, var_idx[self._lost_padding:-self._lost_padding])

pangolin = _Pangolin()
