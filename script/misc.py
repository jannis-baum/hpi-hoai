import os

import matplotlib.pyplot as plt
import numpy as np

from definitions import fig_dir


def savefig_pdf(name):
    os.makedirs(fig_dir, exist_ok=True)
    pdf_path = os.path.join(fig_dir, f'{name}.pdf')

    plt.savefig(pdf_path, bbox_inches='tight')

# root mean squared with additional penalty for under estimation
def asym_rmse(y_true, y_pred, under_estimation_penalty: float):
    diff = y_true - y_pred
    under_estimates = diff > 0

    diff **= 2
    # penalize under estimation
    diff[under_estimates] *= under_estimation_penalty

    # take penalty into account while normalizing to stay in same ball park
    normalized_count = sum(under_estimates) * under_estimation_penalty + sum(~under_estimates)
    return np.sqrt(diff.sum() / normalized_count)
