import os

import matplotlib.pyplot as plt

from definitions import fig_dir


def savefig_pdf(name):
    os.makedirs(fig_dir, exist_ok=True)
    pdf_path = os.path.join(fig_dir, f'{name}.pdf')

    plt.savefig(pdf_path, bbox_inches='tight')
