import os

from script.gene import Gene

root_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(root_dir, 'data')

brca1 = Gene('BRCA1', 'chr17', 43044295, 43170245)
