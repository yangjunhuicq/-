import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

adata = sc.read('input.h5ad')
sc.tl.embedding_density(adata, groupby='disease', basis='umap')
with rc_context({'figure.figsize': (20, 15)}):
    sc.pl.embedding_density(adata, key='umap_density_disease', basis='umap')
plt.savefig('outpre.DensityPlot.pdf')