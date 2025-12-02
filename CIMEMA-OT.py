
import warnings

warnings.filterwarnings("ignore")

import matplotlib.pyplot as plt
import numpy as np
import pertpy as pt
import scanpy as sc
import os
import pandas as pd

plt.rcParams["figure.dpi"] = 250
plt.rcParams["font.size"] = 15
import anndata as ad

paths = "~/CINEMA_OT"
if not os.path.exists(paths):
    os.makedirs(paths)
    print(f"The path {paths} has been created")
os.chdir(paths)


adata1 = sc.read_h5ad('rdat_rds.h5ad')

subset_control = adata1[adata1.obs['Pertubation'].isin(['Control'])]

# %%
pertub=pd.read_csv('Perturbation_module.csv',header=0,index_col=0)


# For one perturbation

ubset_slc25a41= adata1[adata1.obs['Pertubation'].str.contains('SMPD3')]

perturb='SMPD3'

adatas = {}
adatas["SMPD3"]=ubset_slc25a41
adatas["control"]=subset_control
adatac = ad.concat(adatas, label="Label")

cot = pt.tl.Cinemaot()
sc.pp.pca(adatac)
de = cot.causaleffect(
    adatac,
    pert_key="Label",
    control="control",
    return_matching=True,
    thres=0.5,
    smoothness=1e-5,
    eps=1e-3,
    solver="Sinkhorn",
    preweight_label="cell_states",
)
sc.pp.neighbors(adatac, use_rep="cf")
sc.tl.umap(adatac, random_state=1)
sc.pl.umap(adatac, color=["Label", "cell_states"], wspace=0.5)
sc.pp.neighbors(de, use_rep="X_embedding")
sc.tl.umap(de)
sc.pl.umap(de, color=["Label","cell_states"],wspace=0.5)

pfig=cot.plot_vis_matching(
    adatac,
    de,
    pert_key="Label",
    control="control",
    de_label=None,
    source_label="cell_states",
    normalize="row",
    min_val=0.05,return_fig=True
)
fig = pfig.get_figure()
hdata = pfig.collections[0].get_array().data
n_rows = len(pfig.get_yticklabels())
n_cols = len(pfig.get_xticklabels())
matrix = np.array(hdata).reshape(n_rows, n_cols)
y_labels = [t.get_text() for t in pfig.get_yticklabels()]
x_labels = [t.get_text() for t in pfig.get_xticklabels()]
matchMat = pd.DataFrame(matrix, index=y_labels, columns=x_labels)
matchMat.to_csv(perturb+"_match_matrix.csv")

