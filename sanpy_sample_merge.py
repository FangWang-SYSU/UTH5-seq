import os
import numpy as np
import scanpy as sc
import anndata as ad
from scipy import sparse

def read_10x_sample(path, sample_id=None):

    path = os.path.abspath(path)
    if sample_id is None:
        sample_id = os.path.basename(path.rstrip("/"))

    h5 = os.path.join(path, "filtered_feature_bc_matrix.h5")
    mtx_dir = os.path.join(path, "filtered_feature_bc_matrix")

    if os.path.exists(h5):
        adata = sc.read_10x_h5(h5)  # 更快、省内存
    elif os.path.isdir(mtx_dir):
        adata = sc.read_10x_mtx(mtx_dir, var_names="gene_symbols", make_unique=True)
    else:
        # 兼容直接传到 filtered_feature_bc_matrix 目录
        if os.path.basename(path) == "filtered_feature_bc_matrix":
            adata = sc.read_10x_mtx(path, var_names="gene_symbols", make_unique=True)
        else:
            raise FileNotFoundError(f"Not found: {h5} or {mtx_dir}")

    
    adata.obs_names = sample_id + "-" + adata.obs_names.astype(str)
    adata.obs["sample"] = sample_id

   
    if sparse.issparse(adata.X):
        C = adata.X.tocsr()
    else:
        C = sparse.csr_matrix(adata.X)
    if not np.issubdtype(C.dtype, np.integer):
        C.data = np.rint(C.data).astype(np.int32, copy=False)
    else:
        C.data = C.data.astype(np.int32, copy=False)
    adata.layers["counts"] = C

    return adata


def merge_10x_samples(sample_dirs, sample_ids=None):

    ads = []
    for i, p in enumerate(sample_dirs):
        sid = sample_ids[i] if sample_ids is not None else None
        ads.append(read_10x_sample(p, sid))

    A = ad.concat(
        ads,
        join="outer",
        label="sample",
        merge="first",
        fill_value=0,
        index_unique=None,
    )
    A.X = A.layers["counts"]  
    return A