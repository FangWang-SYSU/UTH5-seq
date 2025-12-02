import magic
import scprep
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import os

## input parameter 
inputmat='data_count_mat.csv'
outdir='~outdir'
prefix='test'
os.makedirs(outdir, mode=0o777, exist_ok=True)
os.chdir(outdir)

## counts matrix was as  input ,row was cell ,column was genes
print("---Data loading-----")
chunk = pd.read_csv(inputmat, chunksize=10000,index_col=0)
emt_data = pd.concat(chunk)
print('---normalization-----')
emt_data = scprep.normalize.library_size_normalize(emt_data)
emt_data = scprep.transform.sqrt(emt_data)
print('---imputation-----')
emt_magic = magic.MAGIC(solver='exact',n_jobs=40,random_state=1234)
emt_datapro = emt_magic.fit_transform(emt_data, genes='all_genes')
print('---writing to file-----')
emt_datapro.to_csv(prefix+"_expression_imputation_with_magic.csv")






