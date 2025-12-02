
covariates="nCount_RNA,nFeature_RNA,percent.mt"
python3.10 run_FR_Perturb.py --input-h5ad CAF.h5ad \
	--perturbation-column-name Perturbation \
	--control-perturbation-name Control \
	--compute-pval \
	--out CAFs \
	--covariates nCount_RNA,nFeature_RNA,percent.mt \
	--large-dataset \
	--guide-pooled \
	--perturbation-delimiter :

