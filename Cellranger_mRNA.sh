# for each one mRNA transcriptome
cellranger count \
 --id=index1 \
 --transcriptome=refdata-gex-GRCh38-2020-A \
 --fastqs=input \
 --localcores=30
--sample=index1