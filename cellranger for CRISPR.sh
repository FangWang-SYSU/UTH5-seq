# cellranger for CRISPR data
cellranger count \
 --id=index1 \
 --libraries=libraryFile.csv 
 --transcriptome=refdata-gex-GRCh38-2020-A \
 --feature-ref=GeCKO_v2_library.csv \
 --localcores=12

