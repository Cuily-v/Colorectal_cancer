1. pyscenic grn --num_workers 1 --output ./1.tsv --method grnboost2  ./1.loom hs_hgnc_tfs.txt
2. pyscenic ctx ./1.tsv hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather --annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname ./1.loom --mode "dask_multiprocessing" --output ./1.csv --num_workers 1 --mask_dropouts
3. pyscenic aucell ./1.loom ./1.csv --output  ./1.loom --num_workers 1
