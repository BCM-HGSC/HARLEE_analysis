# HARLEE_discovery

This repository contains the scripts and intermediate data files utilized in the manuscript "A Genocentric Approach to Discovery of Mendelian Disease Genes".

The folder `2018` and its subdirectories contain several scripts used to automate parameter-sweep querying of HARLEE (Hadoop ARchitecture LakE of Exomes). To protect the clinical samples contained within the HARLEE databases utilized for this project, [Impala](https://impala.apache.org) connection URLs and authorization credentials have been removed. Query automation is implemented in R, with Impala connections established through an R wrapper for JDBC – [RJDBC](https://CRAN.R-project.org/package=RJDBC). Query outputs – aggregate gene-level counts of samples with variants matching the specificied filtering parameters – are stored alongside query scripts.

The folder `data_lake_viz` contains downstream analysis and visualization scripts. These are implemented in the R language in [Jupyter](https://jupyter.org) notebooks. The file `data_lake_viz_functions.R` is a small R library of helper functions utilized across the notebooks. The folders `omim`, `uniprot`, and `possible_genes` contain annotation files and lists used for analysis.

