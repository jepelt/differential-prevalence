# Benchmarking DiPPER

This repository contains the code and datasets required to reproduce the
analyses and figures in the manuscript introducing DiPPER (preprint available
[here](https://arxiv.org/abs/2602.05938)).

**Note:** The development version of the **DiPPER** R package itself is hosted
in a separate repository
([https://github.com/jepelt/DiPPER](https://github.com/jepelt/DiPPER)),
which includes an installation guide and a usage example.

## Reproducing the Paper Analyses

All the code and datasets used for benchmarking in the paper are located in
the `R` folder.

### Data processing
Raw data was processed using the `data_tses_16s_030625.R` and
`data_tses_sg_030625.R` scripts. The raw data files for the 16S data can be
downloaded from [Zenodo](https://zenodo.org/records/1146764).

The processed data are stored in the files `data_tses_16s_030625.rds` and
`data_tses_sg_030625.rds` as lists of `TreeSummarizedExperiment` objects.
Further data processing was done using the `data_090625.R` script (which
created an additional file), and the `data_10...R` scripts were used to split
the large data objects into smaller files to enable parallel benchmarking on
the CSC cluster.

### Benchmarking DiPPER
The script `run_dipper_default_110526.R` was used to run the default version of
DiPPER on the benchmarking datasets. These computations were executed on the
CSC (IT Center for Science Ltd) high-performance computing cluster. The results
from running DiPPER are saved in the file `res_dipper_default_110526.rds`.

Versions of DiPPER with alternative hyperprior settings were run using the
`run_dipper_symmetric/s1/s2/s3_110526.R` scripts, and the corresponding
results are stored in the files `res_dipper_symmetric/s1/s2/s3_110526.rds`.

### Running the competing methods
The competing methods were executed locally on a standard computer using
their respective `run_[METHOD_NAME]...R` scripts. As an exception, LDM was run
on the CSC cluster using the script `run_ldm_120525_csc_updated.R`, and the
results are stored in `res_ldm_120526.rds`. Moreover, the "Bayesian" MaAsLin2
was run using the `run_blinrall_301225.r` script (which utilizes the Stan code
`stan_blinrall_301225.stan`).

### Summarizing results and creating figures
The benchmarking results and Figures 3–5 presented in the Results section of
the paper can be generated using the corresponding `fig_[FIGURE_NUMBER]...R`
scripts. The scripts to generate the figures in the Supplementary Material of
the manuscript are named `fig_s[FIGURE_NUMBER]...R`.