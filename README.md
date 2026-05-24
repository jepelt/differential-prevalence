# DiPPER - Differential Prevalence via Probabilistic Estimation in R

DiPPER is a Bayesian hierarchical approach designed for differential prevalence
analysis in microbiome studies. Unlike standard frequentist approaches (based
on e.g., the Wald test), which may fail or yield infinite estimates in boundary
cases (e.g., when a taxon is completely absent in one group), DiPPER produces
robust, finite estimates through hierarchical shrinkage. Furthermore, DiPPER
provides differential prevalence estimates and uncertainty intervals that are
inherently adjusted for multiplicity.

Technically, DiPPER utilizes a common asymmetric Laplace prior (whose variance
and skewness are determined by the data) for the differential prevalence
parameters. This choice is motivated by a natural assumption that for most taxa
in most studies the true differential prevalence effects are likely close to
zero, and the observation that typically, within a given microbiome study, most
of the non-zero prevalence differences have the same direction.

This repository contains the code and datasets (located in the `R` folder) for
reproducing the analyses and figures in the paper introducing DiPPER (see a
pre-print [here](https://arxiv.org/abs/2602.05938)).

## Reproducing the paper analyses

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

## Installation of DiPPER and example usage

You can install the development version of DiPPER from GitHub using:

```r
# install.packages("remotes")
remotes::install_github("jepelt/DiPPER")
```

DiPPER also requires the `cmdstanr` package.

```r
install.packages("cmdstanr",
                 repos = c("(https://stan-dev.r-universe.dev/)",
                           getOption("repos")))

# Set up the C++ toolchain (Windows users may need Rtools)
cmdstanr::check_cmdstan_toolchain(fix = TRUE)

# Install the Stan backend (only needs to be done once)
cmdstanr::install_cmdstan()
```

### Example usage

Below is a simple example workflow using the example data included in the
package.

```r
library(DiPPER)

# Load example data (TreeSummarizedExperiment object)
# This dataset compares (N = 20 + 20) rats on a High/Low fat diet.
data("tse_hintikka")

# Run DiPPER. 
# The first term in the formula (here: Fat) is automatically 
# used as the variable of interest. XOS (xylo-oligosaccharide supplementation)
# is included as a covariate to adjust for.
# Note: When using DiPPER for the first time, it may take around two minutes
# to run the function due to the compilation of the Stan model.

fit <- DiPPER(
  tse = tse_hintikka,
  formula = ~ Fat + XOS,
  tax_rank = "Genus",
  seed = 1
)

# Extract summarized results as a data.frame
res <- summary(fit)

# Create a forest plot (showing only 'significant' taxa)
plot(fit, show_taxa = "significant")
```