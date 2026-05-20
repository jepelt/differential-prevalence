## DiPPER - Differential Prevalence via Probabilistic Estimation in R

DiPPER is a Bayesian hierarchical model designed for differential prevalence analysis in microbiome studies. Unlike standard frequentist approaches (based on e.g., the Wald test), which may fail or yield infinite estimates in boundary cases (i.e., when a taxon is completely absent in one group), DiPPER produces robust, finite estimates through Bayesian regularization. Furthermore, the model provides differential prevalence estimates and uncertainty intervals that are inherently adjusted for multiplicity.

Technically, DiPPER utilizes a common asymmetric Laplace prior (whose variance and skewness are determined by the data) for the differential prevalence parameters. This choice is  motivated by a natural assumption that for most taxa in most studies the true differential prevalence effects are likely close to zero, and the observation that typically, within a given microbiome study, most of the non-zero prevalence differences have the same direction.

This repository contains the code and datasets (located in the R folder) for reproducing the analyses and figures in the paper introducing DiPPER (see a pre-print [here](https://arxiv.org/abs/2602.05938)).


### Example analysis with DiPPER

Please see an example [here](https://github.com/jepelt/DiPPER).
