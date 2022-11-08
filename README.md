The application of catch-effort models to estimate the efficacy of
aerial shooting operations on sambar deer (*Cervus unicolor*)
================

## Overview

This repository contains data and code from:

Ramsey, D.S.L., McMaster, D., and Thomas, E. (2022). The application of
catch-effort models to estimate the efficacy of aerial shooting
operations on sambar deer (*Cervus unicolor*) *Wildlife Research*

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7302402.svg)](https://doi.org/10.5281/zenodo.7302402)

### File descriptions:

-   `r/sambar_removal.r` reads data and fits dynamic catch-effort
    N-mixture models to helicopter aerial shooting data from 10 sites
    using `nimble`.
-   `r/sambar_removal_JAGS.r` Same as above but fits the model using
    JAGS
-   `r/misc_functions.r` contains various functions require by the main
    script.

## Prerequisites

The script require packages `tidyverse`, `nimble`, `MCMCviz`, `abind`,
`jagsUI`.
