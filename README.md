# DAP-PIR: Deterministic Approximation of Posteriors with Pseudo Importance Resampling

## Overview

In genetic fine mapping, identifying causal variants remains a key challenge due to computational limitations. Here, we introduce a novel method, DAP-PIR, which builds upon the Sum of Single Effects (SuSiE) results by integrating pseudo-importance resampling (PIR) to efficiently explore models with high probability. Specifically, DAP-PIR leverages SuSiE's variational approximation to identify high-priority SNPs in independent effects and refines posterior probabilities through the deterministic approximation of posteriors (DAP) based on pseudo importance resampling (PIR) strategy.

## Installation

The `dap` package implements the algorithm of DAP-PIR. You can install the development version of `dap` from GitHub.

``` r
# install devtools
install.packages("devtools")

# install DAP-PIR package
devtools::install_github("wnbo9/DAP-PIR/dap")

# load the package
library(dap)
```

## Usage

Check the [vignettes](https://github.com/wnbo9/dap-pir/tree/main/dap/vignettes/).
