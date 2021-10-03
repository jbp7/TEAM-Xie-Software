# TEAM - Testing on an Aggregation Tree Method

&nbsp;

# Overview 

<p align="justify"> 

`TEAM` is a multiple-testing method that embeds hypothesis testing on an aggregation tree to test hypotheses from fine- to coarse- resolution, while controlling the false discovery rate (FDR). Specifically, we developed `TEAM` as a method to identify where two probability density functions (pdfs) differ. First, `TEAM` partitions the multivariate sample space into bins with the finest resolution. It can accommodate different partitioning schemes. Second, `TEAM` embeds testing on an aggregation tree with user-specified number of layers. The first layer bins are defined by the initial partition, and in each bin, `TEAM` tests if the pdf of one sample is higher than that of the other (e.g. reference). On higher layers, TEAM will gradually aggregate the non-rejected bins and test if the aggregated bins harbor differential pdfs. This fine- to coarse-resolution testing structure not only boosts the testing power, but also pinpoints the regions with differential pdfs at the finest possible resolution. We apply `TEAM` to a flow cytometry study that aims to identify regions of T cell activation, based on multivariate protein marker expression.
  
 </p>

![TEAM](/uploads/7475dbf686561699ebfc809904dd2d0d/TEAM.png)

# Getting started

## Installation

TEAM may be installed from GitLab via the `remotes` package:
  
```{r,eval=FALSE}
install.packages("remotes")
remotes::install_gitlab(
  repo = "jichunxie/Xie-lab-software_TEAM",
  subdir = "TEAM",
  host = "gitlab.oit.duke.edu",
  build_vignettes = TRUE,
  build_manual = TRUE)
```

## Documentation

After installation, a detailed vignette containing examples of how to use TEAM can be viewed by calling:
```{r,eval=FALSE}
browseVignettes("TEAM")
```

## Basic Usage
<p align="justify">
TEAM requires three inputs: a partition on the sample space (see below for details), the maximum number of layers, $L$, for the aggregation tree, and the desired FDR level, $\alpha$.

The procedure can be run as follows:
</p>
 
```{r, eval=FALSE}
TEAM(partition_info=partition,
     L=L,
     alpha=alpha)
```

Please see the vignette for details in constructing the partition and the output of the function.

An example flow cytometry dataset from the External Quality Assurance Program Oversight Laboratory (EQAPOL) program can be loaded from our GitLab repository by calling:

```{r, eval=FALSE}
URL <- "https://gitlab.oit.duke.edu/jichunxie/xie-lab-software_team/-/raw/main/TEAM/data/example_eqapol_data.rda"
download.file(URL,"example_eqapol_data.rda", method="curl")
load("example_eqapol_data.rda")
```

# Paper

`TEAM` accompanies the paper:
   
[Pura J., Li X., Chan C., Xie J. 2021. TEAM: A Multiple Testing Algorithm on the Aggregation Tree for Flow Cytometry Analysis](https://arxiv.org/abs/1906.07757)


# Authors

John Pura, Duke University 

Xuechan Li, Duke University

Cliburn Chan, Duke University

Jixuen Xie, Duke University
