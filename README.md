# TEAM - Testing on an Aggregation Tree Method

&nbsp;

# Overview 

<p align="justify"> 

`TEAM` is a multiple-testing method that embeds hypothesis testing on an aggregation tree to test hypotheses from fine- to coarse- resolution, while controlling the false discovery rate (FDR). Specifically, we developed `TEAM` as a method to identify where two probability density functions (pdfs) differ. First, `TEAM` partitions the multivariate sample space into bins with the finest resolution. It can accommodate different partitioning schemes. Second, `TEAM` embeds testing on an aggregation tree with user-specified number of layers. The first layer bins are defined by the initial partition, and in each bin, `TEAM` tests if the pdf of one sample is higher than that of the other (e.g. reference). On higher layers, TEAM will gradually aggregate the non-rejected bins and test if the aggregated bins harbor differential pdfs. This fine- to coarse-resolution testing structure not only boosts the testing power, but also pinpoints the regions with differential pdfs at the finest possible resolution. We apply `TEAM` to a flow cytometry study that aims to identify regions of T cell activation, based on multivariate protein marker expression.
  
 </p>

![TEAM](/uploads/7475dbf686561699ebfc809904dd2d0d/TEAM.png)

# Getting started

## Installation

If the package TEAM is not installed, install it via CRAN:
  
```{r,eval=FALSE}
install.packages("TEAM")
```
## Documentation

[Insert vignette information]

## Basic Usage
<p align="justify">
TEAM requires three inputs: a partition on the sample space (see below for details), the maximum number of layers, $L$, for the aggregation tree, and the desired FDR level, $\alpha$.

The procedure can be run as follows:
</p>
 
```{r, eval=FALSE}
TEAM(partition_info=partition2D,
     L=L,
     alpha=alpha)
```

An example 11-color flow cytometry dataset from the External Quality Assurance Program Oversight Laboratory (EQAPOL) program can be loaded using the following call:

```{r, eval=FALSE}
data(EQAPOL_Ex)
```

# Paper

[Link to arxiv paper]

# Authors

John Pura, Duke University 
Xuechan Li, Duke University
Cliburn Chan, Duke University
Jixuen Xie, Duke University
