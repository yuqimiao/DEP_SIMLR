setwd("/ifs/scratch/msph/biostat/sw2206/yuqi")
.libPaths("/ifs/scratch/msph/biostat/sw2206/yuqi/R_libs")

BiocManager::install("minfi")
BiocManager::install("sva")
BiocManager::install("wateRmelon")
BiocManager::install("impute")
library(tidyverse)
library(minfi)
library(sva)
library(readr)
library(dplyr)
library(wateRmelon)
