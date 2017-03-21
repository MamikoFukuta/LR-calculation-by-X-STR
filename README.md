LR-calculation-by-X-STR
====

This program can calculate likelihood ratio of kinship test by X-chromosomal short tandem repeats (X-STR) using R.  
 
## Description
A simple calculation program to estimate likelihood ratio (LR) of four types of kinship test with X-STR incorporating linkage, linkage disequilibrium, and mutation.  
Four relationship types are following;  

*    father-daughter (FD) vs. unrelated individuals (UR) without data of mother  
*    full-sister (FS) vs. UR without data of parents  
*    FS vs. maternal half-sisters (MHS) without data of parents  
*    paternal half-sisters (PHS) vs. UR  without data of parents  

Any set of X-STR can be used if you have allele frequency data.  
All program were written by R version 3.1.2.  


## Requirement
*    R  (The free statistical software. Download from http://www.R-project.org.)  
*    Text file data of allele frequency (and haplotype frequency) and recombination rate of X-STR. Example files are uploaded.  
*    Genotype data for kinship test (Also, you can make a sample data by "simulation_sample_generator.R")  

## Code
*    simulation_sample_generator  
*    LR_FDUR  
*    LR_FSUR  
*    LR_FSMHS  
*    LR_PHSUR  


