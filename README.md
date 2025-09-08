# PairwisePaternity
This repository contains R code implementing the Pairwise algorithm for paternity analysis, developed as part of our paper “Pairwise paternity assignment with forward–backward simulations: Refining CERVUS using trio-based likelihood and locus-specific error rates” Traditional parentage methods, such as CERVUS, often rely on simplifying assumptions like homogeneous genetic structure and constant error rates across loci. These assumptions are frequently violated in real-world datasets, leading to reduced accuracy or false assignments. To address these limitations, the Pairwise algorithm introduces case-specific significance thresholds for each father–mother–offspring trio using forward and backward simulations, while also accounting for locus-specific genotyping error rates.

Our simulations demonstrated that this approach significantly reduces false-positive assignments by improving true-negative accuracy, particularly under challenging conditions such as incomplete parental sampling, missing maternal genotypes, or when candidate parents are closely related. While CERVUS generally achieves high true-positive rates when the true parent is present in the candidate pool, the Pairwise algorithm provides a more robust framework in ecological and conservation contexts where sampling may be incomplete or genetic relatedness is high. By offering trio-specific likelihood distributions, the method improves reliability and interpretability of paternity inference.

Overall, the Pairwise algorithm complements existing tools such as CERVUS and COLONY by providing a lightweight, simulation-based framework for case-by-case paternity resolution. It is especially suited for studies requiring flexibility, computational efficiency, and robust control of false assignments. The R scripts in this repository reproduce the analyses in the paper and provide an open framework for further refinement and adaptation to both microsatellite (STR) and SNP-based parentage datasets.

## Workflow

The workflow consists of 6 main functions that process genotype data and estimate parentage:

1. **Historical simulation:**
This function generates the historical population used as the foundation for paternity analysis simulations. It creates multiple generations of individuals with short tandem repeat (STR) genotypes under Mendelian inheritance. The simulation incorporates mutation processes, allele frequency stabilization, and decreasing population sizes across historical generations to establish realistic genetic diversity. The output provides a stable base population with allele frequency distributions that can be used for subsequent forward or backward simulations in the Pairwise algorithm.

## Example: Simulating STR Population

The function `Historical_Simulation()` generates a multi-generational STR-based population, following the genomic data simulation design from our paper.

```R
# Load libraries
install.packages(c("doParallel","foreach","data.table"))  # once
library(doParallel)
library(foreach)
library(data.table)
# Source the function files you need
source("Historical_Simulation.R")

#Run the simulation
set.seed(123)  # reproducibility

Num_SRT <- 15

simulated_STR <- Historical_Simulation(
# STR parameters
Num_SRT=15,
Allele_num=rep(10,Num_SRT), # Number of the allele for each simulated STR
Allele_num_random=F, # If True, randome number of allele will be produced based on Min_allele_num and Max_allele_num
Typing_Error=c(rep(0.01,Num_SRT)), # typing error for each STR
Mutation_Rate=0.00001, # Mutation rate

# general parameters
Overlap_Gen_Num=2, # Maximum of overlapping generations for both parents
Overlap_Gen_Prop=c(0.6,0.3,0.1), # Overlapping generations, with a maximum of Overlap_Gen_Num +1 generations for both parents
Sire_Prop=0.5, #The proportion of males contributing to the next generation
Dam_Prop=0.8, # The proportion of females contributing to the next generation
Num_Offs=c(1,2),# vector to show possible number of offspring that can be seen in simulation
Prop_Num_Offs=c(0.95,0.05), # proportion for each possible number of offspring in Num_Offs
Sex_Prop=c(0.5,0.5),# first male and second female # historical parameters
Num_Hist_Start=2000, # Number of the initial generation (g = 0) of the historical population
Hist_Gen=50, # Number of historical generations
Num_Hist_End=200,# Number of the last generation (g = Hist_Gen) of the historical population

# Population Parameters
Num_Pop_Start=200, # Number of the initial generation of the studied population
Pop_Gen=10, # Number of generations for studied population
Num_Pop_End=200, # Number of the last generation of the studied population
n.cores = parallel::detectCores() - 2
)

# Extract outputs
ped <- simulated_STR$Simulated_Ped
STR <- simulated_STR$Simulated_STR
```


2) Register parallel workers

Historical_Simulation() takes n.cores, but registering the backend is still good practice (harmless if the function spins its own cluster).
3. **Function 2 (name here):** Short explanation.  
4. **Function 3 (name here):** Short explanation.  
5. **Function 4 (name here):** Short explanation.  
6. **Function 5 (name here):** Short explanation.  
7. **Function 6 (name here):** Short explanation. 


