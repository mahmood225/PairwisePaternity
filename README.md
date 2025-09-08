# PairwisePaternity
This repository contains R code implementing the Pairwise algorithm for paternity analysis, developed as part of our paper “Pairwise paternity assignment with forward–backward simulations: Refining CERVUS using trio-based likelihood and locus-specific error rates” Traditional parentage methods, such as CERVUS, often rely on simplifying assumptions like homogeneous genetic structure and constant error rates across loci. These assumptions are frequently violated in real-world datasets, leading to reduced accuracy or false assignments. To address these limitations, the Pairwise algorithm introduces case-specific significance thresholds for each father–mother–offspring trio using forward and backward simulations, while also accounting for locus-specific genotyping error rates.

Our simulations demonstrated that this approach significantly reduces false-positive assignments by improving true-negative accuracy, particularly under challenging conditions such as incomplete parental sampling, missing maternal genotypes, or when candidate parents are closely related. While CERVUS generally achieves high true-positive rates when the true parent is present in the candidate pool, the Pairwise algorithm provides a more robust framework in ecological and conservation contexts where sampling may be incomplete or genetic relatedness is high. By offering trio-specific likelihood distributions, the method improves reliability and interpretability of paternity inference.

Overall, the Pairwise algorithm complements existing tools such as CERVUS and COLONY by providing a lightweight, simulation-based framework for case-by-case paternity resolution. It is especially suited for studies requiring flexibility, computational efficiency, and robust control of false assignments. The R scripts in this repository reproduce the analyses in the paper and provide an open framework for further refinement and adaptation to both microsatellite (STR) and SNP-based parentage datasets.

## Workflow

The workflow consists of 6 main functions that process genotype data and estimate parentage:

1. **Function 1 (Case_Assignment_Simulation.R):** Short explanation.  
2. **Function 2 (name here):** Short explanation.  
3. **Function 3 (name here):** Short explanation.  
4. **Function 4 (name here):** Short explanation.  
5. **Function 5 (name here):** Short explanation.  
6. **Function 6 (name here):** Short explanation. 


