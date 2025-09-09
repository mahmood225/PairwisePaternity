# PairwisePaternity
This repository contains R code implementing the Pairwise algorithm for paternity analysis, developed as part of our paper “Pairwise paternity assignment with forward–backward simulations: Refining CERVUS using trio-based likelihood and locus-specific error rates” Traditional parentage methods, such as CERVUS, often rely on simplifying assumptions like homogeneous genetic structure and constant error rates across loci. These assumptions are frequently violated in real-world datasets, leading to reduced accuracy or false assignments. To address these limitations, the Pairwise algorithm introduces case-specific significance thresholds for each father–mother–offspring trio using forward and backward simulations, while also accounting for locus-specific genotyping error rates.

Our simulations demonstrated that this approach significantly reduces false-positive assignments by improving true-negative accuracy, particularly under challenging conditions such as incomplete parental sampling, missing maternal genotypes, or when candidate parents are closely related. While CERVUS generally achieves high true-positive rates when the true parent is present in the candidate pool, the Pairwise algorithm provides a more robust framework in ecological and conservation contexts where sampling may be incomplete or genetic relatedness is high. By offering trio-specific likelihood distributions, the method improves reliability and interpretability of paternity inference.

Overall, the Pairwise algorithm complements existing tools such as CERVUS and COLONY by providing a lightweight, simulation-based framework for case-by-case paternity resolution. It is especially suited for studies requiring flexibility, computational efficiency, and robust control of false assignments. The R scripts in this repository reproduce the analyses in the paper and provide an open framework for further refinement and adaptation to both microsatellite (STR) and SNP-based parentage datasets.

## Workflow

The workflow consists of 6 main functions that process genotype data and estimate parentage:

### Simulating STR Population

The function `Historical_Simulation()` generates a multi-generational STR-based population, following the genomic data simulation design from our paper.
This function generates the historical population used as the foundation for paternity analysis simulations. It creates multiple generations of individuals with short tandem repeat (STR) genotypes under Mendelian inheritance. The simulation incorporates mutation processes, allele frequency stabilization, and decreasing population sizes across historical generations to establish realistic genetic diversity. The output provides a stable base population with allele frequency distributions that can be used for subsequent forward or backward simulations in the Pairwise algorithm.

```R
# Load libraries
install.packages(c("doParallel","foreach","data.table"))  # once
library(doParallel)
library(foreach)

# Source the function files you need
source("Historical_Simulation.R")

#Run the simulation
set.seed(123)  # reproducibility

simulated_STR <- Historical_Simulation(
# STR parameters
Num_SRT=15, # Number of simulated STR
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

# Add name for simulated STRs
colnames(STR)=paste0("STR",1:15)

# Calculating the allele freq
Allele_freq <- list()
for (i in which(!colnames(STR)%in%c("ID","Sex", "DamID", "SireID", "Date"))) {
  p=table(unlist(strsplit(STR[,i], ""), use.names=FALSE))/sum(table(unlist(strsplit(STR[,i], ""), use.names=FALSE)))
  Allele_freq[[length(Allele_freq)+1]] <- list(p)
}
names(Allele_freq)<-colnames(STR)[which(!colnames(STR)%in%c("ID","Sex", "DamID", "SireID", "Date"))]


Allele_freq[["STR1"]]

```

### Typing Error Calculation

The function `Typing_Error_Calculation()` estimates locus-specific typing errors after the simulated population has been generated.  
In real genetic studies, typing errors are common and their exact rates are usually unknown. Here, an error means that the true paternal or maternal allele at a locus is replaced by a random allele. These errors can appear in the genotypes of fathers, mothers, or offspring.  

To approximate realistic error rates, the function looks for mismatches between parents and offspring across the simulated population (except the last two generations, which are kept for validation). It then counts how often mismatches occur at each STR locus, giving you an estimated error rate per locus that can be used in later paternity tests.

```R
# Source the function files you need
source("./Typing_Error_Calculation.R")

# Combine pedigree and STR genotypes into one data frame
STR <- cbind(ped, STR)

# Estimate locus-specific typing error rates
Typing_Error <- Typing_Error_Calculation(
  STR = STR,
  Min_genotyped_STR = 7,              # Minimum number of STR loci that must be genotyped in an individual to include it in error estimation
  Correct_Genotype_Probability = 0.8, # Proportion of correctly assignable genotypes assumed to be true when calculating errors
  Minimum_Typing_Error = 0.001,       # Lower bound: if the estimated error is smaller, it is set to this minimum value
  Ensure_Correctness = FALSE,         # If TRUE, assumes assigned parents are always correct when calculating error rates
  Mother_Offspring_Calulation = TRUE, # Use mother–offspring comparisons to estimate errors
  Father_Offspring_Calulation = TRUE  # Use father–offspring comparisons to estimate errors
)

Typing_Error
```

### Forward Simulation – Paternity Testing

The function `Forward_Simulation()` is used to check whether an alleged father is the true biological father.  
It works by creating many simulated offspring genotypes through Mendelian inheritance (sampling alleles from the parents).  

- If only the father’s genotype is known, the missing maternal alleles are generated based on population allele frequencies.  
- Each simulated offspring is then adjusted to include realistic features like typing errors or missing loci.  
- For every simulated offspring, the function calculates a LOD score for the alleged father.  
- Repeating this thousands of times builds a distribution of LOD scores that represent the **true father case**.  
- Finally, this distribution is used to test real offspring against a chosen significance level (P-value).  

This makes the paternity test more reliable by providing case-specific thresholds instead of relying on population-wide assumptions.

```r
# Example: Running Forward Simulation for Paternity Test

# Load function
source("Forward_Simulation.R")

# Run forward simulation
forward_results <- Forward_Simulation(
  STR = STR,
  ped = ped,
  Allele_freq = Allele_freq,
  Typing_Error = Typing_Error,
  Num_Test = 5000,                   # Number of simulated offspring
  Min_genotyped_STR = 5,             # Minimum STR loci required per sample
  Missing_parent = TRUE,             # Allow cases where mother is unknown
  n.cores = parallel::detectCores()-1
)

# View first few results
head(forward_results)
```

### Backward Simulation – Paternity Assignment

The function `Case_Assignment_Simulation()` is used when the goal is to **identify the true father among multiple candidate males**.  
Instead of generating offspring from known parents, this approach works in reverse:

- It creates simulated genotypes for possible fathers by sampling alleles from the offspring (and excluding maternal alleles if the mother is known).  
- Missing alleles are filled in using population allele frequencies.  
- Each simulated father is then modified to include realistic features like typing errors or missing loci.  
- For every simulated candidate, the function calculates a LOD score.  
- Repeating this many times builds a distribution of LOD scores for **true fathers**, which is then compared against the observed candidates.  

This simulation helps determine which candidate father is most likely to be the biological parent, while accounting for errors and incomplete information.

```r
# Example: Running Backward Simulation for Paternity Assignment

# Load function
source("./Case_Assignment_Simulation.R")

# Run backward simulation (assignment)
result_Case_sim_assignment <- Case_Assignment_Simulation(
  STR = STR,
  Cases = selected_ID,            # Offspring to test for assignment
  Typing_Error = Typing_Error,    # Locus-specific error rates
  Allele_freq = Allele_freq,      # Allele frequencies
  Num_Test = 5000,                # Number of simulations
  Min_genotyped_STR = 1,          # Minimum number of STR loci per sample

  # Scenarios
  One_Parent_Unknown_Another = TRUE,  # Allow assignment when one parent is missing
  Father_Known_Mother = TRUE,         # Test assignment when mother is known
  Mother_Known_Father = TRUE,         # Test assignment when father is known

  # Parallel computing
  n.cores = parallel::detectCores() - 1
)

# View first few results
head(result_Case_sim_assignment)
```

### Case Study – Paternity Test When Mother Is Unknown

In this example, we test whether the alleged father of **offspring ID 2188** is the true biological father when the mother’s genotype is missing.  
We use the function `LOD_Both_Unknown()` to calculate the likelihood ratio (LOD) score, and compare it to the simulated Δ (delta) distributions from both the **paternity test** and the **assignment test**.

```r
# Load the function
source("./LOD_Both_Unknown.R")

# Extract genotypes for offspring and parents
Offs_STR   <- STR[STR$ID == 2188, colnames(STR) %in% names(Typing_Error)]
FaID       <- ped$SireID[ped$ID == 2188]
Father_STR <- STR[STR$ID == FaID, colnames(STR) %in% names(Typing_Error)]
MoID       <- ped$DamID[ped$ID == 2188]
Mother_STR <- STR[STR$ID == MoID, colnames(STR) %in% names(Typing_Error)]

# Run paternity test (mother unknown)
LOD <- LOD_Both_Unknown(
  Offs_STR,
  Father_STR,
  Allele_freq,
  Typing_Error
)

# Extract simulated Δ distributions
Delta_LOD_paternity  <- result_Case_sim_paternity$Cases_Fa_Un_Mo[,"2188"]
Delta_LOD_assignment <- result_Case_sim_assignment$One_Parent_Unknown_Another[,"2188"]

# Calculate p-values
pval_paternity  <- sum(Delta_LOD_paternity  <= LOD) / length(Delta_LOD_paternity)
pval_assignment <- sum(Delta_LOD_assignment <= LOD) / length(Delta_LOD_assignment)

pval_paternity
pval_assignment

```

**Results**

**Paternity test p-value:** 0.9878

**Assignment test p-value:** 0.477

At a significance threshold of **α = 0.01**, the paternity test clearly supports that the alleged father of individual 2188 is the **true father**.
The assignment test gives a lower value, but since the direct paternity test is highly significant, the conclusion remains that the tested male is indeed the biological father.


### Case Study – Paternity Test When Mother Is Known

In this example, we test whether the alleged father of **offspring ID 2188** is the true biological father **given that the mother’s genotype is known**.  
We use `LOD_One_Parent_Known()` to compute the LOD score for the alleged father while conditioning on the known mother.  
We then compare that LOD to the **case-specific** Δ (delta) distributions from both the **paternity test** and the **assignment test** to obtain p-values.

```r
# Load the function
source("./LOD_One_Parent_Known.R")

# Compute LOD for father given known mother
LOD <- LOD_One_Parent_Known(
  Offs_STR            = Offs_STR,     # offspring STR vector (e.g., ID 2188)
  Alleged_Parent_STR  = Father_STR,   # alleged parent STR vector (father)
  Known_Parent_STR    = Mother_STR,   # known parent STR vector (mother)
  Allele_freq         = Allele_freq,
  Typing_Error        = Typing_Error
)

# Extract case-specific Δ distributions (ID "2188")
Delta_LOD_paternity  <- result_Case_sim_paternity$Cases_Fa_Kn_Mo[,"2188"]
Delta_LOD_assignment <- result_Case_sim_assignment$Cases_Fa_Kn_Mo[,"2188"]

# Compute p-values
pval_paternity  <- sum(Delta_LOD_paternity  <= LOD) / length(Delta_LOD_paternity)
pval_assignment <- sum(Delta_LOD_assignment <= LOD) / length(Delta_LOD_assignment)

pval_paternity
pval_assignment
```
**Results**

**Paternity test p-value:** 0.8436

**Assignment test p-value:** 0.7118

**Interpretation:**
Both p-values are well above typical alpha thresholds for rejecting paternity (e.g., α = 0.01), which means the alleged father of individual 2188 is very likely to be the true father.
However, compared to the previous test case (when the mother was unknown, p-value ≈ 0.9878), these values are lower. This suggests that when the mother is included, the statistical support is still strong but not as extreme.

Overall, since 0.8436 > 0.01, we fail to reject the hypothesis that the alleged father is the true father, and the result supports biological paternity.
The assignment test p-value (0.7118) also indicates consistency but provides slightly weaker evidence compared to the direct paternity test.

### Case Study – Paternity Test for Both Parents Jointly

This test evaluates whether **both alleged parents (mother and father)** are the true biological parents **simultaneously**.  
Using `LOD_Parents_Jointly()`, we compute a joint LOD score that conditions on the genotypes of **both** parents and the offspring.  
We then compare this LOD to the **case-specific** Δ (delta) distribution from the joint paternity simulation.  
> Note: This mode is **paternity-only** (no assignment test), since it evaluates a specific parental pair rather than searching among multiple candidates.

```r
# Load the function
source("./LOD_Parents_Jointly.R")

# Compute joint LOD for mother+father together
LOD <- LOD_Parents_Jointly(
  Offs_STR     = Offs_STR,     # offspring STR vector (e.g., ID 2188)
  Mother_STR   = Mother_STR,   # alleged mother STR vector
  Father_STR   = Father_STR,   # alleged father STR vector
  Allele_freq  = Allele_freq,
  Typing_Error = Typing_Error
)

# Extract the case-specific Δ distribution (ID "2188") for the joint test
Delta_LOD_paternity <- result_Case_sim_paternity$Cases_Both_Jointly[,"2188"]

# Compute the p-value for the joint paternity test
pval_joint <- sum(Delta_LOD_paternity <= LOD) / length(Delta_LOD_paternity)
pval_joint
```

**Results**

**Joint paternity test p-value**: 1.0000

**Interpretation:**
With p = 1.0, the observed joint LOD is at (or beyond) the extreme upper tail of the simulated true-parent distribution.
At a stringent threshold such as α = 0.01, this provides overwhelming support that the tested pair (the alleged mother and father of ID 2188) are indeed the **true biological parents**.


The **same procedure** applies when testing whether the **alleged mother** is the true biological mother, both when the **father is unknown** and when the **father is known**.
