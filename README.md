# PhenotypeGenotypeModeling

This is a C++ implementation of a mathematical model describing the relationship between genotype and phenotype, which evolve under natural selection by an environment.

- The model was presented by
    - Laarits, T., P. Bordalo, and B. Lemos. "Genes under weaker stabilizing selection increase network evolvability and rapid regulatory adaptation to an environmental shift." Journal of Evolutionary Biology 29.8 (2016): 1602-1616.
    Available at: https://doi.org/10.1111/jeb.12897

## A Brief Model Description

- In the model, genotypes are represented by matrices, whilst phenotypes are represented by associated vectors. That is, for a given population, each individual is represented by a matrix and its associated vector.

- The environment which exerts influence in the population is also represented as vector named target.

- Reproduction of the population over generations occurs by random selection of pairs, with the probability of being selected proportional to the fitness of the individual. 

- The fitness is quantified by comparing the phenotype with the target. The closer an individual is to the target, the more it is probable that it will be selected for reproduction.

- New matrices are constructed by picking rows from both parents. After new matrices are built, they may undergo mutations. 

- In the model, mutations are the substitution of one of the matrices' entries.

## Files in the repository

The bash script `compile_and_execute.sh` will compile the file `pgtype_experiment.cpp`using *g++*. The **Eigen** library is a requirement of this implementation.
The script will also run the compiled executable, which in the script's configuration will run 20 isolated experiments, with 8000 generations each.
Only the final mean values of phenotypes and fitness will be saved after each experiment.

This configuration can be changed, you can run a `pgtype.exe -help` command to learn which configurations you can change at execution without having to change the code itself and later recompiling it.

## References

- Laarits, T., P. Bordalo, and B. Lemos. "Genes under weaker stabilizing selection increase network evolvability and rapid regulatory adaptation to an environmental shift." Journal of Evolutionary Biology 29.8 (2016): 1602-1616.
    - Available at: https://doi.org/10.1111/jeb.12897

- Guennebaud, Gaël, and Jacob, Benoit and others. "Eigen v3" (2010)
    - Available at: http://eigen.tuxfamily.org 

