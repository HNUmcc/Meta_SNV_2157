This repository provided all source data and codes for generation of all results in the Meta_SNV_2157 manuscript.

Go beyond “abundance”: cross-cohort single-nucleotide-variant profiling of gut microbiota suggests a novel gut-health assessment approach. Correspondence: Jiachao Zhang (E-mail: zhjch321123@163.com) and Shi Huang (E-mail: shihuang@hku.hk)

we generate three feature tables for building up a predictive index for health status (i.e., GMHI): (1) the sequence abundance profiles for all mutated genomes, (2) SNV rate for all mutated genomes (Eq 1); (3) SNV count for all mutant genes producing SCFAs for each sample.
GMHI based on SNV rate of mutated genomes. 
1. The health- and nonhealthy-enriched markers were identified for each feature table using multiple Wilcoxon rank-sum tests. The compositional data was central-log-ratio transformed prior to statistical tests. Each feature table was further filtered by healthy-enriched markers (MH) and healthy-depleted markers (MN) respectively.
2. In either MH or MN sub-table, we calculate Shannon diversity (Hs or Ns) and richness (Hr or Nr) based on markers for each sample.
3. We further calculated the median Hr (or Nr) from 1% of the top (bottom)-ranked samples from (Hp or Np).
4. We next calculated the “collective” sequence abundance, SNV rate of mutated genomes or SNV frequency of genes producing SCFAs for health-enriched (i.e., psi_H) or nonhealthy-enriched markers (i.e., psi_N) for each metagenomic sample.
psi_H=Hr/Hp×Hs (2)
psi_N=Nr/Np×Ns (3)
5. Calculating a GMHI for each sample.
GMHI=〖log〗_10 ((psi_H)/(psi_N)) (4)

