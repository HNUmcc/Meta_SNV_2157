This repository provided all source data and codes for generation of all results in the Meta_SNV_2157 manuscript.

Single nucleotide variation of gut microbiota suggest health assessment approach. Correspondence: Jiachao Zhang (E-mail: zhjch321123@163.com) and Shi Huang (E-mail: shihuang047@gmail.com)

The pipeline of GMHI on microbiome analyses:
1. Identification of strains with different abundance, SNV rate and SNV pattern of SCFAs gene in different host health status. Healthy-prevalent markers (MH) and healthy-scarce markers (MN) were obtained.
2. Calculating Shannon diversity (Hs and Ns) and richness (Hr and Nr) based on markers for each sample.
3. Identifying median Hr (Nr) from 1% of the top (bottom)-ranked samples (Hp and Np).
4. Collective abundance of markers, and then we can calculate:
psi_H = ((Hr/Hp)*Hs)
psi_N = ((Nr/Np)*Ns)
5. Calculating GMHI
GMHI= log10( psi_H/psi_N)
