setwd("~/MyProjects/GMHI_snv_chenchen/")
source("BetweenGroup.test.R")
#I/O---------------------------------------------------
species_relative_abundance_file <-"data/genome_abd_1711.txt"#"./mutant_genome_SNV.txt"# "./mutant_genome_SNV.txt" #
snv_file <- "data/all_snv_rate_1711.txt"
metadata_file <- "data/metadata.txt"
prefix="species"#"genomeSNV"
output_file = 'GMHI_output.csv'

#-------------------------------
# install and load necessary libraries for data analyses
#-------------------------------
p <- c("reshape2","ggplot2", "dplyr", "ROCR", "parallel", "tidyverse", "foreach","viridis", "doMC", "plyr", "rlang")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="https://cloud.r-project.org/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))

species_profile <- read.table(species_relative_abundance_file, sep = "\t",
                            header = TRUE, row.names = 1, check.names = F) 
species_profile <- species_profile[, order(colnames(species_profile))]

#Load metadata
metadata <- read.table(metadata_file, sep = "\t",  header = TRUE,row.names = 1, check.names = F) 
metadata <- metadata[order(rownames(metadata)), ]

library(tidyverse)
t_species_profile <- data.frame(t(species_profile), check.rows = F,check.names = F)

# identify differentially abundant species features
sp_data_list<-filter_samples_by_sample_ids_in_metadata(t_species_profile, metadata)
sp_out <- BetweenGroup.test(sp_data_list$data, factor(sp_data_list$metadata$group), positive_class = "Healthy", q_cutoff = 1)
sink(paste(prefix, "_wilcox_out.tsv", sep=""))
cat("\t"); write.table(sp_out, quote=FALSE, sep="\t", row.names = TRUE)
sink()
#sp_out_sig <- sp_out %>% filter(IfSig=="Sig")
sp_out_sig <- sp_out %>% filter(non.param.test_p<0.05)

# SNV profile
snv_profile <- read.table(snv_file, sep = "\t", header = TRUE, row.names = 1, check.names = F) 
snv_profile <- snv_profile[, order(colnames(snv_profile))]

# Check if sample ids are consistent between feature table and metadata
t_snv_profile<- t(snv_profile)
snv_data_list<-filter_samples_by_sample_ids_in_metadata(t_snv_profile, metadata)
# Identify differentially abundant SNV features
snv_out <- BetweenGroup.test(snv_data_list$data, factor(snv_data_list$metadata$group), positive_class = "Healthy", q_cutoff = 1)
sink(paste("SNV_genome_wilcox_out.tsv", sep=""))
cat("\t"); write.table(snv_out, quote=FALSE, sep="\t", row.names = TRUE)
sink()
#sp_out_sig <- sp_out %>% filter(IfSig=="Sig")
snv_out_sig <- snv_out %>% filter(non.param.test_p<0.05)

MH_sp_list <- rownames(sp_out_sig)[which(grepl("Healthy_enriched", sp_out_sig[, "IfSigEnr"]))] # Health-prevalent species
MN_sp_list <- rownames(sp_out_sig)[which(grepl("Healthy_depleted", sp_out_sig[, "IfSigEnr"]))] # Health-scarce species

MH_snv_list <- rownames(snv_out_sig)[which(grepl("Healthy_enriched", snv_out_sig[, "IfSigEnr"]))] # Health-prevalent species
MN_snv_list <- rownames(snv_out_sig)[which(grepl("Healthy_depleted", snv_out_sig[, "IfSigEnr"]))] # Health-scarce species

t_species_snv_profile <- t_species_profile * t_snv_profile
spcies_snv_data_list<-filter_samples_by_sample_ids_in_metadata(t_species_snv_profile, metadata)
sp_snv_out <- BetweenGroup.test(spcies_snv_data_list$data, factor(spcies_snv_data_list$metadata$group), 
                             positive_class = "Healthy", q_cutoff = 1)
sink(paste("SNV_weigthed_genomeAbd_wilcox_out.tsv", sep=""))
cat("\t"); write.table(sp_snv_out, quote=FALSE, sep="\t", row.names = TRUE)
sink()
sp_snv_out_sig <- sp_snv_out %>% filter(non.param.test_p<0.05)
MH_sp_snv_list <- rownames(sp_snv_out_sig)[which(grepl("Healthy_enriched", sp_snv_out_sig[, "IfSigEnr"]))] # Health-prevalent species
MN_sp_snv_list <- rownames(sp_snv_out_sig)[which(grepl("Healthy_depleted", sp_snv_out_sig[, "IfSigEnr"]))] # Health-scarce species

GMHI <- function(x, MH_list, MN_list){
  # Extracting Health-prevalent species present in metagenome
  # Extracting Health-scarce species present in metagenome
  MH_x <- x[row.names(x) %in% MH_list, ]
  MN_x <- x[row.names(x) %in% MN_list, ]
  # Diversity among Health-prevalent species
  # Diversity among Health-scarce species
  alpha <- function(x){sum((log(x[x>0]))*(x[x>0]))*(-1)}
  MH_shannon <- apply((MH_x), 2, alpha) 
  MN_shannon <- apply((MN_x), 2, alpha) 
  
  # Richness of Health-prevalent species
  # Richness of Health-scarce species
  R_MH <- apply(MH_x, 2, function(i) (sum(i > 0))) 
  R_MN <- apply(MN_x, 2, function(i) (sum(i > 0)))
  
  # Median RMH from 1% of the top-ranked samples (see Methods)
  Median_percentile <- function(m, percentile, top_rank=TRUE){
    n <- round(max(length(m)) * percentile, 0)
    if(top_rank){
      sub_m <- m[rank(-m, ties.method="random") < n]
    }else{
      sub_m <- m[rank(m, ties.method="random") < n]
    }
    ifelse(median(sub_m)==0, 1, median(sub_m))
  }
  MH_prime <- Median_percentile(R_MH, 0.01, top_rank = TRUE)
  MN_prime <- Median_percentile(R_MN, 0.01, top_rank = FALSE)
  # Median RMH from 1% of the top-ranked samples (see Methods)
  # Median RMN from 1% of the bottom-ranked samples (see Methods)
  #MH_prime <- 7
  #MN_prime <- 31

  # Collective abundance of Health-prevalent species
  # Collective abundance of Health-scarce species
  psi_MH <- ((R_MH/MH_prime)*MH_shannon) 
  psi_MN <- ((R_MN/MN_prime)*MN_shannon)
  GMHI <- data.frame(gmhi=log10((psi_MH+0.00001)/(psi_MN+0.00001))) # 0.00001 added to avoid having the denominator as 0
  GMHI
}

GMHI_sp_df<-GMHI(x=species_profile, MH_sp_list, MN_sp_list)
GMHI_snv_df<-GMHI(x=snv_profile, MH_snv_list, MN_snv_list) 
GMHI_sp_snv_df<-GMHI(x=species_snv_profile, MH_sp_snv_list, MN_sp_snv_list)

y = as.numeric(factor(metadata$group))
gmhi_sp <- GMHI_sp_df$gmhi
gmhi_snv <- GMHI_snv_df$gmhi
gmhi_sp_snv <- GMHI_sp_snv_df$gmhi
outdir <- "./"

prefix=""
library("pROC")
# genome_seq_abd
rocobj <- pROC::plot.roc(y, gmhi_sp, main="", percent=TRUE,ci=TRUE) # print the AUC (will contain the CI)
ciobj <- pROC::ci.se(rocobj, specificities=seq(0, 100, 5)) # over a select set of specificities
dat.ci <- data.frame(x = as.numeric(gsub("%", "", rownames(ciobj))),
                     lower = ciobj[, 1],
                     upper = ciobj[, 3])
p<-ggroc(rocobj) + theme_minimal() + 
  geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=0.7, color = "grey") + 
  coord_equal() + 
  geom_ribbon(data = dat.ci, aes(x = x, ymin = lower, ymax = upper), fill = "steelblue", alpha= 0.2) + 
  ggtitle(capture.output(rocobj$auc)) 
p
ggsave(filename=paste(outdir, prefix, "gmhi_genome_seq_abd.pROC.ci.pdf",sep=""), plot=p, width=4, height=4)

# snv_rate
rocobj <- pROC::plot.roc(y, gmhi_snv, main="", percent=TRUE,ci=TRUE) # print the AUC (will contain the CI)
ciobj <- pROC::ci.se(rocobj, specificities=seq(0, 100, 5)) # over a select set of specificities
dat.ci <- data.frame(x = as.numeric(gsub("%", "", rownames(ciobj))),
                     lower = ciobj[, 1],
                     upper = ciobj[, 3])
p<-ggroc(rocobj) + theme_minimal() + 
  geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=0.7, color = "grey") + 
  coord_equal() + 
  geom_ribbon(data = dat.ci, aes(x = x, ymin = lower, ymax = upper), fill = "steelblue", alpha= 0.2) + 
  ggtitle(capture.output(rocobj$auc)) 
p
ggsave(filename=paste(outdir, prefix, "gmhi_snv_rate.pROC.ci.pdf",sep=""), plot=p, width=4, height=4)

# snv_weighted_abd

rocobj <- pROC::plot.roc(y, gmhi_sp_snv, main="", percent=TRUE,ci=TRUE) # print the AUC (will contain the CI)
ciobj <- pROC::ci.se(rocobj, specificities=seq(0, 100, 5)) # over a select set of specificities
dat.ci <- data.frame(x = as.numeric(gsub("%", "", rownames(ciobj))),
                     lower = ciobj[, 1],
                     upper = ciobj[, 3])
p<-ggroc(rocobj) + theme_minimal() + 
  geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=0.7, color = "grey") + 
  coord_equal() + 
  geom_ribbon(data = dat.ci, aes(x = x, ymin = lower, ymax = upper), fill = "steelblue", alpha= 0.2) + 
  ggtitle(capture.output(rocobj$auc)) 
p
ggsave(filename=paste(outdir, prefix, "gmhi_snv_weighted_abd.pROC.ci.pdf",sep=""), plot=p, width=4, height=4)

#

if (file.exists(output_file)){
  file.remove(output_file)
}

out <- data.frame(gmhi_sp, gmhi_snv, gmhi_sp_snv, snv_data_list$metadata)

write.csv(out, file=output_file) # Saving GMHI results as 'GMHI_output.csv'. User should change the path to the appropriate directory.

# End