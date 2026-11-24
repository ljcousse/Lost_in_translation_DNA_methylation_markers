###############################################
###############################################
##                                           ##
##  Targeted Bisulphite sequencing analysis  ##
##  Sensitivity analyses                     ##
##                                           ##
###############################################
###############################################

# Sensitivity analyses and comparisons of the revision (Supplementary Figure 4, rebuttal).
# Runs on the objects saved at the end of TargetedBisulphite_analysis.R (output_for_sensitivity_analysis/): run that
# script first.

######################
# Parameter settings #
######################

options(warn=1L)

config_file <- Sys.getenv("PROJECT_CONFIG", "config.env")
if (!file.exists(config_file)) {
  stop("Config file not found: ", config_file,
       ". Copy config.example.env to config.env and edit the paths.")
}
readRenviron(config_file)

project_dir <- Sys.getenv("PROJECT_DIR")
if (!dir.exists(project_dir)) stop("PROJECT_DIR does not exist: ", project_dir)

# Objects of the main analysis
infile_sensitivity <- file.path(project_dir, "output_for_sensitivity_analysis", "objects_for_sensitivity_analysis.Rda")
if (!file.exists(infile_sensitivity)){
  stop("Objects of the main analysis not found: ", infile_sensitivity, ". Run TargetedBisulphite_analysis.R first.")
}

# Log file: console output (also shown on screen) + messages/warnings/errors (log only)
log_file <- file.path(project_dir, paste0("TargetedBisulphite_sensitivity_", format(Sys.time(), "%Y%m%d_%H%M"), ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
print(sprintf("Start: %s", format(Sys.time())))
print(sprintf("Log file: %s", log_file))
# Exactly one trailing "/" (paths below are built with paste(project_dir, ..., sep=""))
project_dir <- paste0(sub("/+$", "", project_dir), "/")
setwd(project_dir)

objects_loaded <- load(infile_sensitivity)
print(sprintf("Loaded %d objects of the main analysis: %s", length(objects_loaded), infile_sensitivity))

# Primer CpGs: excluded from the per-amplicon averages (as avAIC_table / Figure 4 of the main analysis)
CpGs_in_primer <- dat_correl$ID[(dat_correl$dist_CpG==0) & !is.na(dat_correl$dist_CpG)]
print(length(CpGs_in_primer))
# expect: 46

#############################
# Loading required packages #
#############################

library(BiSeq)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(limma)
library(biomaRt)
library(openxlsx)
library(stringr)
library(precrec)
library(scales)
library(parallel)
library(survival)
library(biomaRt)
library(survminer)
library(gridExtra)
library(ggplot2)
library(data.table)

##########################################
# Sensitivity NDRG4: without sample 178N #
##########################################

## Sensitivity: NDRG4 without 178N (unmatched; patient fixed effect -> expected no influence; R5, ToDo A4/E5)
#####################

sample_swap <- "178N"
mask_keep <- CRC_annot_sort_BS$Sample_ID != sample_swap
print(sprintf("Removing %s (%d/%d samples kept)", sample_swap, sum(mask_keep), length(mask_keep)))
CRC_annot_no178N <- CRC_annot_sort_BS[mask_keep, ]

formula_CF <- formula(paste("Tissue ~ Patient_ID",sep="")) # NDRG4 base model (formula_CF is redefined later in the main script)
model_CF_no178N <- glm(formula_CF, data = CRC_annot_no178N, family = binomial(link="logit"))
delta_AIC_no178N <- rep(NA, nrow(Mval_NDRG4_sort))
for (i in 1:nrow(Mval_NDRG4_sort)){
  formula_i <- formula(paste("Tissue~Meth_CpG_",i,"+Patient_ID",sep=""))
  model_i <- glm(formula_i, data = CRC_annot_no178N, family = binomial(link="logit"))
  delta_AIC_no178N[i] <- AIC(model_i) - AIC(model_CF_no178N)
}
print(sprintf("Max |Delta AIC difference| with vs without %s: %.4f", sample_swap,
              max(abs(delta_AIC_no178N - dat_NDRG4_ROC_logit_perCpG$delta_AIC), na.rm = T)))
# expect: ~0 (only numerical differences at the floor)


############################################################
# Sensitivity LY75: sequencing run as confounder (R3.M6a) #
############################################################
# Run 239 = additional melanoma samples; all CRC/ccRCC samples are in run 231, so only LY75 is affected.
# Non-nested as the main analysis (ToDo RV7, 2026-11-16): Meth + Run vs Gender + Age + Run, i.e. methylation replaces
# age + sex and Run is in both models; both on the same samples. Main results above unchanged.
levels_run <- levels(melanoma_annot_sort_BS$Run)
print(sprintf("Run levels: %s (reference: %s)", paste(levels_run, collapse = ", "), levels_run[1]))
print(table(melanoma_annot_sort_BS$Run, melanoma_annot_sort_BS$MK_death, dnn = c("Run", "Death melanoma")))
# expect: 2 levels; deaths 5/16 (run 231) vs 6/10 (run 239)

formula_CF_run <- formula("Surv(Followup_time, MK_death) ~ Gender+Age+Run")
mask_run_ref <- melanoma_annot_sort_BS$Run == levels_run[1]

dat_LY75_run <- dat_LY75_ROC_logit_perCpG[,c("ID","Amplicon","delta_AIC")]
dat_LY75_run$diff_Mval_run <- NA
dat_LY75_run$pval_run <- NA
dat_LY75_run$n_samples <- NA
dat_LY75_run$AIC_run <- NA
dat_LY75_run$AIC_CF_run <- NA
dat_LY75_run$delta_AIC_run <- NA
for (i in 1:nrow(Mval_LY75_sort)){
  mval_i <- melanoma_annot_sort_BS[,paste("Meth_CpG_",i,sep="")]

  # Batch check: M-value difference between runs (run 239 - run 231)
  n_ref_i <- sum(!is.na(mval_i[mask_run_ref]))
  n_other_i <- sum(!is.na(mval_i[!mask_run_ref]))
  if (n_ref_i > 1 & n_other_i > 1){
    dat_LY75_run$diff_Mval_run[i] <- mean(mval_i[!mask_run_ref], na.rm = T) - mean(mval_i[mask_run_ref], na.rm = T)
    dat_LY75_run$pval_run[i] <- t.test(mval_i[!mask_run_ref], mval_i[mask_run_ref])$p.value
  }

  # Cox models with Run, on the same samples (non-nested: methylation + Run vs Gender + Age + Run)
  formula_i_run <- formula(paste("Surv(Followup_time, MK_death) ~ Meth_CpG_",i,"+Run",sep=""))
  mask_complete_i <- complete.cases(melanoma_annot_sort_BS[,c(paste("Meth_CpG_",i,sep=""),"Followup_time","MK_death","Gender","Age","Run")])
  cox_run_i <- coxph(formula_i_run, data = melanoma_annot_sort_BS[mask_complete_i,])
  cox_CF_run_i <- coxph(formula_CF_run, data = melanoma_annot_sort_BS[mask_complete_i,])
  dat_LY75_run$n_samples[i] <- sum(mask_complete_i)
  dat_LY75_run$AIC_run[i] <- AIC(cox_run_i)
  dat_LY75_run$AIC_CF_run[i] <- AIC(cox_CF_run_i)
  dat_LY75_run$delta_AIC_run[i] <- dat_LY75_run$AIC_run[i] - dat_LY75_run$AIC_CF_run[i]
}
dat_LY75_run$fdr_run <- p.adjust(dat_LY75_run$pval_run, method = "BH")
write.table(dat_LY75_run,"LY75_AIC_table_run.txt",col.names = T,row.names = F,quote = F, sep="\t")

print(sprintf("LY75 CpGs with M-value difference between runs (FDR < 0.05): %d/%d",
              sum(dat_LY75_run$fdr_run < 0.05, na.rm = T), sum(!is.na(dat_LY75_run$fdr_run))))
print(sprintf("LY75 CpGs Delta AIC < -6: %d without Run, %d with Run (non-nested; of %d)",
              sum(dat_LY75_run$delta_AIC < -6, na.rm = T), sum(dat_LY75_run$delta_AIC_run < -6, na.rm = T), nrow(dat_LY75_run)))
print(sprintf("LY75 min Delta AIC: %.2f without Run, %.2f with Run",
              min(dat_LY75_run$delta_AIC, na.rm = T), min(dat_LY75_run$delta_AIC_run, na.rm = T)))
print(cor.test(dat_LY75_run$delta_AIC, dat_LY75_run$delta_AIC_run, method = "spearman", exact = F))

# Per amplicon (as Figure 4: mean and minimum over the target-sequence CpGs of the amplicon)
mask_target_run <- !(dat_LY75_run$ID %in% CpGs_in_primer)
avAIC_LY75_run <- data.frame(Amplicon = unique(dat_LY75_run$Amplicon), stringsAsFactors = F)
avAIC_LY75_run$avAIC <- sapply(avAIC_LY75_run$Amplicon, function(x){mean(dat_LY75_run$delta_AIC[dat_LY75_run$Amplicon==x & mask_target_run])})
avAIC_LY75_run$avAIC_run <- sapply(avAIC_LY75_run$Amplicon, function(x){mean(dat_LY75_run$delta_AIC_run[dat_LY75_run$Amplicon==x & mask_target_run])})
avAIC_LY75_run$minAIC <- sapply(avAIC_LY75_run$Amplicon, function(x){min(dat_LY75_run$delta_AIC[dat_LY75_run$Amplicon==x & mask_target_run])})
avAIC_LY75_run$minAIC_run <- sapply(avAIC_LY75_run$Amplicon, function(x){min(dat_LY75_run$delta_AIC_run[dat_LY75_run$Amplicon==x & mask_target_run])})
print(avAIC_LY75_run)
write.table(avAIC_LY75_run,"LY75_avAIC_run.txt",col.names = T,row.names = F,quote = F, sep="\t")

# MSP with Run (separation, ToDo A18: Delta AIC stays a bound), non-nested as above, on the same samples
mask_complete_MSP_run <- complete.cases(melanoma_annot_sort_BS[,c("Meth","Followup_time","MK_death","Gender","Age","Run")])
cox_MSP_run <- coxph(formula("Surv(Followup_time, MK_death) ~ Meth+Run"), data = melanoma_annot_sort_BS[mask_complete_MSP_run,])
cox_CF_run <- coxph(formula_CF_run, data = melanoma_annot_sort_BS[mask_complete_MSP_run,])
print(summary(cox_MSP_run))
print(sprintf("Delta AIC for MSP amplicon with Run (non-nested): %s (n = %d)", AIC(cox_MSP_run) - AIC(cox_CF_run), sum(mask_complete_MSP_run)))


###########################################################################
# Sensitivity LY75 + GREM1: unadjusted and clinical covariates (R3.M2)   #
###########################################################################
# Main models (TargetedBisulphite_analysis.R): Meth vs Gender + Age (non-nested). Here, on the same samples per comparison:
#   unadjusted: Meth vs empty model; clinical: Meth + Gender + Age + clinical covariates vs Gender + Age + clinical covariates.
# With 11 (melanoma) / 15 (ccRCC) events these models are small-sample sensitivity analyses only.

sensitivity_cox <- function(dat, n_cpg, dat_perCpG, event, covariates, gene){
  cn_time_event <- c("Followup_time", event)
  formula_null <- formula(paste0("Surv(Followup_time, ", event, ") ~ 1"))
  formula_clin_CF <- formula(paste0("Surv(Followup_time, ", event, ") ~ ", paste(c("Gender", "Age", covariates), collapse = "+")))

  dat_sens <- dat_perCpG[,c("ID","Amplicon","delta_AIC")]
  dat_sens$n_unadj <- NA
  dat_sens$delta_AIC_unadj <- NA
  dat_sens$n_clin <- NA
  dat_sens$delta_AIC_clin <- NA
  dat_sens$rho_age <- NA
  dat_sens$pval_age <- NA
  dat_sens$pval_sex <- NA
  for (i in 1:n_cpg){
    cn_meth_i <- paste("Meth_CpG_", i, sep = "")

    # Unadjusted: Meth vs empty model
    mask_unadj_i <- complete.cases(dat[,c(cn_meth_i, cn_time_event)])
    cox_unadj_i <- coxph(formula(paste0("Surv(Followup_time, ", event, ") ~ ", cn_meth_i)), data = dat[mask_unadj_i,])
    cox_null_i <- coxph(formula_null, data = dat[mask_unadj_i,])
    dat_sens$n_unadj[i] <- sum(mask_unadj_i)
    dat_sens$delta_AIC_unadj[i] <- AIC(cox_unadj_i) - AIC(cox_null_i)

    # Clinical covariates in both models
    mask_clin_i <- complete.cases(dat[,c(cn_meth_i, cn_time_event, "Gender", "Age", covariates)])
    formula_clin_i <- formula(paste0("Surv(Followup_time, ", event, ") ~ ", paste(c(cn_meth_i, "Gender", "Age", covariates), collapse = "+")))
    cox_clin_i <- coxph(formula_clin_i, data = dat[mask_clin_i,])
    cox_clin_CF_i <- coxph(formula_clin_CF, data = dat[mask_clin_i,])
    dat_sens$n_clin[i] <- sum(mask_clin_i)
    dat_sens$delta_AIC_clin[i] <- AIC(cox_clin_i) - AIC(cox_clin_CF_i)

    # Confounding check: methylation vs age and sex
    mask_cov_i <- complete.cases(dat[,c(cn_meth_i, "Gender", "Age")])
    cor_age_i <- cor.test(dat[mask_cov_i, cn_meth_i], dat$Age[mask_cov_i], method = "spearman", exact = F)
    dat_sens$rho_age[i] <- cor_age_i$estimate
    dat_sens$pval_age[i] <- cor_age_i$p.value
    dat_sens$pval_sex[i] <- wilcox.test(dat[mask_cov_i, cn_meth_i] ~ dat$Gender[mask_cov_i], exact = F)$p.value
  }
  write.table(dat_sens, paste0(gene, "_AIC_table_sensitivity.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

  print(sprintf("%s CpGs Delta AIC < -6: %d main (non-nested vs age + sex), %d unadjusted, %d clinical (of %d)", gene,
                sum(dat_sens$delta_AIC < -6, na.rm = T), sum(dat_sens$delta_AIC_unadj < -6, na.rm = T),
                sum(dat_sens$delta_AIC_clin < -6, na.rm = T), nrow(dat_sens)))
  print(sprintf("%s min Delta AIC: %.2f main, %.2f unadjusted, %.2f clinical; clinical models on %d-%d samples", gene,
                min(dat_sens$delta_AIC, na.rm = T), min(dat_sens$delta_AIC_unadj, na.rm = T),
                min(dat_sens$delta_AIC_clin, na.rm = T), min(dat_sens$n_clin), max(dat_sens$n_clin)))
  print(sprintf("%s CpGs correlated with age (p < 0.05): %d, median rho %.2f; differing by sex (p < 0.05): %d", gene,
                sum(dat_sens$pval_age < 0.05), median(dat_sens$rho_age), sum(dat_sens$pval_sex < 0.05)))

  # Per amplicon (as Figure 4: mean and minimum over the target-sequence CpGs of the amplicon)
  mask_target_sens <- !(dat_sens$ID %in% CpGs_in_primer)
  avAIC_sens <- data.frame(Amplicon = unique(dat_sens$Amplicon), stringsAsFactors = F)
  avAIC_sens$avAIC <- sapply(avAIC_sens$Amplicon, function(x){mean(dat_sens$delta_AIC[dat_sens$Amplicon==x & mask_target_sens])})
  avAIC_sens$avAIC_unadj <- sapply(avAIC_sens$Amplicon, function(x){mean(dat_sens$delta_AIC_unadj[dat_sens$Amplicon==x & mask_target_sens])})
  avAIC_sens$avAIC_clin <- sapply(avAIC_sens$Amplicon, function(x){mean(dat_sens$delta_AIC_clin[dat_sens$Amplicon==x & mask_target_sens])})
  avAIC_sens$minAIC <- sapply(avAIC_sens$Amplicon, function(x){min(dat_sens$delta_AIC[dat_sens$Amplicon==x & mask_target_sens])})
  avAIC_sens$minAIC_unadj <- sapply(avAIC_sens$Amplicon, function(x){min(dat_sens$delta_AIC_unadj[dat_sens$Amplicon==x & mask_target_sens])})
  avAIC_sens$minAIC_clin <- sapply(avAIC_sens$Amplicon, function(x){min(dat_sens$delta_AIC_clin[dat_sens$Amplicon==x & mask_target_sens])})
  print(avAIC_sens)
  write.table(avAIC_sens, paste0(gene, "_avAIC_sensitivity.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

  # MSP: unadjusted and clinical
  mask_msp_unadj <- complete.cases(dat[,c("Meth", cn_time_event)])
  cox_msp_unadj <- coxph(formula(paste0("Surv(Followup_time, ", event, ") ~ Meth")), data = dat[mask_msp_unadj,])
  cox_msp_null <- coxph(formula_null, data = dat[mask_msp_unadj,])
  mask_msp_clin <- complete.cases(dat[,c("Meth", cn_time_event, "Gender", "Age", covariates)])
  cox_msp_clin <- coxph(formula(paste0("Surv(Followup_time, ", event, ") ~ ", paste(c("Meth", "Gender", "Age", covariates), collapse = "+"))),
                        data = dat[mask_msp_clin,])
  cox_msp_clin_CF <- coxph(formula_clin_CF, data = dat[mask_msp_clin,])
  print(sprintf("%s MSP Delta AIC: unadjusted %.2f (n = %d), clinical %.2f (n = %d)", gene,
                AIC(cox_msp_unadj) - AIC(cox_msp_null), sum(mask_msp_unadj),
                AIC(cox_msp_clin) - AIC(cox_msp_clin_CF), sum(mask_msp_clin)))
  return(dat_sens)
}

# Clinical covariates: types and missing values
ccRCC_annot_sort_BS$Size <- as.numeric(ccRCC_annot_sort_BS$Size)
ccRCC_annot_sort_BS$Stage <- as.factor(ccRCC_annot_sort_BS$Stage)
print(sprintf("Melanoma Breslow: %d/%d missing", sum(is.na(melanoma_annot_sort_BS$Breslow)), nrow(melanoma_annot_sort_BS)))
print(table(ccRCC_annot_sort_BS$Stage, ccRCC_annot_sort_BS$NK_death, useNA = "ifany", dnn = c("Stage", "Death")))
print(sprintf("ccRCC Size: %d/%d missing", sum(is.na(ccRCC_annot_sort_BS$Size)), nrow(ccRCC_annot_sort_BS)))
# expect: few or no missing values; Stage with few levels (each level costs one parameter)

dat_LY75_sens <- sensitivity_cox(melanoma_annot_sort_BS, nrow(Mval_LY75_sort), dat_LY75_ROC_logit_perCpG, "MK_death", c("Breslow"), "LY75")
dat_GREM1_sens <- sensitivity_cox(ccRCC_annot_sort_BS, nrow(Mval_GREM1_sort), dat_GREM1_ROC_logit_perCpG, "NK_death", c("Stage","Size"), "GREM1")

##################################################################################
# Non-nested (competing models) LY75 + GREM1: methylation only vs clinical base #
##################################################################################
# Original design (GitHub): Surv ~ Meth vs Surv ~ Gender + Age. Delta AIC = AIC(Meth only) - AIC(base), same samples.
# Within a marker the base is (nearly) a constant shift, so the ranking of CpGs equals the unadjusted ranking.

nonnested_cox <- function(dat, n_cpg, dat_sens, event, bases, gene){
  cn_time_event <- c("Followup_time", event)
  dat_nn <- dat_sens[,c("ID","Amplicon","delta_AIC","delta_AIC_unadj")]
  for (base_name in names(bases)){
    formula_base <- formula(paste0("Surv(Followup_time, ", event, ") ~ ", paste(bases[[base_name]], collapse = "+")))
    cn_delta <- paste0("delta_AIC_nn_", base_name)
    cn_n <- paste0("n_nn_", base_name)
    dat_nn[,cn_delta] <- NA
    dat_nn[,cn_n] <- NA
    for (i in 1:n_cpg){
      cn_meth_i <- paste("Meth_CpG_", i, sep = "")
      mask_nn_i <- complete.cases(dat[,c(cn_meth_i, cn_time_event, bases[[base_name]])])
      cox_meth_only_i <- coxph(formula(paste0("Surv(Followup_time, ", event, ") ~ ", cn_meth_i)), data = dat[mask_nn_i,])
      cox_base_i <- coxph(formula_base, data = dat[mask_nn_i,])
      dat_nn[i,cn_n] <- sum(mask_nn_i)
      dat_nn[i,cn_delta] <- AIC(cox_meth_only_i) - AIC(cox_base_i)
    }
    print(sprintf("%s non-nested vs %s: CpGs Delta AIC < -6: %d/%d, min %.2f, samples %d-%d", gene, base_name,
                  sum(dat_nn[,cn_delta] < -6, na.rm = T), nrow(dat_nn), min(dat_nn[,cn_delta], na.rm = T),
                  min(dat_nn[,cn_n]), max(dat_nn[,cn_n])))

    # MSP: methylation status only vs base
    mask_msp_nn <- complete.cases(dat[,c("Meth", cn_time_event, bases[[base_name]])])
    cox_msp_only <- coxph(formula(paste0("Surv(Followup_time, ", event, ") ~ Meth")), data = dat[mask_msp_nn,])
    cox_msp_base <- coxph(formula_base, data = dat[mask_msp_nn,])
    print(sprintf("%s MSP non-nested vs %s: Delta AIC %.2f (n = %d)", gene, base_name,
                  AIC(cox_msp_only) - AIC(cox_msp_base), sum(mask_msp_nn)))
  }
  write.table(dat_nn, paste0(gene, "_AIC_table_nonnested.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

  # Per amplicon (as Figure 4: mean and minimum over the target-sequence CpGs of the amplicon)
  mask_target_nn <- !(dat_nn$ID %in% CpGs_in_primer)
  avAIC_nn <- data.frame(Amplicon = unique(dat_nn$Amplicon), stringsAsFactors = F)
  for (cn_j in grep("^delta_AIC", colnames(dat_nn), value = T)){
    avAIC_nn[,sub("^delta_AIC", "avAIC", cn_j)] <- sapply(avAIC_nn$Amplicon, function(x){mean(dat_nn[dat_nn$Amplicon==x & mask_target_nn, cn_j])})
    avAIC_nn[,sub("^delta_AIC", "minAIC", cn_j)] <- sapply(avAIC_nn$Amplicon, function(x){min(dat_nn[dat_nn$Amplicon==x & mask_target_nn, cn_j])})
  }
  print(avAIC_nn)
  write.table(avAIC_nn, paste0(gene, "_avAIC_nonnested.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
  return(dat_nn)
}

dat_LY75_nn <- nonnested_cox(melanoma_annot_sort_BS, nrow(Mval_LY75_sort), dat_LY75_sens, "MK_death",
                             list(agesex = c("Gender","Age"), clin = c("Gender","Age","Breslow")), "LY75")
dat_GREM1_nn <- nonnested_cox(ccRCC_annot_sort_BS, nrow(Mval_GREM1_sort), dat_GREM1_sens, "NK_death",
                              list(agesex = c("Gender","Age"), clin = c("Gender","Age","Stage","Size")), "GREM1")

# Check (v12): the main analysis is now the non-nested comparison with age + sex
print(sprintf("Max |main - non-nested age + sex|: LY75 %.4f, GREM1 %.4f",
              max(abs(dat_LY75_nn$delta_AIC - dat_LY75_nn$delta_AIC_nn_agesex)),
              max(abs(dat_GREM1_nn$delta_AIC - dat_GREM1_nn$delta_AIC_nn_agesex))))
# expect: 0


#####################################################################
# Sensitivity LY75 + GREM1: read-based models, non-nested (ToDo A14) #
#####################################################################
# Height + density only vs base model (age + sex, or with the clinical covariates), on the samples complete for both.
# Was part of the read-based loops of the main script; rbAIC_nn_agesex equals rbAIC_hybr of avAIC_table_ext.txt.

## Define settings
#####################
ccRCC_annot_sort$Size <- as.numeric(ccRCC_annot_sort$Size)
ccRCC_annot_sort$Stage <- as.factor(ccRCC_annot_sort$Stage)
rb_info <- list(LY75 = list(annot = melanoma_annot_sort, sample_col = "melanoma", skip = "LY75_6_S", clin = "Breslow"),
                GREM1 = list(annot = ccRCC_annot_sort, sample_col = "ccRCC", skip = "GREM1_60", clin = c("Stage","Size")))

## Models
#####################
avAIC_rb_nn <- avAIC_table_ext[,c("AmpliconID","rbAIC_hybr")]
avAIC_rb_nn$rbAIC_nn_agesex <- NA
avAIC_rb_nn$rbAIC_nn_clin <- NA
for (gene in names(rb_info)){
  info <- rb_info[[gene]]
  annot_rb <- info$annot
  for (amplicon in amplicons[grepl(gene, amplicons)]){
    if (amplicon == info$skip){
      next()
    }
    annot_rb$meth_state <- NA
    annot_rb$meth_state2 <- NA
    for (i in 1:nrow(annot_rb)){
      sample_i <- sample_annot$Sample_ID[sample_annot[[info$sample_col]]==annot_rb$Sample_ID[i]]
      annot_rb$meth_state[i] <- meth_state[sample_i,amplicon]
      annot_rb$meth_state2[i] <- meth_state2[sample_i,amplicon]
    }
    for (base_name in c("agesex", "clin")){
      cn_base <- c("Gender","Age")
      if (base_name == "clin"){
        cn_base <- c("Gender","Age",info$clin)
      }
      mask_nn_rb <- complete.cases(annot_rb[,c("meth_state","meth_state2","Followup_time","NK_death",cn_base)])
      cox_rb_only <- coxph(formula("Surv(Followup_time, NK_death) ~ meth_state + meth_state2"), data = annot_rb[mask_nn_rb,])
      cox_rb_base <- coxph(formula(paste0("Surv(Followup_time, NK_death) ~ ", paste(cn_base, collapse = "+"))), data = annot_rb[mask_nn_rb,])
      avAIC_rb_nn[avAIC_rb_nn$AmpliconID==amplicon, paste0("rbAIC_nn_", base_name)] <- AIC(cox_rb_only) - AIC(cox_rb_base)
    }
  }
}
print(avAIC_rb_nn[!is.na(avAIC_rb_nn$rbAIC_nn_agesex),])
print(sprintf("Max |rbAIC_nn_agesex - rbAIC_hybr|: %.4f",
              max(abs(avAIC_rb_nn$rbAIC_nn_agesex - as.numeric(avAIC_rb_nn$rbAIC_hybr)), na.rm = T)))
# expect: 0
write.table(avAIC_rb_nn,"avAIC_rb_nonnested.txt",col.names = T,row.names = F,quote = F, sep="\t")

###############################################################
# Sensitivity primer CpGs: rank-based tests, NDRG4 floor (R9) #
###############################################################

## Sensitivity: rank-based tests and without NDRG4 CpGs at the Delta AIC floor (R9, ToDo A1)
#####################

mask_floor_correl <- dat_correl$ID %in% dat_NDRG4_ROC_logit_perCpG$ID[dat_NDRG4_ROC_logit_perCpG$at_floor]
print(sprintf("CpGs at the NDRG4 floor: %d/%d", sum(mask_floor_correl), nrow(dat_correl)))
dat_correl_nofloor <- dat_correl[!mask_floor_correl, ]
dat_correl_nofloor_not_in_primer <- dat_correl_nofloor[(dat_correl_nofloor$dist != 0) | is.na(dat_correl_nofloor$dist),]
dat_correl_nofloor_filt <- dat_correl_nofloor[(dat_correl_nofloor$CG_qual!=0) & (dat_correl_nofloor$dist!=0) & (!is.na(dat_correl_nofloor$dist)),]

### 0 vs 1 / 0 vs 2 primers with CpGs (manuscript: t-test p = 2.4E-6 and < 2.2E-16)
print(wilcox.test(dat_correl_not_in_primer$deltaAIC[dat_correl_not_in_primer$CG_qual==0],
                  dat_correl_not_in_primer$deltaAIC[dat_correl_not_in_primer$CG_qual==1], exact = F))
print(wilcox.test(dat_correl_not_in_primer$deltaAIC[dat_correl_not_in_primer$CG_qual==0],
                  dat_correl_not_in_primer$deltaAIC[dat_correl_not_in_primer$CG_qual==2], exact = F))
print(t.test(dat_correl_nofloor_not_in_primer$deltaAIC[dat_correl_nofloor_not_in_primer$CG_qual==0],
             dat_correl_nofloor_not_in_primer$deltaAIC[dat_correl_nofloor_not_in_primer$CG_qual==2]))

### Correlations (manuscript: rho = -0.275 and 0.287)
print(cor.test(dat_correl_filt$diff_meth, dat_correl_filt$deltaAIC, method = "spearman", exact = F))
print(cor.test(dat_correl_filt$dist_CpG, dat_correl_filt$deltaAIC, method = "spearman", exact = F))
print(cor.test(dat_correl_nofloor_filt$diff_meth, dat_correl_nofloor_filt$deltaAIC))
print(cor.test(dat_correl_nofloor_filt$dist_CpG, dat_correl_nofloor_filt$deltaAIC))



##################################################################
# Intra-amplicon differences: original permutation test (E18)    #
##################################################################
# Test of the original submission (values of each CpG permuted between samples; null: neither CpG is informative), kept
# for comparison with the permutation test (within samples) of the main analysis (rebuttal, Supplementary Figure 4). Plots in AIC_plots/ and AIC_plots_red/.

source(file.path(project_dir,"Intra_amplicon_diff_functions.R"))


# Create directories
if(!dir.exists("AIC_plots")){dir.create("AIC_plots")}
if(!dir.exists("AIC_plots/CpG_plots")){dir.create("AIC_plots/CpG_plots")}

# all comparisons
cex_AIC_plots <- 1.25
seed_perm <- 42
fdr_treshold_toplot <- 0.10
n_bp_add <- -1

## including primer CpGs
####################

#NDRG4
pval_list_NDRG4 <- perform_AIC_analysis(primer_annot_NDRG4,Mval_NDRG4_sort,CRC_annot_sort,"NDRG4")

#GREM1
pval_list_GREM1 <- perform_AIC_analysis(primer_annot_GREM1,Mval_GREM1_sort,ccRCC_annot_sort,"GREM1")

#LY75
pval_list_LY75 <- perform_AIC_analysis(primer_annot_LY75,Mval_LY75_sort,melanoma_annot_sort,"LY75")

# Full table
part_pval_GREM1 <- cbind(rep(names(pval_list_GREM1)[pval_list_GREM1!="No CpGs"],times=as.numeric(unlist(lapply(pval_list_GREM1,nrow)))),
                         do.call(rbind,pval_list_GREM1[pval_list_GREM1!="No CpGs"]))
colnames(part_pval_GREM1) <- c("Amplicon",colnames(part_pval_GREM1[2:ncol(part_pval_GREM1)]))
rownames(part_pval_GREM1) <- 1:nrow(part_pval_GREM1)
part_pval_NDRG4 <- cbind(rep(names(pval_list_NDRG4)[pval_list_NDRG4!="No CpGs"],times=as.numeric(unlist(lapply(pval_list_NDRG4,nrow)))),
                         do.call(rbind,pval_list_NDRG4[pval_list_NDRG4!="No CpGs"]))
colnames(part_pval_NDRG4) <- c("Amplicon",colnames(part_pval_NDRG4[2:ncol(part_pval_NDRG4)]))
rownames(part_pval_NDRG4) <- 1:nrow(part_pval_NDRG4)
part_pval_LY75 <- cbind(rep(names(pval_list_LY75)[pval_list_LY75!="No CpGs"],times=as.numeric(unlist(lapply(pval_list_LY75,nrow)))),
                        do.call(rbind,pval_list_LY75[pval_list_LY75!="No CpGs"]))
colnames(part_pval_LY75) <- c("Amplicon",colnames(part_pval_LY75[2:ncol(part_pval_LY75)]))
rownames(part_pval_LY75) <- 1:nrow(part_pval_LY75)
                
pval_total_list <- rbind(part_pval_GREM1,part_pval_NDRG4,part_pval_LY75)
rownames(pval_total_list) <- 1:nrow(pval_total_list)
write.table(pval_total_list,file="pval_total_list.txt",col.names=T,row.names=F,quote=F,sep="\t")

unique(pval_total_list$Amplicon[pval_total_list$fdr<0.10])
unique(pval_total_list$Amplicon[pval_total_list$pval<0.10])


# # checks
# dat_check_GREM1 <- dat_GREM1_ROC_logit_perCpG[,c("ID","AIC")]
# dat_check_GREM1_pvals <- do.call(rbind,pval_list_GREM1)[!duplicated(do.call(rbind,pval_list_GREM1)$CpG_1),c("CpG_1","AIC_1","pval","fdr")]
# dat_check_GREM1 <- cbind(dat_check_GREM1,matrix(rep(NA,4*nrow(dat_check_GREM1)),nrow=nrow(dat_check_GREM1)))
# dat_check_GREM1[dat_check_GREM1$ID%in%dat_check_GREM1_pvals$CpG_1,3:6] <- dat_check_GREM1_pvals 
# daT_check_GREM1
# 
# dat_check_NDRG4 <- dat_NDRG4_ROC_logit_perCpG[,c("ID","AIC")]
# dat_check_NDRG4_pvals <- do.call(rbind,pval_list_NDRG4)[!duplicated(do.call(rbind,pval_list_NDRG4)$CpG_1),c("CpG_1","AIC_1","pval","fdr")]
# dat_check_NDRG4 <- cbind(dat_check_NDRG4,matrix(rep(NA,4*nrow(dat_check_NDRG4)),nrow=nrow(dat_check_NDRG4)))
# dat_check_NDRG4[dat_check_NDRG4$ID%in%dat_check_NDRG4_pvals$CpG_1,3:6] <- dat_check_NDRG4_pvals 
# dat_check_NDRG4
# 
# dat_check_LY75 <- dat_LY75_ROC_logit_perCpG[,c("ID","AIC")]
# dat_check_LY75_pvals <- do.call(rbind,pval_list_LY75)[!duplicated(do.call(rbind,pval_list_LY75)$CpG_1),c("CpG_1","AIC_1","pval","fdr")]
# dat_check_LY75 <- cbind(dat_check_LY75,matrix(rep(NA,4*nrow(dat_check_LY75)),nrow=nrow(dat_check_LY75)))
# dat_check_LY75[dat_check_LY75$ID%in%dat_check_LY75_pvals$CpG_1,3:6] <- dat_check_LY75_pvals 
# dat_check_LY75

## excluding primer CpGs
####################

# unlink("AIC_plots_red/*")
if(!dir.exists("AIC_plots_red")){dir.create("AIC_plots_red")}

pval_cutoff_red <- 0.05

reduced_AIC_analysis("GREM1",pval_list_GREM1,pval_cutoff_red,dat_GREM1_ROC_logit_perCpG)

reduced_AIC_analysis("NDRG4",pval_list_NDRG4,pval_cutoff_red,dat_NDRG4_ROC_logit_perCpG)

pval_cutoff_red <- 0.10

reduced_AIC_analysis("LY75",pval_list_LY75,pval_cutoff_red,dat_LY75_ROC_logit_perCpG)


## Comparison with the permutation test within samples (main analysis)
#####################
key_perm <- paste(pval_total_list$CpG_1, pval_total_list$CpG_2, sep = "|")
key_swap <- paste(res_swap$CpG_1, res_swap$CpG_2, sep = "|")
if (sum(key_swap%in%key_perm) != length(key_swap)){
  stop("Pairs of the main analysis missing in the original permutation test")
}
res_swap_perm <- res_swap
res_swap_perm$delta_AIC_test <- as.numeric(setNames(pval_total_list$delta_AIC, key_perm)[key_swap])
res_swap_perm$pval_perm <- as.numeric(setNames(pval_total_list$pval, key_perm)[key_swap])
print(sprintf("Max |Delta AIC main test - original test| (z-scoring must not change the AIC): %.4f",
              max(abs(res_swap_perm$delta_AIC_obs - res_swap_perm$delta_AIC_test))))
# expect: ~0
write.table(res_swap_perm, "intra_amplicon_swap_vs_permutation.txt", col.names = T, row.names = F, quote = F, sep = "\t")
summary_swap_perm <- summary_swap
summary_swap_perm$min_p_perm <- sapply(summary_swap_perm$Amplicon, function(x){min(res_swap_perm$pval_perm[res_swap_perm$Amplicon == x])})
summary_swap_perm$n_sig_perm <- sapply(summary_swap_perm$Amplicon, function(x){sum(res_swap_perm$pval_perm[res_swap_perm$Amplicon == x] < res_swap_perm$cutoff[res_swap_perm$Amplicon == x])})
summary_swap_perm$star_perm <- summary_swap_perm$n_sig_perm > 0
print(summary_swap_perm)
write.table(summary_swap_perm, "intra_amplicon_swap_vs_permutation_summary.txt", col.names = T, row.names = F, quote = F, sep = "\t")
for (gene in unique(summary_swap_perm$Gene)){
  x <- summary_swap_perm[summary_swap_perm$Gene == gene,]
  print(sprintf("%s: amplicons with a star - original test %d, main test %d (of %d); pairs significant - original %d, main %d",
                gene, sum(x$star_perm), sum(x$star_swap), nrow(x), sum(x$n_sig_perm), sum(x$n_sig_swap)))
}

########################################################################
# Sensitivity LY75: intra-amplicon permutation test with sequencing run #
########################################################################

## LY75: permutation test (within samples) with sequencing run in both CpG models (R3.M6a, Louis 2026-11-17)
#####################
# Surv ~ CpG_a + Run vs Surv ~ CpG_b + Run, same target-CpG pairs, permutations, seed and cut-off as the main test.
# Decision rule (Louis): if the LY75 intra-amplicon results do not change, the main analysis stays without run.
RNGkind("L'Ecuyer-CMRG")
set.seed(42)
info_LY75_run <- pair_info$LY75
info_LY75_run$cn <- c(info_LY75_run$cn, "Run")
pairs_LY75 <- res_swap[res_swap$Gene == "LY75",]
print(sprintf("LY75: permutation test with Run for %d target-CpG pairs", nrow(pairs_LY75)))
res_swap_LY75_run <- NULL
for (r in 1:nrow(pairs_LY75)){
  dat_r <- build_pair_data(info_LY75_run, "LY75", pairs_LY75$CpG_1[r], pairs_LY75$CpG_2[r])
  dat_r$CpG_a <- as.numeric(scale(dat_r$CpG_a))
  dat_r$CpG_b <- as.numeric(scale(dat_r$CpG_b))
  obs_r <- aic_pair_diff(dat_r, "LY75", covariates = "+Run")
  null_r <- unlist(mclapply(1:n_swap, swap_pair, dat = dat_r, gene = "LY75", covariates = "+Run", mc.cores = 40))
  mask_ok <- !is.na(null_r)
  pval_swap_run <- mean(abs(null_r[mask_ok]) >= abs(obs_r) - 1e-9)
  res_swap_LY75_run <- rbind(res_swap_LY75_run, data.frame(Amplicon = pairs_LY75$Amplicon[r], CpG_1 = pairs_LY75$CpG_1[r],
                                                           CpG_2 = pairs_LY75$CpG_2[r], n_samples = nrow(dat_r),
                                                           delta_AIC_obs = pairs_LY75$delta_AIC_obs[r], delta_AIC_obs_run = obs_r,
                                                           pval_swap = pairs_LY75$pval_swap[r], pval_swap_run = pval_swap_run,
                                                           n_swap_ok_run = sum(mask_ok), stringsAsFactors = F))
}
RNGkind(rng_kind_ori)
res_swap_LY75_run$fdr_swap_amplicon_run <- ave(res_swap_LY75_run$pval_swap_run, res_swap_LY75_run$Amplicon,
                                               FUN = function(p){p.adjust(p, method = "BH")})
write.table(res_swap_LY75_run, "intra_amplicon_swap_test_LY75_run.txt", col.names = T, row.names = F, quote = F, sep = "\t")

# Summary: main vs with Run (counts, chance level, BH within amplicon, amplicons with a star, agreement of p-values)
n_pairs_LY75 <- nrow(res_swap_LY75_run)
for (cutoff in c(0.05, 0.10)){
  mask_main <- res_swap_LY75_run$pval_swap < cutoff
  mask_run <- res_swap_LY75_run$pval_swap_run < cutoff
  print(sprintf("LY75 pairs p < %.2f: main %d, with Run %d (of %d; %.0f expected by chance); in both %d",
                cutoff, sum(mask_main), sum(mask_run), n_pairs_LY75, n_pairs_LY75 * cutoff, sum(mask_main & mask_run)))
  print(sprintf("LY75 amplicons with a pair p < %.2f: main %s | with Run %s", cutoff,
                paste(sort(unique(res_swap_LY75_run$Amplicon[mask_main])), collapse = ", "),
                paste(sort(unique(res_swap_LY75_run$Amplicon[mask_run])), collapse = ", ")))
}
print(sprintf("LY75 min BH within amplicon: main %.3f, with Run %.3f; pairs BH < 0.05 with Run: %d",
              min(pairs_LY75$fdr_swap_amplicon), min(res_swap_LY75_run$fdr_swap_amplicon_run),
              sum(res_swap_LY75_run$fdr_swap_amplicon_run < 0.05)))
print(cor.test(res_swap_LY75_run$pval_swap, res_swap_LY75_run$pval_swap_run, method = "spearman", exact = F))
print(sprintf("LY75 samples per pair with Run: %d-%d",
              min(res_swap_LY75_run$n_samples), max(res_swap_LY75_run$n_samples)))
# expect: main 6 pairs at 0.05 / 21 at 0.10 (v18), no BH < 0.05


##############################################
# Simulation of the intra-amplicon tests (E18) #
##############################################
# Supplementary Figure 4: original permutation test (between samples) vs revised permutation test (within samples). Run last: the script sets its own seed
# and RNG kind. Output in simulation_intra_amplicon/ (~10-20 min on 40 cores).
source(file.path(project_dir,"Simulation_intra_amplicon_tests.R"))

#############
# Close log #
#############

print(sprintf("End: %s", format(Sys.time())))
print(sessionInfo())
sink(type = "message")
sink()
close(log_con)
