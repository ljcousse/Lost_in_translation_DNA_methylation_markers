#################################################
#################################################
##                                             ##
##  Simulation_intra_amplicon_tests.R          ##
##                                             ##
#################################################
#################################################
# Supplementary Figure (ToDo E18, R-E18 in the rebuttal): behaviour of the two tests for a difference in clinical value
# between two CpGs of the same amplicon, in the prognostic setting (Cox regression, n = 26 as in the ccRCC and melanoma
# cohorts).
#  - Permutation test (original submission): the observed Delta AIC = AIC(CpG 2) - AIC(CpG 1) is compared with the Delta
#    AIC after permuting the values of each CpG independently (null = neither CpG is informative).
#  - Swap test (revision): both CpGs z-scored, per sample the two values swapped with probability 0.5 (null = both CpGs
#    equally informative); p = fraction of swapped data sets with |Delta AIC| >= observed |Delta AIC|.
# Data: a latent methylation factor drives survival (log hazard ratio beta per SD); each CpG = r * latent + noise, so a
# CpG with a higher r carries more clinical value. Equal r = equal clinical value.
# Outputs (PROJECT_DIR/simulation_intra_amplicon/): simulation_results.txt (rate of p < 0.05 per setting and test),
# simulation_example_null.txt, Supplementary_Figure_simulation.pdf.

##################
# Load libraries #
##################
library(survival)
library(parallel)
library(ggplot2)
library(patchwork)

##################
# Initialisation #
##################

## Define settings
#####################
config_file <- Sys.getenv("PROJECT_CONFIG", "config.env")
if (!file.exists(config_file)){
  stop("Config file not found: ", config_file, ". Copy config.example.env to config.env and edit the paths.")
}
readRenviron(config_file)
project_dir <- Sys.getenv("PROJECT_DIR")
if (!dir.exists(project_dir)){
  stop("PROJECT_DIR does not exist: ", project_dir)
}

outdir <- file.path(project_dir, "simulation_intra_amplicon")
dir.create(outdir, showWarnings = F, recursive = T)
outfile_results <- file.path(outdir, "simulation_results.txt")
outfile_example <- file.path(outdir, "simulation_example_null.txt")
outfile_figure <- file.path(outdir, "Supplementary_Figure_simulation.pdf")

n_samples <- 26        # as the ccRCC and melanoma cohorts
n_rep <- 500           # simulated data sets per setting
n_null <- 1000         # permutations / swaps per data set
p_cutoff <- 0.05
n_cores <- 40          # set to 1 on Windows (mclapply forks)
seed_sim <- 42

# Settings: equal clinical value (false-positive rate) and unequal clinical value (power)
settings <- rbind(data.frame(scenario = "Equal clinical value", beta = c(0, 0.5, 1, 1.5), r1 = 0.85, r2 = 0.85),
                  data.frame(scenario = "Different clinical value", beta = 1.5, r1 = 0.95, r2 = c(0.85, 0.70, 0.50)))
settings$setting <- 1:nrow(settings)
print(settings)

## Functions
#####################
simulate_data <- function(beta, r1, r2){
  latent <- rnorm(n_samples)
  cpg_1 <- r1 * latent + sqrt(1 - r1^2) * rnorm(n_samples)
  cpg_2 <- r2 * latent + sqrt(1 - r2^2) * rnorm(n_samples)
  time_event <- rexp(n_samples, exp(beta * latent) / 1000)
  time_censor <- runif(n_samples, 0, 2500)
  dat <- data.frame(time = pmin(time_event, time_censor), event = as.numeric(time_event <= time_censor),
                    cpg_1 = as.numeric(scale(cpg_1)), cpg_2 = as.numeric(scale(cpg_2)))
  return(dat)
}

get_aic <- function(dat, x){
  model <- suppressWarnings(coxph(Surv(dat$time, dat$event) ~ x))
  return(AIC(model))
}

get_delta_aic <- function(dat, x_1, x_2){
  # Delta AIC = AIC(model CpG 2) - AIC(model CpG 1), as in Intra_amplicon_diff_functions
  return(get_aic(dat, x_2) - get_aic(dat, x_1))
}

get_null_permutation <- function(dat){
  null_perm <- rep(NA, n_null)
  for (b in 1:n_null){
    null_perm[b] <- get_delta_aic(dat, sample(dat$cpg_1), sample(dat$cpg_2))
  }
  return(null_perm)
}

get_null_swap <- function(dat){
  null_swap <- rep(NA, n_null)
  for (b in 1:n_null){
    mask_swap <- runif(n_samples) < 0.5
    x_1 <- ifelse(mask_swap, dat$cpg_2, dat$cpg_1)
    x_2 <- ifelse(mask_swap, dat$cpg_1, dat$cpg_2)
    null_swap[b] <- get_delta_aic(dat, x_1, x_2)
  }
  return(null_swap)
}

run_one_data_set <- function(i, beta, r1, r2){
  dat <- simulate_data(beta, r1, r2)
  obs <- get_delta_aic(dat, dat$cpg_1, dat$cpg_2)
  null_perm <- get_null_permutation(dat)
  null_swap <- get_null_swap(dat)
  # two-sided p-values exactly as in the analysis scripts
  pval_perm <- 2 * min(mean(obs > null_perm), mean(obs < null_perm))
  pval_swap <- mean(abs(null_swap) >= abs(obs) - 1e-9)
  return(c(n_events = sum(dat$event), pval_perm = pval_perm, pval_swap = pval_swap))
}

##############
# Simulation #
##############
RNGkind("L'Ecuyer-CMRG")
set.seed(seed_sim)

results <- NULL
for (s in 1:nrow(settings)){
  print(sprintf("Setting %d/%d: %s, beta %.1f, r1 %.2f, r2 %.2f (%s)", s, nrow(settings), settings$scenario[s],
                settings$beta[s], settings$r1[s], settings$r2[s], format(Sys.time())))
  res_s <- mclapply(1:n_rep, run_one_data_set, beta = settings$beta[s], r1 = settings$r1[s], r2 = settings$r2[s],
                    mc.cores = n_cores)
  res_s <- do.call(rbind, res_s)
  for (test in c("pval_perm", "pval_swap")){
    n_sign <- sum(res_s[,test] < p_cutoff)
    ci <- binom.test(n_sign, n_rep)$conf.int
    results <- rbind(results, data.frame(settings[s,], test = test, rate_sign = n_sign / n_rep,
                                         ci_low = ci[1], ci_high = ci[2], mean_events = mean(res_s[,"n_events"])))
  }
}
results$test <- ifelse(results$test == "pval_perm", "Permutation test (original)", "Swap test (revision)")
print(results)
write.table(results, outfile_results, col.names = T, row.names = F, quote = F, sep = "\t")
# expect (earlier pilot, 100 data sets): equal value, beta 0: ~0.03 both; beta 1.0: perm ~0.26, swap ~0.04;
# beta 1.5: perm ~0.42, swap ~0.03; r 0.95 vs 0.50: swap ~0.28

## Example data set: null distributions
#####################
set.seed(seed_sim)
dat_example <- simulate_data(beta = 1.5, r1 = 0.85, r2 = 0.85)
obs_example <- get_delta_aic(dat_example, dat_example$cpg_1, dat_example$cpg_2)
dat_null <- rbind(data.frame(test = "Permutation test (original)", delta_AIC = get_null_permutation(dat_example)),
                  data.frame(test = "Swap test (revision)", delta_AIC = get_null_swap(dat_example)))
print(sprintf("Example (two equally informative CpGs, beta 1.5): observed Delta AIC %.2f, %d events",
              obs_example, sum(dat_example$event)))
write.table(dat_null, outfile_example, col.names = T, row.names = F, quote = F, sep = "\t")

########
# Plot #
########
colours_test <- c("Permutation test (original)" = "#eb6834", "Swap test (revision)" = "#2a78d6")

p1 <- ggplot(dat_null, aes(x = delta_AIC, fill = test)) +
  geom_histogram(bins = 60, alpha = 0.6, position = "identity") +
  geom_vline(xintercept = c(-1, 1) * abs(obs_example), linetype = "dashed") +
  scale_fill_manual(values = colours_test) +
  labs(x = expression(Delta*"AIC between the two CpGs"), y = "Number of permuted / swapped data sets", fill = NULL,
       title = "A. Null distributions, two equally informative CpGs") +
  theme_bw() +
  theme(legend.position = "top")

mask_equal <- results$scenario == "Equal clinical value"
p2 <- ggplot(results[mask_equal,], aes(x = beta, y = rate_sign, colour = test)) +
  geom_hline(yintercept = p_cutoff, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2.5, position = position_dodge(width = 0.1)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.1, position = position_dodge(width = 0.1)) +
  scale_colour_manual(values = colours_test) +
  labs(x = "Prognostic effect of both CpGs (log hazard ratio per SD of the latent factor)",
       y = sprintf("Fraction p < %.2f (false positives)", p_cutoff), colour = NULL,
       title = "B. Equal clinical value") +
  theme_bw() +
  theme(legend.position = "none")

mask_diff <- results$scenario == "Different clinical value"
p3 <- ggplot(results[mask_diff,], aes(x = r2, y = rate_sign, colour = test)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.03)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.03, position = position_dodge(width = 0.03)) +
  scale_colour_manual(values = colours_test) +
  labs(x = "Correlation of the second CpG with the latent factor (first CpG: 0.95)",
       y = sprintf("Fraction p < %.2f (power)", p_cutoff), colour = NULL,
       title = "C. Different clinical value (log hazard ratio 1.5)") +
  theme_bw() +
  theme(legend.position = "none")

cairo_pdf(outfile_figure, width = 8, height = 11) # cairo: Delta glyph (ToDo E20)
print(p1 + p2 + p3 + plot_layout(ncol = 1))
dev.off()
print(sprintf("Done: %s", format(Sys.time())))
