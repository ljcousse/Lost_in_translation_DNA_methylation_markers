# Lost in translation: DNA methylation markers

Code for the manuscript *"Lost in translation: How CpG site selection and assay design determine the clinical value of cancer DNA methylation markers"* (Lommen, Coussement et al., [[journal, year, DOI]]).

We profiled three established cancer DNA methylation markers at single-CpG resolution with targeted (amplicon) bisulfite sequencing: NDRG4 (colorectal cancer versus normal colon, diagnostic), GREM1 (clear-cell renal cell carcinoma, prognostic) and LY75 (melanoma, prognostic). For every CpG, amplicon, read-based measure and MSP assay, we quantified the clinical value with the Akaike Information Criterion (ΔAIC compared with a baseline model). CpGs within the same amplicon were compared with a permutation test in which the values of the two CpGs are permuted within samples, so that its null hypothesis is that both CpGs have equal clinical value.

## Repository contents

| File | Purpose |
|---|---|
| `TargetedBisulphiteSeq.sh` | Pre-processing: GRCh37 reference genome (Ensembl), adapter and quality trimming (Trim Galore), quality control (FastQC, MultiQC), bisulfite alignment (Bismark/Bowtie2), amplicon assignment (calls `AmpliconAwareMapping.R`), methylation calling (`bismark_methylation_extractor`, `bismark2bedGraph`) and read-level methylation patterns (wgbs_tools `bam2pat`) |
| `AmpliconAwareMapping.R` | Assigns every read pair to the amplicon it originates from, so that overlapping amplicons are analysed separately (called by `TargetedBisulphiteSeq.sh`) |
| `TargetedBisulphite_analysis.R` | Main analysis: M-values per CpG, ΔAIC per CpG, amplicon and MSP assay, CpGs covered by two amplicons, primer-CpG analyses (including local CpG density), intra-amplicon permutation test (within samples), read-based measures (relative co-methylation height and density) and the figures and tables of the manuscript. Saves the objects needed by the sensitivity analyses |
| `TargetedBisulphite_sensitivity.R` | Sensitivity analyses on the objects saved by the main analysis: NDRG4 without sample 178N, LY75 with sequencing run, GREM1/LY75 unadjusted and with clinical covariates, rank-based primer-CpG tests, the original permutation test (values permuted between samples; comparison), the LY75 permutation test with sequencing run and the simulation |
| `Intra_amplicon_diff_functions.R` | Functions for the intra-amplicon comparisons (original permutation test) and the amplicon plots of Figure 5 (sourced by both analysis scripts) |
| `Plots_AIC_function.R` | Plot function for the ΔAIC figures (Figure 2; sourced by the main analysis) |
| `Simulation_intra_amplicon_tests.R` | Simulation behind Supplementary Figure 4: false-positive rate and power of the original permutation test (between samples) and the revised permutation test (within samples) (sourced by the sensitivity script; can also be run on its own) |
| `conda_environment.yml` | Conda environment with the tools and package versions used |
| `config.example.env` | Template for the local paths (copy to `config.env`) |

## Installation

```bash
git clone https://github.com/ljcousse/Lost_in_translation_DNA_methylation_markers.git
cd Lost_in_translation_DNA_methylation_markers

conda env create -f conda_environment.yml
conda activate lost_in_translation
```

`wgbs_tools` is not available via conda. Install it from https://github.com/nloyfer/wgbs_tools into `PROJECT_DIR/wgbs_tools`; `TargetedBisulphiteSeq.sh` initialises the GRCh37 reference (`wgbstools init_genome ... hg37_87`).

The analysis was run on a Linux server (AlmaLinux 9) with R 4.2.3 and Bioconductor 3.16. The exact versions are in `conda_environment.yml`.

## Configuration

Local paths are not hard-coded. Copy the template and fill in your paths:

```bash
cp config.example.env config.env
```

`config.env` contains:
- `PROJECT_DIR`: working directory with the input files; all intermediate files and output are written here;
- `FASTQ_DIR_RUN231`, `FASTQ_DIR_RUN239`: FASTQ files of the two sequencing runs (pre-processing only).

All scripts read `config.env` from the directory they are started in; use another file with `PROJECT_CONFIG=/path/to/file.env`. Keep `config.env` local (it is listed in `.gitignore`).

The scripts expect the repository files in `PROJECT_DIR` (the pre-processing calls `AmpliconAwareMapping.R` there and the analysis scripts source their function files from there): clone the repository into `PROJECT_DIR`, or copy the scripts there, and start every step from `PROJECT_DIR`.

## Running the analysis

1. **Pre-processing:** `bash TargetedBisulphiteSeq.sh` (FASTQ → trimmed reads → Bismark alignments → reads assigned per amplicon → CpG methylation calls in `meth/` and wgbs_tools `.pat` files in `output_wgbstools/`).
2. **Main analysis:** `Rscript TargetedBisulphite_analysis.R`. Output goes to `PROJECT_DIR`; a time-stamped log (`TargetedBisulphite_analysis_<date>_<time>.log`) records the console output and session information. At the end, the objects needed by the sensitivity analyses are saved in `output_for_sensitivity_analysis/objects_for_sensitivity_analysis.Rda`.
3. **Sensitivity analyses:** `Rscript TargetedBisulphite_sensitivity.R` (needs `output_for_sensitivity_analysis/objects_for_sensitivity_analysis.Rda`; own log `TargetedBisulphite_sensitivity_<date>_<time>.log`). The simulation writes to `simulation_intra_amplicon/`.

Both permutation tests and the simulation use `parallel::mclapply` with 40 cores, the amplicon assignment with 20; change `mc.cores` (simulation: `n_cores`) for smaller machines. `mclapply` forks, so on Windows use 1 core. The original permutation test (sensitivity script) is the slowest step.

### Main outputs

| Output | Content |
|---|---|
| `NDRG4_AIC.pdf`, `GREM1_AIC.pdf`, `LY75_AIC.pdf`, `LegendFirstPanel.jpg` | ΔAIC per CpG along the gene region (Figure 2) |
| `waterfall_plots/` | MESA plots (Figure 3, Supplementary Information 1) |
| `AIC_plots/Figure4_DeltaAIC_per_amplicon.pdf`, `AIC_plots/Figure4_data.txt` | ΔAIC per amplicon (Figure 4) |
| `AIC_plots_swap/` | Intra-amplicon differences, permutation test (Figure 5) |
| `intra_amplicon_swap_test.txt`, `intra_amplicon_swap_summary.txt` | Permutation test (within samples) per pair of target-sequence CpGs and per amplicon |
| `supplementary_table_1.txt` | Primer and amplicon annotation with primer-CpG methylation (Supplementary Table 1) |
| `res_dupl_CpG.txt` | CpGs covered by two amplicons, paired t-test on M-values with Benjamini-Hochberg correction (Supplementary Table 3) |
| `NDRG4_AIC_table.txt`, `GREM1_AIC_table.txt`, `LY75_AIC_table.txt` | ΔAIC per CpG (Supplementary Table 4) |
| `*_AIC_table_nonnested.txt`, `*_avAIC_nonnested.txt`, `LY75_AIC_table_run.txt`, `LY75_avAIC_run.txt`, `intra_amplicon_swap_test_LY75_run.txt` | Sensitivity analyses (Supplementary Table 5; sensitivity script) |
| `avAIC_table_ext.txt` | Mean and minimum CpG-based ΔAIC and read-based ΔAIC per amplicon (Supplementary Table 6) |
| `violin_deltaAIC_CpG_*.jpg`, `AIC_plots/cpg_density_profile.pdf` | Primer-CpG analyses (Supplementary Figure 1) and CpG density per gene region |
| `suppl_fig_3/` | Supplementary Figure 3 |
| `simulation_intra_amplicon/` | Simulation of the intra-amplicon tests (Supplementary Figure 4; sensitivity script) |

### Statistical design (short)

- **NDRG4 (diagnostic, paired tumour/normal):** logistic regression. Tissue ~ CpG methylation (M-value) + patient ID, compared with tissue ~ patient ID.
- **GREM1 and LY75 (prognostic):** Cox regression. Methylation alone is compared with a baseline model with age and sex (competing models, same samples). A negative ΔAIC means that the methylation model explains the outcome better than the baseline model.
- **Intra-amplicon comparisons (permutation test within samples):** for every pair of target-sequence CpGs in an amplicon, the ΔAIC between the two CpG models (as above) is compared with its distribution when the two M-values (z-scored) are exchanged per sample with probability 0.5 (2,000 permutations). Null hypothesis: both CpGs have equal clinical value; amplicons with a pair at p < 0.05 are marked in Figures 4 and 5.
- **Sensitivity analyses:** methylation versus an empty model and versus age, sex and clinical covariates (Breslow thickness for melanoma; stage and tumour size for ccRCC), and for LY75 with sequencing run in both models.

## Data

The targeted bisulfite sequencing data (raw and processed) are available in NCBI GEO under accession [GSE338748](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE338748). [[Status: private with a reviewer token until publication.]]

The clinical data (sample sheet, clinical follow-up) are not included in this repository because they are patient data. [[Access: e.g. available from the corresponding author upon reasonable request, subject to approval by the MUMC ethics committee.]]

The analysis expects these input files in `PROJECT_DIR`:
- `SampleSheet MUMC run.xlsx` (sample sheet of the sequencing runs);
- `SupplementaryTable1.xlsx` (primer and amplicon annotation; Supplementary Table 1 of the paper);
- `Klinische data location III CRC, ccRCC, melanoma- compleet 17022022.xlsx` and `..._edited.xlsx` (clinical data CRC, ccRCC and melanoma);
- `Clinical follow-up data additional melanoma samples - Location-III.xlsx` (clinical data of the melanoma samples of sequencing run 239).

## License

© 2026 Louis Coussement. All rights reserved. This code is provided for review and reproducibility purposes only. No permission is granted to use, copy, modify, or distribute it without the author's written consent. The license will be changed to the MIT License upon acceptance of the manuscript.

## Citation

If you use this code, please cite:

Lommen K\*, Coussement L\*, Vandersmissen J, Herman JG, Van Criekinge W, Smits KM, De Meyer T\*\*, van Engeland M\*\*. Lost in translation: How CpG site selection and assay design determine the clinical value of cancer DNA methylation markers. [[Journal, year, DOI]].

\* shared first authorship; \*\* shared last authorship
