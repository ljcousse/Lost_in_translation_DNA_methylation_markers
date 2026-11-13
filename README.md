# Lost in translation: DNA methylation markers

[![DOI](https://zenodo.org/badge/DOI/[[ZENODO_DOI]].svg)](https://doi.org/[[ZENODO_DOI]])

Code for the manuscript *"Lost in translation: How CpG site selection and assay design determine the clinical value of cancer DNA methylation markers"* (Lommen, Coussement et al., [[journal, year, DOI]]).

We profiled three established cancer DNA methylation markers at single-CpG resolution with targeted (amplicon) bisulfite sequencing: NDRG4 (colorectal cancer, diagnostic), GREM1 (clear-cell renal cell carcinoma, prognostic) and LY75 (melanoma, prognostic). For every CpG, amplicon, read-based measure and MSP assay, we quantified the clinical value with the Akaike Information Criterion (ΔAIC compared with a baseline model). We then compared CpGs within the same amplicon with a permutation test.

## Repository contents

| File | Purpose |
|---|---|
| `TargetedBisulphiteSeq.bash` | Pre-processing: quality control (FastQC, MultiQC), adapter and quality trimming (Trim Galore), bisulfite alignment and methylation calling (Bismark/Bowtie2, GRCh37), read-level methylation patterns (wgbs_tools `bam2pat`) |
| `AmpliconAwareMapping.R` | Assigns every read pair to the amplicon it originates from, so that overlapping amplicons are analysed separately |
| `TargetedBisulphite_analysis.R` | Main analysis: M-values per CpG, ΔAIC per CpG, amplicon and MSP assay, read-based measures (relative co-methylation height and density), CpGs covered by two amplicons, primer-CpG analyses, intra-amplicon comparisons and the tables/figures of the manuscript |
| `Intra_amplicon_diff_functions.R` | Functions for the intra-amplicon permutation test (sourced by the main analysis) |
| `Plots_AIC_function.R` | Plot function for the ΔAIC figures (sourced by the main analysis) |
| `scripts/cpg_density.R`, `scripts/plot_cpg_density.R` | CpG density per amplicon region (revision, reviewer 3) |
| `conda_environment.yml` | Conda environment with the tools and package versions used |
| `config.example.env` | Template for the local paths (copy to `config.env`) |

## Installation

```bash
git clone https://github.com/ljcousse/Lost_in_translation_DNA_methylation_markers.git
cd Lost_in_translation_DNA_methylation_markers

conda env create -f conda_environment.yml
conda activate lost_in_translation
```

`wgbs_tools` is not available via conda. Install it from https://github.com/nloyfer/wgbs_tools and initialise the GRCh37 reference (`wgbstools init_genome ... hg37_87`; see `TargetedBisulphiteSeq.bash`).

The analysis was run on a Linux server (AlmaLinux 9) with R 4.2.3 and Bioconductor 3.16. The exact versions are in `conda_environment.yml`.

## Configuration

Local paths are not hard-coded. Copy the template and fill in your paths:

```bash
cp config.example.env config.env
```

`TargetedBisulphite_analysis.R` reads `PROJECT_DIR` from `config.env`. Use another file with `PROJECT_CONFIG=/path/to/file.env`. Keep `config.env` local (it is listed in `.gitignore`).

## Running the analysis

1. **Pre-processing:** `bash TargetedBisulphiteSeq.bash` (FASTQ → trimmed reads → Bismark alignments and methylation calls → wgbs_tools `.pat` files).
2. **Amplicon assignment:** `Rscript AmpliconAwareMapping.R`.
3. **Analysis:** `Rscript TargetedBisulphite_analysis.R`. Output goes to `PROJECT_DIR`, and a time-stamped log file records the console output and session information. The intra-amplicon permutation test, bootstrap and swap test use `parallel::mclapply` with 40 cores; change `mc.cores` for smaller machines.

### Main outputs

| Output | Content |
|---|---|
| `NDRG4_AIC_table.txt`, `GREM1_AIC_table.txt`, `LY75_AIC_table.txt` | ΔAIC per CpG (Supplementary Table 4) |
| `avAIC_table.txt`, `avAIC_table_ext.txt` | Mean and minimum ΔAIC per amplicon and read-based ΔAIC (Table 3) |
| `pval_total_list.txt` | Intra-amplicon comparisons of all CpG pairs (Figure 4) |
| `*_AIC_table_sensitivity.txt`, `*_AIC_table_nonnested.txt` | Sensitivity analyses: unadjusted and with clinical covariates (Breslow thickness; stage and tumour size) |
| `intra_amplicon_bootstrap_stability.txt`, `intra_amplicon_swap_test.txt` | Stability of the intra-amplicon differences |
| `AIC_plots/` | ΔAIC figures, including the Table 3 figure |

### Statistical design (short)

- **NDRG4 (diagnostic, paired tumour/normal):** logistic regression. Tissue ~ CpG methylation + patient ID, compared with tissue ~ patient ID.
- **GREM1 and LY75 (prognostic):** Cox regression. Methylation alone is compared with a baseline model with age and sex (competing models, same samples). A negative ΔAIC means that the methylation model explains the outcome better than the baseline model.

## Data

The targeted bisulfite sequencing data (raw and processed) are available in NCBI GEO under accession [GSE338748](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE338748). [[Status: private with a reviewer token until publication.]]

The clinical data (sample sheet, clinical follow-up) are not included in this repository because they are patient data. [[Access: e.g. available from the corresponding author upon reasonable request, subject to approval by the MUMC ethics committee.]]

The analysis expects these input files in `PROJECT_DIR`:
- `SampleSheet MUMC run.xlsx`;
- `SupplementaryTable1.xlsx` (primer and amplicon annotation; Supplementary Table 1 of the paper);
- the clinical data files.

## License

© 2026 Louis Coussement. All rights reserved. This code is provided for review and reproducibility purposes only. No permission is granted to use, copy, modify, or distribute it without the author's written consent. The license will be changed to the MIT License upon acceptance of the manuscript.

## Citation

If you use this code, please cite:

Lommen K\*, Coussement L\*, Vandersmissen J, Herman JG, Van Criekinge W, Smits KM, De Meyer T\*\*, van Engeland M\*\*. Lost in translation: How CpG site selection and assay design determine the clinical value of cancer DNA methylation markers. [[Journal, year, DOI]].

\* shared first authorship; \*\* shared last authorship

Code archive: [[Zenodo DOI]].
