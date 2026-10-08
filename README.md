# GrASE Simulation RNA-Seq Analysis Pipeline

This repository contains the drivers that reproduce the GrASE benchmark: download
the simulated RNA-seq data of Love et al., align it, run GrASE, and produce every
figure and table reported.

**Deposited inputs and results.** The splicing graphs, bipartitions, per-gene GTFs
and DEXSeq GFFs for human GENCODE v28 and v34, the simulation ground truth, and the
significant-bipartition tables are at
[doi:10.5281/zenodo.23141907](https://doi.org/10.5281/zenodo.23141907). Stages that
depend only on the annotation can be skipped by downloading them. The GrASE R
package itself is archived at
[doi:10.5281/zenodo.23169081](https://doi.org/10.5281/zenodo.23169081).

Note the two Zenodo sources are different things: `data/wget.sh` fetches the
simulated FASTQ archives of Love et al. (records 1291375, 1291404, 1291443), while
the record above holds the annotation-derived inputs and our results.

## Prerequisites

Ensure the following tools are installed and available in your PATH:
*   `wget`
*   `tar`
*   `STAR` (v2.7.10b)
*   `parallel`
*   `awk`

## 1. Data Download and Preparation

### Download Data
Run the `wget.sh` script to download the simulated datasets of Love et al. from
Zenodo. Edit it first: most of the twelve archives are commented out.
```bash
bash data/wget.sh
```

### Extract Archives
Extract the downloaded tar files.
```bash
bash tar.sh
```

### Organize Files
Create symbolic links to organize the files into `group1` (containing `sample_01`) and `group2` (containing `sample_02`).
```bash
bash ln.sh
```

## 2. STAR Pipeline

The pipeline is located in the `STAR/` directory.

### Step 1: Generate Genome Index (Pass 1)
Generate the initial genome index.
```bash
bash STAR/00_genomegenerate.sh
```
*Note: Ensure `GRCh38.p13.genome.fa` and `gencode.v34.annotation.gtf` are present in the `STAR` directory or update the paths in the script.*

### Step 2: Run STAR Pass 1
Run the first pass of STAR alignment for both groups.
**Important:** You may need to modify `STAR/02_pass1.sh` to point to the correct sample filenames for each group (`sample_01` for `group1`, `sample_02` for `group2`).

```bash
# Example for a single sample
bash STAR/02_pass1.sh group1 out_1
```
To run in parallel using `group.list`:
```bash
parallel -j 12 bash STAR/02_pass1.sh group1 {} :::: group.list
parallel -j 12 bash STAR/02_pass1.sh group2 {} :::: group.list
```

### Step 3: Process Splice Junctions
Merge and filter the splice junctions detected in Pass 1.

1.  **Merge raw SJ tabs** for each group:
    ```bash
    bash STAR/04_merge_sj.sh group1
    bash STAR/04_merge_sj.sh group2
    ```
2.  **Filter for high confidence junctions**:
    ```bash
    bash STAR/05_filter_sj.sh group1
    bash STAR/05_filter_sj.sh group2
    ```

3.  **Merge all cell types**:
    ```bash
    bash STAR/06_merge_all_celltype.sh
    ```

### Step 4: Generate Genome Index (Pass 2)
Generate the second genome index using the merged splice junctions from all samples.
```bash
bash STAR/07_genomegenerate.sh
```

### Step 5: Run STAR Pass 2
Run the second pass of alignment.
**Important:** Similar to Step 2, ensure `STAR/08_pass2.sh` points to the correct sample filenames.

```bash
parallel -j 12 bash STAR/08_pass2.sh group1 {} :::: group.list
parallel -j 12 bash STAR/08_pass2.sh group2 {} :::: group.list
```

## Notes
*   **File Paths**: The scripts currently contain hardcoded paths (e.g., `$HOME/Love_simulation/...`). You may need to adjust these variables in the scripts to match your directory structure.
*   **Sample Naming**: `02_pass1.sh` and `08_pass2.sh` may have hardcoded strings for `sample_01` or `sample_02`. Verify these match your data before running.

## 3. GrASE analysis pipeline

Drivers in `scripts/` are numbered in execution order. Each driver can be rerun on
its own provided its inputs exist.

| stage | driver | what it does |
|-------|--------|--------------|
| 10 | `10_counts_nc2_multinomial.sh` | exonic-part counts for the n_choose_2 and multinomial splits |
| 20 | `20_tests_bipartition_minreads.sh` | the four bipartition exontest arms, with the pre-testing read floor |
| 21 | `21_tests_nc2_multinomial.sh` | exontest for n_choose_2 and multinomial |
| 22 | `22_tests_multinomial.sh` | multinomial arms only |
| 30 | `30_merge_downstream.sh` | everything downstream of the merged (exon + split-read) exontest; also builds the gene-level tables |
| 41 | `41_gt_junctions_stranded.sh` | junction-level structural ground truth, stranded rMATS runs |
| 50 | `50_eval_nc2.sh` | evaluate the n_choose_2 results |
| 51 | `51_eval_multinomial.sh` | evaluate the multinomial results |
| 60 | `60_modelcomp_internal.sh` | model comparison, internal arm |
| 61 | `61_modelcomp_tsstts.sh` | model comparison, TSS/TTS arm |
| 62 | `62_modelcomp_gtrule_stratified.sh` | precision and recall per model (EBapprox, EBmap, MLE, wilcoxon), split by simulation category |
| 70 | `70_cross_majiq_modulize.sh` | voila modulize: MAJIQ event-type classification |
| 73 | `73_cross_confusion.sh` | precision and recall for the three comparison structures (bipartition, n_choose_2, multinomial), with null-gene confusion counts; gtI by design, see below |
| 80 | `80_plots_model_comparison.sh` | within-bipartition dispersion and effect-size scatter, EBapprox vs EBmap vs MLE |
| 81 | `81_plots_delta_sweep.sh` | null-gene false positives and DTE/DTU precision as the `lfc_diff_net` threshold is swept |
| 84 | `84_plots_roc_fp.sh` | partial ROC in transcript and gene space, and the FP-location panel |
| 85 | `85_tables_manuscript.sh` | native-unit, structural-reach, transcript-level, gene-level and TSS/TTS attribution tables |
| 86 | `86_plots_dte_dtu.sh` | precision-recall and partial ROC with DTE and DTU as rows, one file per universe |


