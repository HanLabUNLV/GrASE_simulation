# GrASE Simulation RNA-Seq Analysis Pipeline

This repository contains scripts to download simulated RNA-seq data and run the STAR alignment pipeline.

## Prerequisites

Ensure the following tools are installed and available in your PATH:
*   `wget`
*   `tar`
*   `STAR` (v2.7.10b)
*   `parallel`
*   `awk`

## 1. Data Download and Preparation

### Download Data
Run the `wget.sh` script to download the simulated datasets from Zenodo.
```bash
bash wget.sh
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
