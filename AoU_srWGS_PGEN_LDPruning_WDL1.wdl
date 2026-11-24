version 1.0

# ============================================================
# AoU srWGS PGEN LD pruning + genome-wide Step 1 dataset prep
#
# Source: AoU v9 srWGS, ACaf-threshold PGEN trio per chromosome.
# User supplies a PLINK-compatible --keep file (FID IID; no header).
#
# Per chromosome, sequentially:
#   1. Download chromosome PGEN/PVAR/PSAM
#   2. LD-prune using --indep-pairwise 200kb 0.5 (step defaults to 1)
#   3. Extract variants listed in *.prune.in into a chromosome PGEN trio
#   4. Remove downloaded original PGEN trio before moving on
#
# Then merge extracted chromosome trios with --pmerge-list.
# This workflow intentionally does not create ancestry or sample keep-lists.
# ============================================================

workflow AoU_srWGS_PGEN_LDPruning_Step1Prep {
  input {
    String pgen_gcs_dir = "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/"
    File keep_list
    String plink2_gcs_path = "gs://cloned-all-circadian-patients-for-fitbit-and-prs-wb-perky-onion/plink2"

    Int window_kb = 200
    Float r2_threshold = 0.5
    Boolean include_chrX = false
    String output_prefix = "AoU_srWGS_1pct_LDpruned"

    Int task_cpu = 8
    Int task_memory_gb = 192
    Int task_disk_gb = 500
    Int plink_memory_mb = 160000
  }


  call PruneAndExtractSequential {
    input:
      pgen_gcs_dir = pgen_gcs_dir,
      keep_list = keep_list,
      plink2_gcs_path = plink2_gcs_path,
      window_kb = window_kb,
      r2_threshold = r2_threshold,
      include_chrX = include_chrX,
      plink_memory_mb = plink_memory_mb,
      task_cpu = task_cpu,
      task_memory_gb = task_memory_gb,
      task_disk_gb = task_disk_gb
  }

  call MergePrunedChromosomes {
    input:
      prune_bundle = PruneAndExtractSequential.prune_bundle,
      plink2_gcs_path = plink2_gcs_path,
      output_prefix = output_prefix,
      plink_memory_mb = plink_memory_mb,
      task_cpu = task_cpu,
      task_memory_gb = task_memory_gb,
      task_disk_gb = task_disk_gb
  }

  output {
    File merged_step1_pgen = MergePrunedChromosomes.merged_pgen
    File merged_step1_pvar = MergePrunedChromosomes.merged_pvar
    File merged_step1_psam = MergePrunedChromosomes.merged_psam
    File merged_qc_summary = MergePrunedChromosomes.merged_qc_summary
    File merge_log = MergePrunedChromosomes.merge_log
    File all_prune_lists_and_chromosome_pgens = PruneAndExtractSequential.prune_bundle
    File per_chromosome_variant_counts = PruneAndExtractSequential.per_chromosome_variant_counts
    File pruning_settings_summary = PruneAndExtractSequential.pruning_settings_summary
    File complete_results_archive = MergePrunedChromosomes.complete_results_archive
  }
}


# ------------------------------------------------------------
# Task 1: Sequential pruning and extraction (one chromosome at a time)
# ------------------------------------------------------------
task PruneAndExtractSequential {
  input {
    String pgen_gcs_dir
    File keep_list
    String plink2_gcs_path
    Int window_kb
    Float r2_threshold
    Boolean include_chrX
    Int plink_memory_mb
    Int task_cpu
    Int task_memory_gb
    Int task_disk_gb
  }


  command <<<
    set -euo pipefail

    mkdir -p original chromosome_pruned prune_lists logs qc

    gsutil cp "~{plink2_gcs_path}" ./plink2
    chmod +x ./plink2

    KEEP_FILE="~{keep_list}"

    chromosomes=( {1..22} )
    if [[ "~{include_chrX}" == "true" ]]; then
      chromosomes+=(X)
    fi

    printf 'chromosome\tprune_in_count\tprune_out_count\textracted_variant_count\n' > qc/per_chromosome_variant_counts.tsv

    for chr in "${chromosomes[@]}"; do
      echo "---- Processing chromosome ${chr} ----"
      prefix="original/acaf_threshold.chr${chr}"

      gsutil cp "~{pgen_gcs_dir}acaf_threshold.chr${chr}.pgen" "${prefix}.pgen"
      gsutil cp "~{pgen_gcs_dir}acaf_threshold.chr${chr}.pvar" "${prefix}.pvar"
      gsutil cp "~{pgen_gcs_dir}acaf_threshold.chr${chr}.psam" "${prefix}.psam"

      # Requested PLINK2 two-argument form. With a kb window, PLINK2
      # defaults the variant-count step to 1 (and requires it to be 1).
      ./plink2 \
        --pfile "${prefix}" \
        --keep "${KEEP_FILE}" \
        --indep-pairwise "~{window_kb}kb" "~{r2_threshold}" \
        --memory "~{plink_memory_mb}" \
        --threads "~{task_cpu}" \
        --out "prune_lists/chr${chr}" \
        > "logs/chr${chr}_pruning.log" 2>&1

      prune_in="prune_lists/chr${chr}.prune.in"
      prune_out="prune_lists/chr${chr}.prune.out"
      if [[ ! -s "${prune_in}" ]]; then
        echo "ERROR: Missing/empty ${prune_in}" >&2
        exit 1
      fi
      if [[ ! -f "${prune_out}" ]]; then
        echo "ERROR: Missing ${prune_out}" >&2
        exit 1
      fi

      # Apply the retained-variant list to the original chromosome PGEN.
      ./plink2 \
        --pfile "${prefix}" \
        --keep "${KEEP_FILE}" \
        --extract "${prune_in}" \
        --make-pgen \
        --memory "~{plink_memory_mb}" \
        --threads "~{task_cpu}" \
        --out "chromosome_pruned/chr${chr}" \
        > "logs/chr${chr}_extraction.log" 2>&1

      in_count=$(wc -l < "${prune_in}" | tr -d '[:space:]')
      out_count=$(wc -l < "${prune_out}" | tr -d '[:space:]')
      extracted_count=$(awk 'BEGIN{n=0} !/^#/ {n++} END{print n}' "chromosome_pruned/chr${chr}.pvar")
      printf '%s\t%s\t%s\t%s\n' "${chr}" "${in_count}" "${out_count}" "${extracted_count}" >> qc/per_chromosome_variant_counts.tsv

      # Remove the original chromosome files to reduce local disk occupancy.
      rm -f "${prefix}.pgen" "${prefix}.pvar" "${prefix}.psam"
    done


    total_in=$(awk 'NR>1 {s+=$2} END{print s+0}' qc/per_chromosome_variant_counts.tsv)
    total_out=$(awk 'NR>1 {s+=$3} END{print s+0}' qc/per_chromosome_variant_counts.tsv)
    total_extracted=$(awk 'NR>1 {s+=$4} END{print s+0}' qc/per_chromosome_variant_counts.tsv)

    {
      echo "LD pruning command: --indep-pairwise ~{window_kb}kb ~{r2_threshold}"
      echo "Step size: PLINK2 default 1 variant (required with kb window)"
      echo "Chromosome X included: ~{include_chrX}"
      echo "Total prune.in variants: ${total_in}"
      echo "Total prune.out variants: ${total_out}"
      echo "Total extracted variants (sum across chromosomes): ${total_extracted}"
      echo "Desired Step 1 marker count: 300000-400000 (QC target only; not enforced)"
    } > qc/pruning_settings_summary.txt

    # Bundle all prune lists, chromosome-level extracted PGEN trios, logs, and QC.
    tar -czf prune_extract_bundle.tar.gz chromosome_pruned prune_lists logs qc
  >>>

  output {
    File prune_bundle = "prune_extract_bundle.tar.gz"
    File per_chromosome_variant_counts = "qc/per_chromosome_variant_counts.tsv"
    File pruning_settings_summary = "qc/pruning_settings_summary.txt"
  }

  runtime {
    cpu: task_cpu
    memory: task_memory_gb + " GB"
    disks: "local-disk " + task_disk_gb + " SSD"
  }
}


# ------------------------------------------------------------
# Task 2: Merge per-chromosome extracted PGENs into one Step 1 set
# ------------------------------------------------------------
task MergePrunedChromosomes {
  input {
    File prune_bundle
    String plink2_gcs_path
    String output_prefix
    Int plink_memory_mb
    Int task_cpu
    Int task_memory_gb
    Int task_disk_gb
  }

  command <<<
    set -euo pipefail

    mkdir -p merge_work
    tar -xzf "~{prune_bundle}" -C merge_work
    cd merge_work

    gsutil cp "~{plink2_gcs_path}" ./plink2
    chmod +x ./plink2

    # Prefix list: PLINK2 --pmerge-list defaults to PGEN/PVAR/PSAM prefixes.
    find chromosome_pruned -maxdepth 1 -type f -name 'chr*.pgen' -print \
      | sort -V \
      | sed 's/\.pgen$//' \
      > pmerge_list.txt

    if [[ ! -s pmerge_list.txt ]]; then
      echo "ERROR: No extracted chromosome PGEN prefixes found for merge" >&2
      exit 1
    fi

    ./plink2 \
      --pmerge-list pmerge_list.txt pfile \
      --memory "~{plink_memory_mb}" \
      --threads "~{task_cpu}" \
      --out "~{output_prefix}" \
      > "~{output_prefix}_merge.log" 2>&1

    if [[ ! -s "~{output_prefix}.pgen" || ! -s "~{output_prefix}.pvar" || ! -s "~{output_prefix}.psam" ]]; then
      echo "ERROR: Expected merged PGEN trio not produced" >&2
      exit 1
    fi

    n_variants=$(awk 'BEGIN{n=0} !/^#/ {n++} END{print n}' "~{output_prefix}.pvar")
    n_samples=$(awk 'NR>1 {n++} END{print n+0}' "~{output_prefix}.psam")
    {
      echo "Merged prefix: ~{output_prefix}"
      echo "Merged variant count: ${n_variants}"
      echo "Merged sample count: ${n_samples}"
      echo "Target variant count: 300000-400000 (QC target only)"
    } > "~{output_prefix}_merged_qc_summary.txt"

    # Final archive contains all chromosome intermediates and the merged Step 1 trio.
    tar -czf "~{output_prefix}_complete_results.tar.gz" \
      chromosome_pruned prune_lists logs qc pmerge_list \
      "~{output_prefix}.pgen" "~{output_prefix}.pvar" "~{output_prefix}.psam" \
      "~{output_prefix}_merge.log" "~{output_prefix}_merged_qc_summary.txt"
  >>>

  output {
    File merged_pgen = "merge_work/~{output_prefix}.pgen"
    File merged_pvar = "merge_work/~{output_prefix}.pvar"
    File merged_psam = "merge_work/~{output_prefix}.psam"
    File merged_qc_summary = "merge_work/~{output_prefix}_merged_qc_summary.txt"
    File merge_log = "merge_work/~{output_prefix}_merge.log"
    File complete_results_archive = "merge_work/~{output_prefix}_complete_results.tar.gz"
  }

  runtime {
    cpu: task_cpu
    memory: task_memory_gb + " GB"
    disks: "local-disk " + task_disk_gb + " SSD"
  }
}
