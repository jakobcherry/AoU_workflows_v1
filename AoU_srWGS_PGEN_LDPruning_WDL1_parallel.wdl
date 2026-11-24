version 1.0

# ============================================================
# AoU v9 srWGS PGEN LD pruning + Step 1 genotype preparation
#
# Scatter architecture: chromosomes 1-22 are independent tasks and
# are submitted in parallel (actual concurrency is controlled by
# Terra/Cromwell quotas and backend scheduling).
#
# Per chromosome:
#   1. Download that chromosome's PGEN/PVAR/PSAM trio.
#   2. Apply the user-provided PLINK --keep list.
#   3. Create .prune.in and .prune.out with:
#        --indep-pairwise 200kb 0.5
#      (PLINK2 defaults the omitted step to 1 for kb windows.)
#   4. Extract .prune.in variants from the original chromosome PGEN.
#   5. Emit the pruned chromosome trio, lists, logs, and QC.
#
# After all chromosome tasks succeed:
#   6. Merge all chromosome PGEN trios into one genome-wide PGEN
#      fileset for REGENIE Step 1 using --pmerge-list.
#
# No ancestry assignment, relatedness filtering, sex-discrepancy
# filtering, or keep-list creation is performed in this WDL.
# The user supplies their own keep list.
# ============================================================

workflow AoU_srWGS_PGEN_LDPruning_Step1Prep_Parallel {
  input {
    String pgen_gcs_dir = "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/"
    File keep_list
    String plink2_gcs_path = "gs://cloned-all-circadian-patients-for-fitbit-and-prs-wb-perky-onion/plink2"

    Int window_kb = 200
    Float r2_threshold = 0.5
    String output_prefix = "AoU_srWGS_1pct_LDpruned"

    Int chromosome_cpu = 8
    Int chromosome_memory_gb = 192
    Int chromosome_disk_gb = 500
    Int plink_memory_mb = 160000

    Int merge_cpu = 8
    Int merge_memory_gb = 192
    Int merge_disk_gb = 500
  }


  Array[String] chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"]

  scatter (chr in chromosomes) {
    call PruneExtractOneChromosome {
      input:
        chromosome = chr,
        pgen_gcs_dir = pgen_gcs_dir,
        keep_list = keep_list,
        plink2_gcs_path = plink2_gcs_path,
        window_kb = window_kb,
        r2_threshold = r2_threshold,
        plink_memory_mb = plink_memory_mb,
        task_cpu = chromosome_cpu,
        task_memory_gb = chromosome_memory_gb,
        task_disk_gb = chromosome_disk_gb
    }
  }


  call MergePrunedChromosomesParallel {
    input:
      chromosome_pgens = PruneExtractOneChromosome.pruned_pgen,
      chromosome_pvars = PruneExtractOneChromosome.pruned_pvar,
      chromosome_psams = PruneExtractOneChromosome.pruned_psam,
      prune_ins = PruneExtractOneChromosome.prune_in,
      prune_outs = PruneExtractOneChromosome.prune_out,
      pruning_logs = PruneExtractOneChromosome.pruning_log,
      extraction_logs = PruneExtractOneChromosome.extraction_log,
      chromosome_qc = PruneExtractOneChromosome.chromosome_qc,
      window_kb = window_kb,
      r2_threshold = r2_threshold,
      plink2_gcs_path = plink2_gcs_path,
      output_prefix = output_prefix,
      plink_memory_mb = plink_memory_mb,
      task_cpu = merge_cpu,
      task_memory_gb = merge_memory_gb,
      task_disk_gb = merge_disk_gb
  }


  output {
    File merged_step1_pgen = MergePrunedChromosomesParallel.merged_pgen
    File merged_step1_pvar = MergePrunedChromosomesParallel.merged_pvar
    File merged_step1_psam = MergePrunedChromosomesParallel.merged_psam
    File merged_qc_summary = MergePrunedChromosomesParallel.merged_qc_summary
    File per_chromosome_variant_counts = MergePrunedChromosomesParallel.per_chromosome_variant_counts
    File merge_log = MergePrunedChromosomesParallel.merge_log
    File complete_results_archive = MergePrunedChromosomesParallel.complete_results_archive

    Array[File] chromosome_pruned_pgens = PruneExtractOneChromosome.pruned_pgen
    Array[File] chromosome_pruned_pvars = PruneExtractOneChromosome.pruned_pvar
    Array[File] chromosome_pruned_psams = PruneExtractOneChromosome.pruned_psam
    Array[File] chromosome_prune_in_lists = PruneExtractOneChromosome.prune_in
    Array[File] chromosome_prune_out_lists = PruneExtractOneChromosome.prune_out
    Array[File] chromosome_pruning_logs = PruneExtractOneChromosome.pruning_log
    Array[File] chromosome_extraction_logs = PruneExtractOneChromosome.extraction_log
    Array[File] chromosome_qc_files = PruneExtractOneChromosome.chromosome_qc
  }
}


# ------------------------------------------------------------
# Scatter task: one chromosome per task
# ------------------------------------------------------------
task PruneExtractOneChromosome {
  input {
    String chromosome
    String pgen_gcs_dir
    File keep_list
    String plink2_gcs_path
    Int window_kb
    Float r2_threshold
    Int plink_memory_mb
    Int task_cpu
    Int task_memory_gb
    Int task_disk_gb
  }

  command <<<
    set -euo pipefail

    mkdir -p input logs qc
    prefix="input/acaf_threshold.chr~{chromosome}"
    out_prefix="chr~{chromosome}_ldpruned"

    # Stage PLINK2 and the chromosome-specific PGEN trio.
    gsutil cp "~{plink2_gcs_path}" ./plink2
    chmod +x ./plink2
    gsutil cp "~{pgen_gcs_dir}acaf_threshold.chr~{chromosome}.pgen" "${prefix}.pgen"
    gsutil cp "~{pgen_gcs_dir}acaf_threshold.chr~{chromosome}.pvar" "${prefix}.pvar"
    gsutil cp "~{pgen_gcs_dir}acaf_threshold.chr~{chromosome}.psam" "${prefix}.psam"

    # Step 1: create LD-pruning inclusion/exclusion lists.
    # Keep-list is user-provided; no sample selection logic is performed here.
    ./plink2 \
      --pfile "${prefix}" \
      --keep "~{keep_list}" \
      --indep-pairwise "~{window_kb}kb" "~{r2_threshold}" \
      --memory "~{plink_memory_mb}" \
      --threads "~{task_cpu}" \
      --out "${out_prefix}" \
      > "logs/${out_prefix}_pruning.log" 2>&1

    prune_in="${out_prefix}.prune.in"
    prune_out="${out_prefix}.prune.out"
    [[ -s "${prune_in}" ]] || { echo "ERROR: Missing/empty ${prune_in}" >&2; exit 1; }
    [[ -f "${prune_out}" ]] || { echo "ERROR: Missing ${prune_out}" >&2; exit 1; }

    # Step 2: extract retained SNPs from the original chromosome PGEN.
    ./plink2 \
      --pfile "${prefix}" \
      --keep "~{keep_list}" \
      --extract "${prune_in}" \
      --make-pgen \
      --memory "~{plink_memory_mb}" \
      --threads "~{task_cpu}" \
      --out "${out_prefix}" \
      > "logs/${out_prefix}_extraction.log" 2>&1

    n_prune_in=$(wc -l < "${prune_in}" | tr -d '[:space:]')
    n_prune_out=$(wc -l < "${prune_out}" | tr -d '[:space:]')
    n_extracted=$(awk 'BEGIN{n=0} !/^#/ {n++} END{print n}' "${out_prefix}.pvar")
    n_samples=$(awk 'NR>1 {n++} END{print n+0}' "${out_prefix}.psam")

    printf 'chromosome\tprune_in_count\tprune_out_count\textracted_variant_count\tsample_count\n' > "qc/chr~{chromosome}_variant_counts.tsv"
    printf '~{chromosome}\t%s\t%s\t%s\t%s\n' "${n_prune_in}" "${n_prune_out}" "${n_extracted}" "${n_samples}" >> "qc/chr~{chromosome}_variant_counts.tsv"

    # Clean up original chromosome trio after extraction to free task-local disk.
    rm -f "${prefix}.pgen" "${prefix}.pvar" "${prefix}.psam"
  >>>

  output {
    File pruned_pgen = "chr~{chromosome}_ldpruned.pgen"
    File pruned_pvar = "chr~{chromosome}_ldpruned.pvar"
    File pruned_psam = "chr~{chromosome}_ldpruned.psam"
    File prune_in = "chr~{chromosome}_ldpruned.prune.in"
    File prune_out = "chr~{chromosome}_ldpruned.prune.out"
    File pruning_log = "logs/chr~{chromosome}_ldpruned_pruning.log"
    File extraction_log = "logs/chr~{chromosome}_ldpruned_extraction.log"
    File chromosome_qc = "qc/chr~{chromosome}_variant_counts.tsv"
  }

  runtime {
    cpu: task_cpu
    memory: task_memory_gb + " GB"
    disks: "local-disk " + task_disk_gb + " SSD"
  }
}


# ------------------------------------------------------------
# Merge task: waits for all scatter tasks, then concatenates
# chromosome-separated pruned PGEN filesets into a Step 1 set.
# ------------------------------------------------------------
task MergePrunedChromosomesParallel {
  input {
    Array[File] chromosome_pgens
    Array[File] chromosome_pvars
    Array[File] chromosome_psams
    Array[File] prune_ins
    Array[File] prune_outs
    Array[File] pruning_logs
    Array[File] extraction_logs
    Array[File] chromosome_qc
    Int window_kb
    Float r2_threshold
    String plink2_gcs_path
    String output_prefix
    Int plink_memory_mb
    Int task_cpu
    Int task_memory_gb
    Int task_disk_gb
  }

  command <<<
    set -euo pipefail
    mkdir -p merge_inputs logs qc
    gsutil cp "~{plink2_gcs_path}" ./plink2
    chmod +x ./plink2

    # Files are localized into this task. Build a pmerge-list from PGEN
    # basenames; each prefix must have matching .pvar and .psam alongside it.
    printf '%s\n' chr*_ldpruned.pgen | sort -V | sed 's/\.pgen$//' > pmerge_list.txt
    [[ -s pmerge_list.txt ]] || { echo "ERROR: Empty pmerge list" >&2; exit 1; }

    # Check every prefix has the companion variant/sample files before merging.
    while IFS= read -r prefix; do
      [[ -s "${prefix}.pgen" && -s "${prefix}.pvar" && -s "${prefix}.psam" ]] || {
        echo "ERROR: Incomplete PGEN trio for prefix ${prefix}" >&2; exit 1;
      }
    done < pmerge_list.txt

    ./plink2 \
      --pmerge-list pmerge_list.txt pfile \
      --make-pgen \
      --memory "~{plink_memory_mb}" \
      --threads "~{task_cpu}" \
      --out "~{output_prefix}" \
      > "logs/~{output_prefix}_merge.log" 2>&1

    [[ -s "~{output_prefix}.pgen" && -s "~{output_prefix}.pvar" && -s "~{output_prefix}.psam" ]] || {
      echo "ERROR: Merged PGEN trio not generated" >&2; exit 1;
    }

    n_variants=$(awk 'BEGIN{n=0} !/^#/ {n++} END{print n}' "~{output_prefix}.pvar")
    n_samples=$(awk 'NR>1 {n++} END{print n+0}' "~{output_prefix}.psam")
    {
      echo "Output prefix: ~{output_prefix}"
      echo "Merged variant count: ${n_variants}"
      echo "Merged sample count: ${n_samples}"
      echo "Target variant count: 300000-400000 (QC target only; not enforced)"
      echo "Chromosome inputs merged: $(wc -l < pmerge_list.txt)"
      echo "LD pruning: --indep-pairwise ~{window_kb}kb ~{r2_threshold}"
    } > "qc/~{output_prefix}_merged_qc_summary.txt"

    # Collate chromosome QC into one TSV.
    printf 'chromosome\tprune_in_count\tprune_out_count\textracted_variant_count\tsample_count\n' > "qc/per_chromosome_variant_counts.tsv"
    for f in chr*_variant_counts.tsv; do
      [[ -e "$f" ]] || continue
      tail -n +2 "$f" >> "qc/per_chromosome_variant_counts.tsv"
    done

    # Make one archive with all intermediate files and final merged files.
    tar -czf "~{output_prefix}_complete_results.tar.gz" \
      chr*_ldpruned.pgen chr*_ldpruned.pvar chr*_ldpruned.psam \
      chr*_ldpruned.prune.in chr*_ldpruned.prune.out \
      chr*_ldpruned_pruning.log chr*_ldpruned_extraction.log \
      chr*_variant_counts.tsv pmerge_list.txt logs qc \
      "~{output_prefix}.pgen" "~{output_prefix}.pvar" "~{output_prefix}.psam"
  >>>

  output {
    File merged_pgen = "~{output_prefix}.pgen"
    File merged_pvar = "~{output_prefix}.pvar"
    File merged_psam = "~{output_prefix}.psam"
    File merged_qc_summary = "qc/~{output_prefix}_merged_qc_summary.txt"
    File per_chromosome_variant_counts = "qc/per_chromosome_variant_counts.tsv"
    File merge_log = "logs/~{output_prefix}_merge.log"
    File complete_results_archive = "~{output_prefix}_complete_results.tar.gz"
  }

  runtime {
    cpu: task_cpu
    memory: task_memory_gb + " GB"
    disks: "local-disk " + task_disk_gb + " SSD"
  }
}
