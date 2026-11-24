version 1.0

# ============================================================
# AoU v9 srWGS PGEN LD pruning + Step 1 genotype preparation
#
# 22 chromosome-specific tasks run in parallel. Cromwell stages
# every File input; this workflow does not call gsutil in a task.
#
# Per chromosome:
#   1. Use the staged chromosome PGEN/PVAR/PSAM trio.
#   2. Apply the supplied participant keep-list.
#   3. Generate .prune.in/.prune.out using:
#        --indep-pairwise 200kb 0.5
#   4. Extract retained variants into a chromosome PGEN trio.
#
# After all 22 chromosome tasks succeed:
#   5. Merge the chromosome PGEN trios into one genome-wide PGEN.
#
# No ancestry assignment, relatedness filtering, sex-discrepancy
# filtering, or keep-list creation is performed here.
# ============================================================

workflow AoU_srWGS_PGEN_LDPruning_Step1Prep_Parallel {
  input {
    File keep_list

    File plink2_file = "gs://cloned-all-circadian-patients-for-fitbit-and-prs-wb-perky-onion/plink2"

    Array[File] pgen_files = [
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr1.pgen",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr2.pgen",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr3.pgen",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr4.pgen",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr5.pgen",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr6.pgen",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr7.pgen",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr8.pgen",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr9.pgen",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr10.pgen",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr11.pgen",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr12.pgen",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr13.pgen",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr14.pgen",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr15.pgen",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr16.pgen",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr17.pgen",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr18.pgen",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr19.pgen",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr20.pgen",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr21.pgen",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr22.pgen"
    ]

    Array[File] pvar_files = [
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr1.pvar",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr2.pvar",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr3.pvar",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr4.pvar",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr5.pvar",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr6.pvar",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr7.pvar",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr8.pvar",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr9.pvar",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr10.pvar",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr11.pvar",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr12.pvar",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr13.pvar",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr14.pvar",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr15.pvar",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr16.pvar",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr17.pvar",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr18.pvar",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr19.pvar",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr20.pvar",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr21.pvar",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr22.pvar"
    ]

    Array[File] psam_files = [
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr1.psam",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr2.psam",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr3.psam",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr4.psam",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr5.psam",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr6.psam",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr7.psam",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr8.psam",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr9.psam",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr10.psam",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr11.psam",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr12.psam",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr13.psam",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr14.psam",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr15.psam",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr16.psam",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr17.psam",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr18.psam",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr19.psam",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr20.psam",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr21.psam",
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr22.psam"
    ]

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

  Array[String] chromosomes = [
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
    "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"
  ]

  scatter (i in range(22)) {
    call PruneExtractOneChromosome {
      input:
        chromosome = chromosomes[i],
        pgen_file = pgen_files[i],
        pvar_file = pvar_files[i],
        psam_file = psam_files[i],
        keep_list = keep_list,
        plink2_file = plink2_file,
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
      plink2_file = plink2_file,
      output_prefix = output_prefix,
      window_kb = window_kb,
      r2_threshold = r2_threshold,
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
# One task per chromosome
# ------------------------------------------------------------
task PruneExtractOneChromosome {
  input {
    String chromosome
    File pgen_file
    File pvar_file
    File psam_file
    File keep_list
    File plink2_file
    Int window_kb
    Float r2_threshold
    Int plink_memory_mb
    Int task_cpu
    Int task_memory_gb
    Int task_disk_gb
  }

  command <<<
    set -euo pipefail

    exec > >(tee "chr~{chromosome}_task.log") 2>&1

    echo "============================================================"
    echo "AoU PGEN LD PRUNING"
    echo "Chromosome: ~{chromosome}"
    echo "Start: $(date -u)"
    echo "============================================================"

    mkdir -p input logs qc

    # Cromwell has localized these File inputs. Copy the executable
    # to a task-local path and link the large genotype files.
    cp "~{plink2_file}" ./plink2
    chmod 750 ./plink2

    prefix="input/acaf_threshold.chr~{chromosome}"
    out_prefix="chr~{chromosome}_ldpruned"

    ln -s "~{pgen_file}" "${prefix}.pgen"
    ln -s "~{pvar_file}" "${prefix}.pvar"
    ln -s "~{psam_file}" "${prefix}.psam"

    echo "PLINK2 version:"
    ./plink2 --version

    echo "Staged input files:"
    ls -lh ./plink2 "~{keep_list}" \
      "${prefix}.pgen" "${prefix}.pvar" "${prefix}.psam"

    [[ -s "~{keep_list}" ]] || { echo "ERROR: keep-list is missing or empty" >&2; exit 1; }
    [[ -s "${prefix}.pgen" && -s "${prefix}.pvar" && -s "${prefix}.psam" ]] || {
      echo "ERROR: One or more staged chromosome files are missing/empty" >&2
      exit 1
    }

    # Step 1: LD-pruning variant lists.
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

    # Step 2: extract retained variants and the selected participants.
    ./plink2 \
      --pfile "${prefix}" \
      --keep "~{keep_list}" \
      --extract "${prune_in}" \
      --make-pgen \
      --memory "~{plink_memory_mb}" \
      --threads "~{task_cpu}" \
      --out "${out_prefix}" \
      > "logs/${out_prefix}_extraction.log" 2>&1

    [[ -s "${out_prefix}.pgen" && -s "${out_prefix}.pvar" && -s "${out_prefix}.psam" ]] || {
      echo "ERROR: PLINK2 did not produce a complete output PGEN trio" >&2
      tail -100 "logs/${out_prefix}_extraction.log" >&2 || true
      exit 1
    }

    n_prune_in=$(wc -l < "${prune_in}" | tr -d '[:space:]')
    n_prune_out=$(wc -l < "${prune_out}" | tr -d '[:space:]')
    n_extracted=$(awk 'BEGIN{n=0} !/^#/ {n++} END{print n}' "${out_prefix}.pvar")
    n_samples=$(awk 'NR>1 {n++} END{print n+0}' "${out_prefix}.psam")

    printf 'chromosome\\tprune_in_count\\tprune_out_count\\textracted_variant_count\\tsample_count\\n' > "qc/chr~{chromosome}_variant_counts.tsv"
    printf '~{chromosome}\\t%s\\t%s\\t%s\\t%s\\n' \
      "${n_prune_in}" "${n_prune_out}" "${n_extracted}" "${n_samples}" \
      >> "qc/chr~{chromosome}_variant_counts.tsv"

    echo "Chromosome ~{chromosome} completed."
    echo "Prune-in variants: ${n_prune_in}"
    echo "Prune-out variants: ${n_prune_out}"
    echo "Extracted variants: ${n_extracted}"
    echo "Samples: ${n_samples}"
    echo "End: $(date -u)"

    # Remove only the symlinks; Cromwell-managed staged source files
    # are not modified.
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
    File task_log = "chr~{chromosome}_task.log"
  }

  runtime {
    cpu: task_cpu
    memory: task_memory_gb + " GB"
    disks: "local-disk " + task_disk_gb + " SSD"
  }
}

# ------------------------------------------------------------
# Merge all 22 chromosome-specific pruned PGEN trios
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
    File plink2_file
    String output_prefix
    Int window_kb
    Float r2_threshold
    Int plink_memory_mb
    Int task_cpu
    Int task_memory_gb
    Int task_disk_gb
  }

  command <<<
    set -euo pipefail

    exec > >(tee "merge_task.log") 2>&1

    mkdir -p logs qc merge_inputs

    cp "~{plink2_file}" ./plink2
    chmod 750 ./plink2

    echo "PLINK2 version:"
    ./plink2 --version

    # Link every staged input into this task's working directory.
    # Basenames are chromosome-specific and therefore unique.
    for f in ~{sep=' ' chromosome_pgens}; do
      ln -s "$f" "./$(basename "$f")"
    done
    for f in ~{sep=' ' chromosome_pvars}; do
      ln -s "$f" "./$(basename "$f")"
    done
    for f in ~{sep=' ' chromosome_psams}; do
      ln -s "$f" "./$(basename "$f")"
    done
    for f in ~{sep=' ' prune_ins}; do
      ln -s "$f" "./$(basename "$f")"
    done
    for f in ~{sep=' ' prune_outs}; do
      ln -s "$f" "./$(basename "$f")"
    done
    for f in ~{sep=' ' pruning_logs}; do
      ln -s "$f" "./merge_inputs/$(basename "$f")"
    done
    for f in ~{sep=' ' extraction_logs}; do
      ln -s "$f" "./merge_inputs/$(basename "$f")"
    done
    for f in ~{sep=' ' chromosome_qc}; do
      ln -s "$f" "./merge_inputs/$(basename "$f")"
    done

    # Explicit chromosome-order list avoids relying on glob ordering.
    printf '%s\\n' \
      chr1_ldpruned \
      chr2_ldpruned \
      chr3_ldpruned \
      chr4_ldpruned \
      chr5_ldpruned \
      chr6_ldpruned \
      chr7_ldpruned \
      chr8_ldpruned \
      chr9_ldpruned \
      chr10_ldpruned \
      chr11_ldpruned \
      chr12_ldpruned \
      chr13_ldpruned \
      chr14_ldpruned \
      chr15_ldpruned \
      chr16_ldpruned \
      chr17_ldpruned \
      chr18_ldpruned \
      chr19_ldpruned \
      chr20_ldpruned \
      chr21_ldpruned \
      chr22_ldpruned > pmerge_list.txt

    [[ $(wc -l < pmerge_list.txt) -eq 22 ]] || {
      echo "ERROR: Expected 22 chromosome prefixes in pmerge list" >&2
      exit 1
    }

    while IFS= read -r prefix; do
      [[ -s "${prefix}.pgen" && -s "${prefix}.pvar" && -s "${prefix}.psam" ]] || {
        echo "ERROR: Incomplete PGEN trio for ${prefix}" >&2
        exit 1
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
      echo "ERROR: Merged PGEN trio was not generated" >&2
      tail -100 "logs/~{output_prefix}_merge.log" >&2 || true
      exit 1
    }

    n_variants=$(awk 'BEGIN{n=0} !/^#/ {n++} END{print n}' "~{output_prefix}.pvar")
    n_samples=$(awk 'NR>1 {n++} END{print n+0}' "~{output_prefix}.psam")

    {
      echo "Output prefix: ~{output_prefix}"
      echo "Merged variant count: ${n_variants}"
      echo "Merged sample count: ${n_samples}"
      echo "Chromosome inputs merged: 22"
      echo "LD pruning: --indep-pairwise ~{window_kb}kb ~{r2_threshold}"
    } > "qc/~{output_prefix}_merged_qc_summary.txt"

    printf 'chromosome\\tprune_in_count\\tprune_out_count\\textracted_variant_count\\tsample_count\\n' > "qc/per_chromosome_variant_counts.tsv"
    for f in merge_inputs/chr*_variant_counts.tsv; do
      [[ -e "$f" ]] || continue
      tail -n +2 "$f" >> "qc/per_chromosome_variant_counts.tsv"
    done

    tar -czf "~{output_prefix}_complete_results.tar.gz" \
      chr*_ldpruned.pgen chr*_ldpruned.pvar chr*_ldpruned.psam \
      chr*_ldpruned.prune.in chr*_ldpruned.prune.out \
      merge_inputs logs qc pmerge_list.txt merge_task.log \
      "~{output_prefix}.pgen" "~{output_prefix}.pvar" "~{output_prefix}.psam"

    echo "Merge completed."
    echo "Merged variants: ${n_variants}"
    echo "Merged samples: ${n_samples}"
  >>>

  output {
    File merged_pgen = "~{output_prefix}.pgen"
    File merged_pvar = "~{output_prefix}.pvar"
    File merged_psam = "~{output_prefix}.psam"
    File merged_qc_summary = "qc/~{output_prefix}_merged_qc_summary.txt"
    File per_chromosome_variant_counts = "qc/per_chromosome_variant_counts.tsv"
    File merge_log = "logs/~{output_prefix}_merge.log"
    File complete_results_archive = "~{output_prefix}_complete_results.tar.gz"
    File merge_task_log = "merge_task.log"
  }

  runtime {
    cpu: task_cpu
    memory: task_memory_gb + " GB"
    disks: "local-disk " + task_disk_gb + " SSD"
  }
}
