version 1.0

# ============================================================
# AoU v9 srWGS PGEN LD pruning + Step 1 genotype preparation
# CHROMOSOME 1 PILOT
#
# Changes from version 1.0:
#   - Chromosome 1 only; no scatter and no merge task.
#   - More stringent LD pruning: 500kb, r2 = 0.2.
#   - Reduced initial resource allocation for benchmarking.
#   - Retains the same EUR participant keep-list.
#
# Outputs:
#   - chr1_ldpruned.pgen
#   - chr1_ldpruned.pvar
#   - chr1_ldpruned.psam
#   - chr1_ldpruned.prune.in
#   - chr1_ldpruned.prune.out
#   - QC counts and logs
# ============================================================

workflow AoU_srWGS_chr1_LDPruning_Pilot {

  input {
    File keep_list

    File plink2_file =
      "gs://cloned-all-circadian-patients-for-fitbit-and-prs-wb-perky-onion/plink2"

    File pgen_file =
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr1.pgen"

    File pvar_file =
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr1.pvar"

    File psam_file =
      "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr1.psam"

    Int window_kb = 500
    Float r2_threshold = 0.2

    String output_prefix = "AoU_srWGS_chr1_LDpruned_500kb_r2_0.2"

    # Initial pilot resources.
    Int task_cpu = 4
    Int task_memory_gb = 64
    Int task_disk_gb = 300
    Int plink_memory_mb = 48000
  }

  call PruneExtractChr1 {
    input:
      pgen_file = pgen_file,
      pvar_file = pvar_file,
      psam_file = psam_file,
      keep_list = keep_list,
      plink2_file = plink2_file,
      window_kb = window_kb,
      r2_threshold = r2_threshold,
      output_prefix = output_prefix,
      plink_memory_mb = plink_memory_mb,
      task_cpu = task_cpu,
      task_memory_gb = task_memory_gb,
      task_disk_gb = task_disk_gb
  }

  output {
    File pruned_pgen = PruneExtractChr1.pruned_pgen
    File pruned_pvar = PruneExtractChr1.pruned_pvar
    File pruned_psam = PruneExtractChr1.pruned_psam

    File prune_in = PruneExtractChr1.prune_in
    File prune_out = PruneExtractChr1.prune_out

    File variant_counts = PruneExtractChr1.variant_counts
    File pruning_log = PruneExtractChr1.pruning_log
    File extraction_log = PruneExtractChr1.extraction_log
    File task_log = PruneExtractChr1.task_log
  }
}

# ------------------------------------------------------------
# Chromosome 1 pruning and extraction
# ------------------------------------------------------------

task PruneExtractChr1 {

  input {
    File pgen_file
    File pvar_file
    File psam_file
    File keep_list
    File plink2_file

    Int window_kb
    Float r2_threshold
    String output_prefix

    Int plink_memory_mb
    Int task_cpu
    Int task_memory_gb
    Int task_disk_gb
  }

  command <<<
    set -euo pipefail

    exec > >(tee "chr1_task.log") 2>&1

    echo "============================================================"
    echo "AoU v9 srWGS CHROMOSOME 1 LD PRUNING PILOT"
    echo "Start: $(date -u)"
    echo "============================================================"

    mkdir -p input logs qc

    # Cromwell stages all File inputs.
    # Copy the executable locally; link the large genotype files.
    cp "~{plink2_file}" ./plink2
    chmod 750 ./plink2

    prefix="input/acaf_threshold.chr1"
    out_prefix="~{output_prefix}"

    ln -s "~{pgen_file}" "${prefix}.pgen"
    ln -s "~{pvar_file}" "${prefix}.pvar"
    ln -s "~{psam_file}" "${prefix}.psam"

    echo "PLINK2 version:"
    ./plink2 --version

    echo "Staged input files:"
    ls -lh ./plink2 "~{keep_list}" \
      "${prefix}.pgen" \
      "${prefix}.pvar" \
      "${prefix}.psam"

    [[ -s "~{keep_list}" ]] || {
      echo "ERROR: Keep-list is missing or empty" >&2
      exit 1
    }

    [[ -s "${prefix}.pgen" && \
       -s "${prefix}.pvar" && \
       -s "${prefix}.psam" ]] || {
      echo "ERROR: One or more chromosome input files are missing" >&2
      exit 1
    }

    echo ""
    echo "Pruning parameters:"
    echo "Window: ~{window_kb}kb"
    echo "R2 threshold: ~{r2_threshold}"
    echo "PLINK memory: ~{plink_memory_mb} MB"
    echo "Threads: ~{task_cpu}"
    echo ""

    # --------------------------------------------------------
    # Step 1: Generate LD-pruned variant lists.
    # --------------------------------------------------------

    echo "Starting LD pruning: $(date -u)"

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

    [[ -s "${prune_in}" ]] || {
      echo "ERROR: Missing or empty prune.in file" >&2
      tail -100 "logs/${out_prefix}_pruning.log" >&2 || true
      exit 1
    }

    [[ -f "${prune_out}" ]] || {
      echo "ERROR: Missing prune.out file" >&2
      exit 1
    }

    n_prune_in=$(wc -l < "${prune_in}" | tr -d '[:space:]')
    n_prune_out=$(wc -l < "${prune_out}" | tr -d '[:space:]')

    echo "Pruning completed: $(date -u)"
    echo "Retained variants: ${n_prune_in}"
    echo "Excluded variants: ${n_prune_out}"

    # --------------------------------------------------------
    # Step 2: Extract retained variants and selected samples.
    # --------------------------------------------------------

    echo ""
    echo "Starting genotype extraction: $(date -u)"

    ./plink2 \
      --pfile "${prefix}" \
      --keep "~{keep_list}" \
      --extract "${prune_in}" \
      --make-pgen \
      --memory "~{plink_memory_mb}" \
      --threads "~{task_cpu}" \
      --out "${out_prefix}" \
      > "logs/${out_prefix}_extraction.log" 2>&1

    [[ -s "${out_prefix}.pgen" && \
       -s "${out_prefix}.pvar" && \
       -s "${out_prefix}.psam" ]] || {
      echo "ERROR: Incomplete extracted PGEN trio" >&2
      tail -100 "logs/${out_prefix}_extraction.log" >&2 || true
      exit 1
    }

    # --------------------------------------------------------
    # Step 3: QC counts.
    # --------------------------------------------------------

    n_extracted=$(awk 'BEGIN{n=0} !/^#/ {n++} END{print n}' \
      "${out_prefix}.pvar")

    n_samples=$(awk 'NR>1 {n++} END{print n+0}' \
      "${out_prefix}.psam")

    printf 'chromosome\twindow_kb\tr2_threshold\tprune_in_count\tprune_out_count\textracted_variant_count\tsample_count\n' \
      > "qc/chr1_variant_counts.tsv"

    printf '1\t~{window_kb}\t~{r2_threshold}\t%s\t%s\t%s\t%s\n' \
      "${n_prune_in}" \
      "${n_prune_out}" \
      "${n_extracted}" \
      "${n_samples}" \
      >> "qc/chr1_variant_counts.tsv"

    echo ""
    echo "============================================================"
    echo "CHROMOSOME 1 PILOT COMPLETED"
    echo "============================================================"
    echo "Retained variants: ${n_prune_in}"
    echo "Excluded variants: ${n_prune_out}"
    echo "Extracted variants: ${n_extracted}"
    echo "Samples: ${n_samples}"
    echo "End: $(date -u)"

    echo ""
    echo "Output file sizes:"
    ls -lh \
      "${out_prefix}.pgen" \
      "${out_prefix}.pvar" \
      "${out_prefix}.psam" \
      "${prune_in}" \
      "${prune_out}"

    # Remove only the local symlinks.
    rm -f "${prefix}.pgen" \
          "${prefix}.pvar" \
          "${prefix}.psam"

  >>>

  output {
    File pruned_pgen = "~{output_prefix}.pgen"
    File pruned_pvar = "~{output_prefix}.pvar"
    File pruned_psam = "~{output_prefix}.psam"

    File prune_in = "~{output_prefix}.prune.in"
    File prune_out = "~{output_prefix}.prune.out"

    File variant_counts = "qc/chr1_variant_counts.tsv"

    File pruning_log =
      "logs/~{output_prefix}_pruning.log"

    File extraction_log =
      "logs/~{output_prefix}_extraction.log"

    File task_log = "chr1_task.log"
  }

  runtime {
    cpu: task_cpu
    memory: task_memory_gb + " GB"
    disks: "local-disk " + task_disk_gb + " SSD"
  }
}
