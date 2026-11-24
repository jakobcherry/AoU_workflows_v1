
version 1.0

workflow AoU_srWGS_chr1_LDPruning_ListsOnly {

  input {
    File keep_list

    File plink2_file = "gs://cloned-all-circadian-patients-for-fitbit-and-prs-wb-perky-onion/plink2"

    File pgen_file = "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr1.pgen"

    File pvar_file = "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr1.pvar"

    File psam_file = "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr1.psam"

    Int window_kb = 500
    Float r2_threshold = 0.2

    String output_prefix = "AoU_srWGS_chr1_LDprune_500kb_r2_0.2"

    Int task_cpu = 4
    Int task_memory_gb = 64
    Int task_disk_gb = 300
    Int plink_memory_mb = 48000
  }

  call GenerateChr1PruneLists {
    input:
      keep_list = keep_list,
      plink2_file = plink2_file,
      pgen_file = pgen_file,
      pvar_file = pvar_file,
      psam_file = psam_file,
      window_kb = window_kb,
      r2_threshold = r2_threshold,
      output_prefix = output_prefix,
      task_cpu = task_cpu,
      task_memory_gb = task_memory_gb,
      task_disk_gb = task_disk_gb,
      plink_memory_mb = plink_memory_mb
  }

  output {
    File prune_in = GenerateChr1PruneLists.prune_in
    File prune_out = GenerateChr1PruneLists.prune_out
    File variant_counts = GenerateChr1PruneLists.variant_counts
    File pruning_log = GenerateChr1PruneLists.pruning_log
    File task_log = GenerateChr1PruneLists.task_log
  }
}

task GenerateChr1PruneLists {

  input {
    File keep_list
    File plink2_file
    File pgen_file
    File pvar_file
    File psam_file

    Int window_kb
    Float r2_threshold
    String output_prefix

    Int task_cpu
    Int task_memory_gb
    Int task_disk_gb
    Int plink_memory_mb
  }

  command <<<
    set -euo pipefail

    exec > >(tee "chr1_pruning_task.log") 2>&1

    echo "============================================================"
    echo "AoU v9 srWGS CHR1 LD PRUNING"
    echo "PRUNE LISTS AND COUNTS ONLY"
    echo "Start: $(date -u)"
    echo "============================================================"

    mkdir -p input logs qc

    cp "~{plink2_file}" ./plink2
    chmod 750 ./plink2

    prefix="input/acaf_threshold.chr1"
    out_prefix="~{output_prefix}"

    ln -s "~{pgen_file}" "${prefix}.pgen"
    ln -s "~{pvar_file}" "${prefix}.pvar"
    ln -s "~{psam_file}" "${prefix}.psam"

    echo "PLINK2 version:"
    ./plink2 --version

    echo "Input files:"
    ls -lh ./plink2 "~{keep_list}" \
      "${prefix}.pgen" "${prefix}.pvar" "${prefix}.psam"

    [[ -s "~{keep_list}" ]] || {
      echo "ERROR: Keep-list missing or empty" >&2
      exit 1
    }

    [[ -s "${prefix}.pgen" && \
       -s "${prefix}.pvar" && \
       -s "${prefix}.psam" ]] || {
      echo "ERROR: Missing chromosome input file" >&2
      exit 1
    }

    echo ""
    echo "LD pruning parameters:"
    echo "Window: ~{window_kb} kb"
    echo "R2 threshold: ~{r2_threshold}"
    echo "Threads: ~{task_cpu}"
    echo "PLINK memory: ~{plink_memory_mb} MB"

    # Generate .prune.in and .prune.out only.
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
      echo "ERROR: Missing or empty prune.in" >&2
      tail -100 "logs/${out_prefix}_pruning.log" >&2 || true
      exit 1
    }

    [[ -f "${prune_out}" ]] || {
      echo "ERROR: Missing prune.out" >&2
      tail -100 "logs/${out_prefix}_pruning.log" >&2 || true
      exit 1
    }

    n_prune_in=$(wc -l < "${prune_in}" | tr -d '[:space:]')
    n_prune_out=$(wc -l < "${prune_out}" | tr -d '[:space:]')

    n_samples=$(awk 'NR>1 {n++} END{print n+0}' \
      "${prefix}.psam")

    printf 'chromosome\twindow_kb\tr2_threshold\tprune_in_count\tprune_out_count\tsample_count\n' \
      > "qc/chr1_variant_counts.tsv"

    printf '1\t~{window_kb}\t~{r2_threshold}\t%s\t%s\t%s\n' \
      "${n_prune_in}" \
      "${n_prune_out}" \
      "${n_samples}" \
      >> "qc/chr1_variant_counts.tsv"

    echo ""
    echo "============================================================"
    echo "CHR1 LD PRUNING COMPLETED"
    echo "============================================================"
    echo "Retained variants: ${n_prune_in}"
    echo "Excluded variants: ${n_prune_out}"
    echo "Samples in input PSAM: ${n_samples}"
    echo "End: $(date -u)"

    echo "Output files:"
    ls -lh "${prune_in}" "${prune_out}" \
      "qc/chr1_variant_counts.tsv"

    rm -f "${prefix}.pgen" \
          "${prefix}.pvar" \
          "${prefix}.psam"

  >>>

  output {
    File prune_in = "~{output_prefix}.prune.in"
    File prune_out = "~{output_prefix}.prune.out"
    File variant_counts = "qc/chr1_variant_counts.tsv"
    File pruning_log = "logs/~{output_prefix}_pruning.log"
    File task_log = "chr1_pruning_task.log"
  }

  runtime {
    cpu: task_cpu
    memory: task_memory_gb + " GB"
    disks: "local-disk " + task_disk_gb + " SSD"
  }
}
