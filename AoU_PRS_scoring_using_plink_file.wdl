version 1.0

workflow AoU_chronotype_PRS {

input {

    # ============================================================
    # Existing AoU v9 PLINK BED files
    # ============================================================

    Array[File] bedfiles = [
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr1.bed",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr2.bed",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr3.bed",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr4.bed",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr5.bed",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr6.bed",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr7.bed",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr8.bed",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr9.bed",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr10.bed",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr11.bed",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr12.bed",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr13.bed",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr14.bed",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr15.bed",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr16.bed",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr17.bed",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr18.bed",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr19.bed",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr20.bed",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr21.bed",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr22.bed"
    ]

    # ============================================================
    # Existing AoU v9 PLINK BIM files
    # ============================================================

    Array[File] bimfiles = [
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr1.bim",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr2.bim",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr3.bim",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr4.bim",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr5.bim",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr6.bim",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr7.bim",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr8.bim",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr9.bim",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr10.bim",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr11.bim",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr12.bim",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr13.bim",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr14.bim",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr15.bim",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr16.bim",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr17.bim",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr18.bim",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr19.bim",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr20.bim",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr21.bim",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr22.bim"
    ]

    # ============================================================
    # Existing AoU v9 PLINK FAM files
    # ============================================================

    Array[File] famfiles = [
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr1.fam",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr2.fam",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr3.fam",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr4.fam",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr5.fam",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr6.fam",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr7.fam",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr8.fam",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr9.fam",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr10.fam",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr11.fam",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr12.fam",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr13.fam",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr14.fam",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr15.fam",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr16.fam",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr17.fam",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr18.fam",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr19.fam",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr20.fam",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr21.fam",
        "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold.chr22.fam"
    ]

    # ============================================================
    # PRS WEIGHT FILE
    # ============================================================

    File prs_file

    # ============================================================
    # Resources
    # ============================================================

    Int cpu = 16
    Int mem = 120
}

# ================================================================
# Scatter chromosomes 1-22
# ================================================================

scatter (i in range(22)) {

    call RunChromosomePRS {
        input:
            chromosome = i + 1,
            bedfile = bedfiles[i],
            bimfile = bimfiles[i],
            famfile = famfiles[i],
            prs_file = prs_file,
            cpu = cpu,
            mem = mem
    }
}

# ================================================================
# Workflow outputs
# ================================================================

output {
    Array[File] sscore = RunChromosomePRS.sscore
    Array[File] log = RunChromosomePRS.log
}

}

# ====================================================================

# TASK

# ====================================================================

task RunChromosomePRS {

input {

    Int chromosome

    File bedfile
    File bimfile
    File famfile

    File prs_file

    Int cpu
    Int mem
}

command <<<

    set -euo pipefail

    # ============================================================
    # Basic configuration
    # ============================================================

    LOG="chr~{chromosome}.log"

    # GCS location for continuously updated diagnostic log.
    # The timestamp/PID makes simultaneous runs distinct.
    DEBUG_ID="chr~{chromosome}_$(date -u +%Y%m%dT%H%M%SZ)_$$"

    DEBUG_BUCKET="gs://cloned-all-circadian-patients-for-fitbit-and-prs-wb-perky-onion/PRS_diagnostics/${DEBUG_ID}"

    # ============================================================
    # Capture ALL stdout/stderr in the task log
    # ============================================================

    exec > >(tee -a "$LOG") 2>&1

    echo ""
    echo "============================================================"
    echo "AoU CHRONOTYPE PRS DIAGNOSTIC RUN"
    echo "============================================================"
    echo "Chromosome: ~{chromosome}"
    echo "Start time: $(date -u)"
    echo "Hostname: $(hostname)"
    echo "PID: $$"
    echo "CPU requested: ~{cpu}"
    echo "Memory requested: ~{mem} GB"
    echo "Diagnostic GCS location: $DEBUG_BUCKET"
    echo "============================================================"
    echo ""

    # ============================================================
    # Function to upload current log to GCS
    #
    # This is important because Cromwell only delocalizes outputs
    # after the task finishes. If the VM is killed, this gives us
    # a copy of the log from shortly before termination.
    # ============================================================

    upload_diagnostic_log() {
        gsutil -q cp "$LOG" "${DEBUG_BUCKET}/chr~{chromosome}.log" \
            2>/dev/null || true
    }

    # ============================================================
    # Upload log periodically while task is running
    # ============================================================

    (
        while true; do
            sleep 30
            upload_diagnostic_log
        done
    ) &

    MONITOR_UPLOAD_PID=$!

    # Make sure the background uploader is stopped when the task
    # exits normally or receives a catchable termination signal.
    cleanup() {
        echo ""
        echo "============================================================"
        echo "CLEANUP"
        echo "============================================================"
        echo "Time: $(date -u)"

        kill "$MONITOR_UPLOAD_PID" 2>/dev/null || true

        echo ""
        echo "Final memory:"
        free -h || true

        echo ""
        echo "Final disk:"
        df -h || true

        echo ""
        echo "Final process list:"
        ps aux --sort=-%mem | head -20 || true

        upload_diagnostic_log || true

        echo ""
        echo "Diagnostic log location:"
        echo "${DEBUG_BUCKET}/chr~{chromosome}.log"
        echo ""
    }

    trap cleanup EXIT

    # ============================================================
    # SYSTEM INFORMATION
    # ============================================================

    echo ""
    echo "============================================================"
    echo "SYSTEM INFORMATION"
    echo "============================================================"

    echo ""
    echo "--- Date ---"
    date -u

    echo ""
    echo "--- Host ---"
    hostname
    uname -a

    echo ""
    echo "--- CPU ---"
    echo "nproc:"
    nproc

    echo ""
    echo "CPU information:"
    grep -E '^(processor|model name|cpu MHz)' /proc/cpuinfo | head -100 || true

    echo ""
    echo "--- MEMORY ---"
    free -h

    echo ""
    echo "Selected /proc/meminfo:"
    grep -E '^(MemTotal|MemFree|MemAvailable|Buffers|Cached|SwapTotal|SwapFree)' /proc/meminfo

    echo ""
    echo "--- DISK ---"
    df -h

    echo ""
    echo "--- DISK /mnt/disks/cromwell_root ---"
    df -h /mnt/disks/cromwell_root

    # ============================================================
    # WORKING DIRECTORY
    # ============================================================

    echo ""
    echo "============================================================"
    echo "WORKING DIRECTORY"
    echo "============================================================"

    pwd

    echo ""
    echo "Directory contents:"
    ls -lah

    echo ""
    echo "All staged files:"
    find . -maxdepth 6 -type f -ls

    # ============================================================
    # INPUT FILES
    # ============================================================

    echo ""
    echo "============================================================"
    echo "INPUT FILES"
    echo "============================================================"

    echo ""
    echo "BED:"
    echo "~{bedfile}"
    ls -lh ~{bedfile}

    echo ""
    echo "BIM:"
    echo "~{bimfile}"
    ls -lh ~{bimfile}

    echo ""
    echo "FAM:"
    echo "~{famfile}"
    ls -lh ~{famfile}

    echo ""
    echo "PRS:"
    echo "~{prs_file}"
    ls -lh ~{prs_file}

    echo ""
    echo "Input file checksums:"
    sha256sum \
        ~{bedfile} \
        ~{bimfile} \
        ~{famfile} \
        ~{prs_file}

    # ============================================================
    # PLINK2
    # ============================================================

    echo ""
    echo "============================================================"
    echo "PLINK2 INFORMATION"
    echo "============================================================"

    echo ""
    echo "PLINK2 version:"
    plink2 --version

    echo ""
    echo "PLINK2 location:"
    which plink2

    # ============================================================
    # BIM / FAM / PRS INSPECTION
    # ============================================================

    echo ""
    echo "============================================================"
    echo "INPUT CONTENT CHECKS"
    echo "============================================================"

    echo ""
    echo "BIM variant count:"
    wc -l ~{bimfile}

    echo ""
    echo "FAM sample count:"
    wc -l ~{famfile}

    echo ""
    echo "PRS file type:"
    file ~{prs_file}

    echo ""
    echo "PRS compressed file size:"
    ls -lh ~{prs_file}

    echo ""
    echo "First 5 lines of PRS file:"
    zcat ~{prs_file} | head -5 || true

    echo ""
    echo "PRS row count:"
    zcat ~{prs_file} | wc -l

    # ============================================================
    # PRE-PLINK RESOURCE SNAPSHOT
    # ============================================================

    echo ""
    echo "============================================================"
    echo "PRE-PLINK RESOURCE SNAPSHOT"
    echo "============================================================"

    date -u

    echo ""
    echo "Memory:"
    free -h

    echo ""
    echo "Disk:"
    df -h /mnt/disks/cromwell_root

    echo ""
    echo "Processes:"
    ps aux --sort=-%mem | head -20 || true

    upload_diagnostic_log

    # ============================================================
    # RESOURCE MONITOR
    #
    # Runs independently while PLINK is executing.
    # Every 30 seconds we record:
    #   - memory
    #   - disk
    #   - CPU/process information
    #
    # This is the critical diagnostic section.
    # ============================================================

    (
        while true; do

            echo ""
            echo "############################################################"
            echo "RESOURCE HEARTBEAT"
            echo "Time: $(date -u)"
            echo "############################################################"

            echo ""
            echo "--- Memory ---"
            free -h

            echo ""
            echo "--- Disk ---"
            df -h /mnt/disks/cromwell_root

            echo ""
            echo "--- Load ---"
            uptime

            echo ""
            echo "--- PLINK processes ---"
            ps -eo pid,ppid,%cpu,%mem,rss,vsz,etime,cmd \
                --sort=-%mem | head -20 || true

            echo ""
            echo "--- Process tree ---"
            ps aux --forest | head -50 || true

        done
    ) &

    RESOURCE_MONITOR_PID=$!

    # ============================================================
    # RUN PLINK2
    # ============================================================

    echo ""
    echo "============================================================"
    echo "STARTING PLINK2 SCORING"
    echo "============================================================"
    echo "Start time: $(date -u)"
    echo "Threads: ~{cpu}"
    echo "Memory: ~{mem}000 MB"
    echo ""

    set +e

    plink2 \
        --memory ~{mem}000 \
        --threads ~{cpu} \
        --bfile ~{sub(bedfile, "\.bed$", "")} \
        --score ~{prs_file} \
            1 2 3 \
            header-read \
            cols=+scoresums \
        --out chr~{chromosome}

    PLINK_RC=$?

    set -e

    # ============================================================
    # STOP RESOURCE MONITOR
    # ============================================================

    kill "$RESOURCE_MONITOR_PID" 2>/dev/null || true

    echo ""
    echo "============================================================"
    echo "PLINK2 FINISHED"
    echo "============================================================"
    echo "End time: $(date -u)"
    echo "PLINK exit code: $PLINK_RC"
    echo ""

    # ============================================================
    # POST-PLINK RESOURCE SNAPSHOT
    # ============================================================

    echo "============================================================"
    echo "POST-PLINK RESOURCE SNAPSHOT"
    echo "============================================================"

    echo ""
    echo "--- Memory ---"
    free -h

    echo ""
    echo "--- Disk ---"
    df -h /mnt/disks/cromwell_root

    echo ""
    echo "--- Processes ---"
    ps aux --sort=-%mem | head -20 || true

    # ============================================================
    # CHECK OUTPUT
    # ============================================================

    echo ""
    echo "============================================================"
    echo "OUTPUT CHECK"
    echo "============================================================"

    if [[ -f "chr~{chromosome}.sscore" ]]; then
        echo "SSCORE EXISTS"
        ls -lh "chr~{chromosome}.sscore"

        echo ""
        echo "SSCORE line count:"
        wc -l "chr~{chromosome}.sscore"

        echo ""
        echo "First 10 lines:"
        head -10 "chr~{chromosome}.sscore"
    else
        echo "WARNING: chr~{chromosome}.sscore DOES NOT EXIST"
    fi

    echo ""
    echo "All chromosome output files:"
    ls -lh chr~{chromosome}* 2>/dev/null || true

    # ============================================================
    # FINAL DIAGNOSTIC UPLOAD
    # ============================================================

    echo ""
    echo "============================================================"
    echo "FINAL DIAGNOSTIC UPLOAD"
    echo "============================================================"

    upload_diagnostic_log

    echo ""
    echo "Diagnostic log:"
    echo "${DEBUG_BUCKET}/chr~{chromosome}.log"

    echo ""
    echo "============================================================"
    echo "END CHROMOSOME ~{chromosome}"
    echo "============================================================"

    # Preserve the actual PLINK exit code.
    exit "$PLINK_RC"

>>>

output {

    File sscore = "chr" + chromosome + ".sscore"

    File log = "chr" + chromosome + ".log"
}

runtime {

    docker:
        "ghcr.io/chrchang/plink-ng:latest"

    cpu:
        cpu

    memory:
        mem + " GB"

    disks:
        "local-disk 500 SSD"
}


}
