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
        #
        # This remains a SINGLE File input so the WDL can be reused
        # with any PRS weight file.
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

        # Single reusable PRS weight file
        File prs_file

        Int cpu
        Int mem
    }

    command <<<

        set -euo pipefail

        # ============================================================
        # Persistent task log
        # ============================================================

        exec > "chr~{chromosome}.log" 2>&1

        echo ""
        echo "============================================================"
        echo "AoU CHRONOTYPE PRS"
        echo "============================================================"
        echo "Chromosome: ~{chromosome}"
        echo "============================================================"
        echo ""

        # ============================================================
        # Localized input files
        # ============================================================

        echo "Localized PLINK files:"
        echo "BED: ~{bedfile}"
        echo "BIM: ~{bimfile}"
        echo "FAM: ~{famfile}"
        echo ""

        echo "PRS weight file:"
        echo "~{prs_file}"
        echo ""

        echo "Working directory:"
        pwd
        echo ""

        echo "Staged files:"
        find . -maxdepth 6 -type f -ls

        # ============================================================
        # Check PLINK files
        # ============================================================

        echo ""
        echo "Checking PLINK files..."

        ls -lh \
            ~{bedfile} \
            ~{bimfile} \
            ~{famfile}

        # ============================================================
        # Check PRS file
        # ============================================================

        echo ""
        echo "Checking PRS weight file..."

        ls -lh ~{prs_file}

        # ============================================================
        # Check PLINK2
        # ============================================================

        echo ""
        echo "PLINK2 version:"
        plink2 --version

        # ============================================================
        # Count chromosome variants
        # ============================================================

        echo ""
        echo "Chromosome ~{chromosome} variant count:"

        wc -l ~{bimfile}

        # ============================================================
        # Inspect PRS weight file
        # ============================================================

        echo ""
        echo "PRS weight file type:"
        file ~{prs_file}

        echo ""
        echo "First 5 lines of PRS weight file:"

        zcat ~{prs_file} | head -5

        # ============================================================
        # Run PLINK2 scoring
        # ============================================================

        echo ""
        echo "============================================================"
        echo "Starting PLINK2 scoring"
        echo "============================================================"
        echo ""

        plink2 \
            --bfile ~{sub(bedfile, "\.bed$", "")} \
            --score ~{prs_file} \
                1 2 3 \
                header-read \
                cols=+scoresums \
            --out chr~{chromosome}

        # ============================================================
        # Check output
        # ============================================================

        echo ""
        echo "============================================================"
        echo "CHR ~{chromosome} COMPLETE"
        echo "============================================================"
        echo ""

        ls -lh chr~{chromosome}*

        echo ""
        echo "First 10 lines of score file:"

        head -10 chr~{chromosome}.sscore

        echo ""
        echo "============================================================"
        echo "END CHROMOSOME ~{chromosome}"
        echo "============================================================"

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