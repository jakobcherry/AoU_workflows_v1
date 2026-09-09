version 1.0

workflow AoU_chronotype_PRS {

    input {

        # ============================================================
        # Existing AoU PLINK files
        # ============================================================

        String plink_bucket =
            "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold"

        # ============================================================
        # PRS weight file
        # ============================================================

        File prs_file

        # ============================================================
        # Output bucket
        # ============================================================

        String output_bucket =
            "gs://cloned-all-circadian-patients-for-fitbit-and-prs-wb-perky-onion/Chronotype_even_PRS_v1"

        # ============================================================
        # Resources
        # ============================================================

        Int cpu = 16
        Int mem = 120
    }


    # ================================================================
    # CHROMOSOMES 1-22
    # ================================================================

    scatter (chromosome in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]) {

        call RunChromosomePRS {

            input:

                chromosome =
                    chromosome,

                plink_prefix =
                    plink_bucket + ".chr" + chromosome,

                prs_file =
                    prs_file,

                output_prefix =
                    output_bucket + "/chr" + chromosome + "/chr" + chromosome,

                cpu =
                    cpu,

                mem =
                    mem
        }
    }


    output {

        Array[File] sscore =
            RunChromosomePRS.sscore

    }
}


# ====================================================================
# TASK
# ====================================================================

task RunChromosomePRS {

    input {

        Int chromosome

        String plink_prefix

        File prs_file

        String output_prefix

        Int cpu

        Int mem
    }


    # =================================================================
    # COMMAND
    # =================================================================

    command <<<
        set -euo pipefail

        echo ""
        echo "============================================================"
        echo "AoU CHRONOTYPE PRS"
        echo "============================================================"
        echo "Chromosome: ~{chromosome}"
        echo ""
        echo "PLINK prefix:"
        echo "~{plink_prefix}"
        echo ""
        echo "PRS file:"
        echo "~{prs_file}"
        echo ""
        echo "Output:"
        echo "~{output_prefix}"
        echo "============================================================"
        echo ""


        # ============================================================
        # Check PLINK installation
        # ============================================================

        plink2 --version


        # ============================================================
        # Check input files
        # ============================================================

        echo "Checking PLINK files..."

        ls -lh \
            ~{plink_prefix}.bed \
            ~{plink_prefix}.bim \
            ~{plink_prefix}.fam


        echo ""
        echo "Checking PRS file..."
        ls -lh ~{prs_file}


        # ============================================================
        # Count variants
        # ============================================================

        echo ""
        echo "PLINK variant count:"

        wc -l ~{plink_prefix}.bim


        # ============================================================
        # Run PRS scoring
        # ============================================================

        echo ""
        echo "Starting PLINK2 scoring..."
        echo ""

        plink2 \
            --bfile ~{plink_prefix} \
            --score ~{prs_file} \
                1 2 3 \
                header-read \
                cols=+scoresums \
            --out ~{output_prefix}


        # ============================================================
        # Check output
        # ============================================================

        echo ""
        echo "============================================================"
        echo "CHR ~{chromosome} COMPLETE"
        echo "============================================================"

        ls -lh ~{output_prefix}*


        echo ""
        echo "First 10 lines of score file:"
        head -10 ~{output_prefix}.sscore

    >>>


    # =================================================================
    # OUTPUT
    # =================================================================

    output {

        File sscore =
            output_prefix + ".sscore"

    }


    # =================================================================
    # RUNTIME
    # =================================================================

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