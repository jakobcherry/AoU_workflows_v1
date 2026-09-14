version 1.0

# ============================================================
# GENERIC AoU PRS SCORING WORKFLOW
#
# Inputs:
#   1. Original PRS CSV/bgzip format:
#      chr,bp,rs_number,effect_allele,noneffect_allele,weight
#
#   2. PLINK2-compatible PRS format:
#      ID EFFECT_ALLELE WEIGHT
#
# The workflow:
#   - scatters across chromosomes 1-22
#   - localizes the AoU v9 PLINK BED/BIM/FAM files
#   - matches the original PRS to AoU variants using:
#       chromosome + position + unordered allele pair
#   - creates an AoU-compatible ID/EFFECT_ALLELE/WEIGHT file
#   - scores directly against the original AoU BED
#   - does NOT create a temporary PRS-specific BED
#
# Outputs:
#   - one .sscore per chromosome
#   - one AoU-compatible score file per chromosome
#   - one variant matching/QC table per chromosome
#   - one extraction list per chromosome
#   - one matching summary per chromosome
#   - PLINK2 log per chromosome
#   - task log per chromosome
# ============================================================


workflow AoU_PRS {

    input {

        # ========================================================
        # PLINK2 executable
        # ========================================================

        File plink2_file = "gs://cloned-all-circadian-patients-for-fitbit-and-prs-wb-perky-onion/plink2"


        # ========================================================
        # PRS input
        # ========================================================

        File prs_file


        # ========================================================
        # AoU v9 PLINK BED files
        # ========================================================

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


        # ========================================================
        # AoU v9 PLINK BIM files
        # ========================================================

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


        # ========================================================
        # AoU v9 PLINK FAM files
        # ========================================================

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


        # ========================================================
        # Resources
        # ========================================================

        Int cpu = 16
        Int mem = 120
    }


    # ============================================================
    # SCATTER ACROSS CHROMOSOMES 1-22
    # ============================================================

    scatter (i in range(22)) {

        call RunChromosomePRS {

            input:
                chromosome = i + 1,
                prs_file = prs_file,
                plink2_file = plink2_file,
                bedfile = bedfiles[i],
                bimfile = bimfiles[i],
                famfile = famfiles[i],
                cpu = cpu,
                mem = mem
        }
    }


    # ============================================================
    # WORKFLOW OUTPUTS
    # ============================================================

    output {

        Array[File] sscore =
            RunChromosomePRS.sscore

        Array[File] score_file =
            RunChromosomePRS.score_file

        Array[File] extract_file =
            RunChromosomePRS.extract_file

        Array[File] match_table =
            RunChromosomePRS.match_table

        Array[File] match_summary =
            RunChromosomePRS.match_summary

        Array[File] log =
            RunChromosomePRS.log

        Array[File] task_log =
            RunChromosomePRS.task_log
    }
}


# ==================================================================
# CHROMOSOME PRS TASK
# ==================================================================

task RunChromosomePRS {

    input {

        Int chromosome

        File prs_file

        File plink2_file

        File bedfile
        File bimfile
        File famfile

        Int cpu
        Int mem
    }


    # ==============================================================
    # COMMAND
    # ==============================================================

    command <<<
        set -euo pipefail

        # ==========================================================
        # TASK LOG
        # ==========================================================

        exec > >(tee "chr~{chromosome}_task.log") 2>&1

        echo "============================================================"
        echo "AoU PRS SCORING"
        echo "============================================================"
        echo "Chromosome: ~{chromosome}"
        echo "Start:      $(date -u)"
        echo ""


        # ==========================================================
        # STAGED FILES
        # ==========================================================

        echo "============================================================"
        echo "STAGED FILES"
        echo "============================================================"

        find . -maxdepth 2 -type f -printf '%p\t%k KB\n' | sort

        echo ""
        echo "BED:"
        ls -lh ~{bedfile}

        echo ""
        echo "BIM:"
        ls -lh ~{bimfile}

        echo ""
        echo "FAM:"
        ls -lh ~{famfile}

        echo ""
        echo "PRS:"
        ls -lh ~{prs_file}

        echo ""


        # ==========================================================
        # PLINK2
        # ==========================================================

        chmod 750 ~{plink2_file}

        echo "============================================================"
        echo "PLINK2 VERSION"
        echo "============================================================"

        ~{plink2_file} --version

        echo ""


        # ==========================================================
        # BASIC AoU DATASET COUNTS
        # ==========================================================

        BIM_N=$(wc -l < ~{bimfile})
        FAM_N=$(wc -l < ~{famfile})

        echo "AoU BIM variants: $BIM_N"
        echo "AoU FAM samples:  $FAM_N"
        echo ""


        # ==========================================================
        # DECOMPRESS PRS IF NECESSARY
        #
        # This supports both:
        #   .gz / .bgz
        #   uncompressed text
        # ==========================================================

        echo "============================================================"
        echo "PREPARING PRS INPUT"
        echo "============================================================"

        if gzip -t ~{prs_file} >/dev/null 2>&1; then
            echo "PRS input is gzip/bgzip compressed."
            gzip -dc ~{prs_file} > prs_input.txt
        else
            echo "PRS input is not gzip compressed."
            cp ~{prs_file} prs_input.txt
        fi

        FIRST_LINE=$(head -1 prs_input.txt | tr -d '\r')

        echo "PRS header:"
        echo "$FIRST_LINE"
        echo ""


        # ==========================================================
        # DETECT PRS FORMAT
        # ==========================================================

        PRS_FORMAT=""

        if [[ "$FIRST_LINE" == "chr,bp,rs_number,effect_allele,noneffect_allele,weight" ]]; then

            PRS_FORMAT="six_column"

        elif printf '%s\n' "$FIRST_LINE" |
            awk 'BEGIN {
                FS="[,\t ]+"
            }
            {
                if (
                    toupper($1) == "ID" &&
                    toupper($2) == "EFFECT_ALLELE" &&
                    toupper($3) == "WEIGHT"
                ) {
                    exit 0
                }
                exit 1
            }'; then

            PRS_FORMAT="three_column"

        else

            echo "ERROR: PRS format was not recognized."
            echo ""
            echo "Expected either:"
            echo "  chr,bp,rs_number,effect_allele,noneffect_allele,weight"
            echo ""
            echo "or:"
            echo "  ID EFFECT_ALLELE WEIGHT"
            echo ""
            exit 1
        fi

        echo "Detected PRS format: $PRS_FORMAT"
        echo ""


        # ==========================================================
        # PREPARE CHROMOSOME-SPECIFIC PRS
        # ==========================================================

        rm -f prs_chr.tsv
        rm -f "chr~{chromosome}_aou_prs_score.tsv"
        rm -f "chr~{chromosome}_prs_extract.txt"
        rm -f "chr~{chromosome}_prs_aou_matches.tsv"
        rm -f "chr~{chromosome}_prs_match_summary.tsv"


        if [[ "$PRS_FORMAT" == "six_column" ]]; then

            echo "Extracting chromosome ~{chromosome} from six-column PRS."

            awk -F',' -v OFS='\t' -v target="~{chromosome}" '
                NR == 1 {
                    next
                }

                {
                    chr = $1
                    gsub(/^chr/, "", chr)

                    if (chr == target) {
                        print chr, $2, $3, $4, $5, $6
                    }
                }
            ' prs_input.txt > prs_chr.tsv


        elif [[ "$PRS_FORMAT" == "three_column" ]]; then

            echo "Preparing three-column PRS."

            awk -F'[,\t ]+' -v OFS='\t' '
                NR == 1 {
                    next
                }

                NF >= 3 {
                    print $1, $2, $3
                }
            ' prs_input.txt > prs_chr.tsv

        fi


        PRS_N=$(wc -l < prs_chr.tsv)

        echo ""
        echo "PRS records for chromosome ~{chromosome}: $PRS_N"
        echo ""


        if [[ "$PRS_N" -eq 0 ]]; then
            echo "ERROR: No PRS records found for chromosome ~{chromosome}."
            exit 1
        fi


        # ==========================================================
        # MATCH ORIGINAL SIX-COLUMN PRS TO AoU BIM
        #
        # Matching key:
        #   base-pair position + unordered allele pair
        #
        # This allows the AoU A1/A2 ordering to differ from the PRS
        # while preserving the PRS effect allele.
        # ==========================================================

        if [[ "$PRS_FORMAT" == "six_column" ]]; then

            echo "============================================================"
            echo "MATCHING PRS TO AoU BIM"
            echo "============================================================"

            awk -F'\t' -v OFS='\t' \
                -v target="~{chromosome}" \
                -v score_file="chr~{chromosome}_aou_prs_score.tsv" \
                -v extract_file="chr~{chromosome}_prs_extract.txt" \
                -v match_file="chr~{chromosome}_prs_aou_matches.tsv" \
                -v summary_file="chr~{chromosome}_prs_match_summary.tsv" '

                function allele_key(a, b) {
                    a = toupper(a)
                    b = toupper(b)

                    if (a <= b) {
                        return a ":" b
                    }

                    return b ":" a
                }


                BEGIN {
                    print "ID", "EFFECT_ALLELE", "WEIGHT" > score_file

                    print "ID" > extract_file

                    print \
                        "chr", \
                        "bp", \
                        "rs_number", \
                        "effect_allele", \
                        "noneffect_allele", \
                        "weight", \
                        "aou_bim_id", \
                        "aou_a1", \
                        "aou_a2", \
                        "effect_orientation" \
                        > match_file
                }


                # ------------------------------------------------
                # FIRST FILE = AoU BIM
                # ------------------------------------------------

                FNR == NR {

                    chr = $1
                    id = $2
                    bp = $4
                    a1 = toupper($5)
                    a2 = toupper($6)

                    gsub(/^chr/, "", chr)

                    if (chr == target) {

                        key = bp ":" allele_key(a1, a2)

                        bim_count[key]++

                        bim_id[key] = id
                        bim_a1[key] = a1
                        bim_a2[key] = a2
                    }

                    next
                }


                # ------------------------------------------------
                # SECOND FILE = CHROMOSOME PRS
                # ------------------------------------------------

                {

                    chr = $1
                    bp = $2
                    prs_id = $3
                    effect = toupper($4)
                    noneffect = toupper($5)
                    weight = $6

                    gsub(/^chr/, "", chr)

                    if (chr != target) {
                        next
                    }

                    key = bp ":" allele_key(effect, noneffect)

                    if (bim_count[key] == 1) {

                        id = bim_id[key]
                        a1 = bim_a1[key]
                        a2 = bim_a2[key]

                        if (effect == a1) {
                            orientation = "A1"
                            a1_count++
                        }
                        else if (effect == a2) {
                            orientation = "A2"
                            a2_count++
                        }
                        else {
                            allele_mismatch++
                            next
                        }

                        print id, effect, weight >> score_file

                        print id >> extract_file

                        print \
                            chr, \
                            bp, \
                            prs_id, \
                            effect, \
                            noneffect, \
                            weight, \
                            id, \
                            a1, \
                            a2, \
                            orientation \
                            >> match_file

                        matched++

                    }
                    else if (bim_count[key] > 1) {

                        ambiguous++

                    }
                    else {

                        unmatched++
                    }
                }


                END {

                    print \
                        "metric", \
                        "count" \
                        > summary_file

                    print \
                        "PRS_variants", \
                        matched + unmatched + ambiguous \
                        >> summary_file

                    print \
                        "matched", \
                        matched \
                        >> summary_file

                    print \
                        "effect_allele_A1", \
                        a1_count \
                        >> summary_file

                    print \
                        "effect_allele_A2", \
                        a2_count \
                        >> summary_file

                    print \
                        "unmatched", \
                        unmatched \
                        >> summary_file

                    print \
                        "ambiguous", \
                        ambiguous \
                        >> summary_file

                    print \
                        "allele_mismatch", \
                        allele_mismatch \
                        >> summary_file
                }

            ' ~{bimfile} prs_chr.tsv


        # ==========================================================
        # MATCH ALREADY-FORMATTED THREE-COLUMN PRS TO AoU BIM
        # ==========================================================

        else

            echo "============================================================"
            echo "MATCHING FORMATTED PRS TO AoU BIM"
            echo "============================================================"

            awk -F'\t' -v OFS='\t' \
                -v target="~{chromosome}" \
                -v score_file="chr~{chromosome}_aou_prs_score.tsv" \
                -v extract_file="chr~{chromosome}_prs_extract.txt" \
                -v match_file="chr~{chromosome}_prs_aou_matches.tsv" \
                -v summary_file="chr~{chromosome}_prs_match_summary.tsv" '

                function normalize_id(x) {
                    y = x

                    if (y ~ /^chr/) {
                        y = substr(y, 4)
                    }

                    return y
                }


                BEGIN {
                    print "ID", "EFFECT_ALLELE", "WEIGHT" > score_file

                    print "ID" > extract_file

                    print \
                        "input_id", \
                        "aou_bim_id", \
                        "effect_allele", \
                        "aou_a1", \
                        "aou_a2", \
                        "weight", \
                        "effect_orientation" \
                        > match_file
                }


                # ------------------------------------------------
                # FIRST FILE = AoU BIM
                # ------------------------------------------------

                FNR == NR {

                    chr = $1
                    id = $2
                    a1 = toupper($5)
                    a2 = toupper($6)

                    gsub(/^chr/, "", chr)

                    if (chr == target) {

                        normalized = normalize_id(id)

                        bim_id[id] = id
                        bim_id[normalized] = id

                        bim_a1[id] = a1
                        bim_a2[id] = a2

                        bim_count[normalized]++
                    }

                    next
                }


                # ------------------------------------------------
                # SECOND FILE = FORMATTED PRS
                # ------------------------------------------------

                {

                    input_id = $1
                    effect = toupper($2)
                    weight = $3

                    normalized = normalize_id(input_id)

                    if (bim_count[normalized] == 1) {

                        id = bim_id[normalized]

                        a1 = bim_a1[id]
                        a2 = bim_a2[id]

                        if (effect == a1) {

                            orientation = "A1"
                            a1_count++

                        }
                        else if (effect == a2) {

                            orientation = "A2"
                            a2_count++

                        }
                        else {

                            allele_mismatch++
                            next
                        }

                        print id, effect, weight >> score_file

                        print id >> extract_file

                        print \
                            input_id, \
                            id, \
                            effect, \
                            a1, \
                            a2, \
                            weight, \
                            orientation \
                            >> match_file

                        matched++

                    }
                    else if (bim_count[normalized] > 1) {

                        ambiguous++

                    }
                    else {

                        unmatched++
                    }
                }


                END {

                    print \
                        "metric", \
                        "count" \
                        > summary_file

                    print \
                        "PRS_variants", \
                        matched + unmatched + ambiguous \
                        >> summary_file

                    print \
                        "matched", \
                        matched \
                        >> summary_file

                    print \
                        "effect_allele_A1", \
                        a1_count \
                        >> summary_file

                    print \
                        "effect_allele_A2", \
                        a2_count \
                        >> summary_file

                    print \
                        "unmatched", \
                        unmatched \
                        >> summary_file

                    print \
                        "ambiguous", \
                        ambiguous \
                        >> summary_file

                    print \
                        "allele_mismatch", \
                        allele_mismatch \
                        >> summary_file
                }

            ' ~{bimfile} prs_chr.tsv

        fi


        # ==========================================================
        # MATCHING QC
        # ==========================================================

        echo ""
        echo "============================================================"
        echo "MATCHING QC"
        echo "============================================================"

        cat "chr~{chromosome}_prs_match_summary.tsv"

        MATCHED=$(awk -F'\t' '
            $1 == "matched" {
                print $2
            }
        ' "chr~{chromosome}_prs_match_summary.tsv")

        UNMATCHED=$(awk -F'\t' '
            $1 == "unmatched" {
                print $2
            }
        ' "chr~{chromosome}_prs_match_summary.tsv")

        AMBIGUOUS=$(awk -F'\t' '
            $1 == "ambiguous" {
                print $2
            }
        ' "chr~{chromosome}_prs_match_summary.tsv")

        ALLELE_MISMATCH=$(awk -F'\t' '
            $1 == "allele_mismatch" {
                print $2
            }
        ' "chr~{chromosome}_prs_match_summary.tsv")


        echo ""
        echo "Matched:          $MATCHED"
        echo "Unmatched:        $UNMATCHED"
        echo "Ambiguous:        $AMBIGUOUS"
        echo "Allele mismatch:  $ALLELE_MISMATCH"
        echo ""


        if [[ "$MATCHED" -eq 0 ]]; then
            echo "ERROR: No PRS variants matched AoU BIM for chromosome ~{chromosome}."
            exit 1
        fi


        if [[ "$AMBIGUOUS" -gt 0 ]]; then
            echo "ERROR: Ambiguous AoU BIM matches detected."
            echo "Refusing to score ambiguous variants."
            exit 1
        fi


        # ==========================================================
        # CHECK SCORE FILE
        # ==========================================================

        SCORE_N=$(($(wc -l < "chr~{chromosome}_aou_prs_score.tsv") - 1))
        EXTRACT_N=$(wc -l < "chr~{chromosome}_prs_extract.txt")

        echo "Score file variants:    $SCORE_N"
        echo "Extract list variants:  $EXTRACT_N"
        echo "Matched variants:       $MATCHED"
        echo ""


        if [[ "$SCORE_N" -ne "$MATCHED" ]]; then
            echo "ERROR: Score file count does not equal match count."
            exit 1
        fi


        if [[ "$EXTRACT_N" -ne "$MATCHED" ]]; then
            echo "ERROR: Extract list count does not equal match count."
            exit 1
        fi


        # ==========================================================
        # CHECK FOR DUPLICATE AoU IDs
        # ==========================================================

        DUPLICATES=$(sort "chr~{chromosome}_prs_extract.txt" |
            uniq -d |
            wc -l)

        echo "Duplicate AoU IDs:      $DUPLICATES"
        echo ""

        if [[ "$DUPLICATES" -gt 0 ]]; then
            echo "ERROR: Duplicate AoU variant IDs detected."
            exit 1
        fi


        # ==========================================================
        # SHOW FIRST SCORE RECORDS
        # ==========================================================

        echo "First 10 score records:"
        head -10 "chr~{chromosome}_aou_prs_score.tsv"

        echo ""


        # ==========================================================
        # SCORE DIRECTLY AGAINST ORIGINAL AoU BED
        #
        # IMPORTANT:
        # No --make-bed.
        # No temporary BED/BIM/FAM is created.
        # The original localized AoU dataset is used directly.
        # ==========================================================

        echo "============================================================"
        echo "PLINK2 PRS SCORING"
        echo "============================================================"

        echo "Chromosome: ~{chromosome}"
        echo "Variants:   $MATCHED"
        echo "Threads:    ~{cpu}"
        echo "Memory:     ~{mem} GB"
        echo "Start:      $(date -u)"
        echo ""


        ~{plink2_file} \
            --memory ~{mem}000 \
            --threads ~{cpu} \
            --bfile ~{sub(bedfile, "\\.bed$", "")} \
            --extract "chr~{chromosome}_prs_extract.txt" \
            --score "chr~{chromosome}_aou_prs_score.tsv" \
                1 2 3 \
                header-read \
                cols=+scoresums \
            --out "chr~{chromosome}"


        # ==========================================================
        # VERIFY FINAL SCORE
        # ==========================================================

        echo ""
        echo "============================================================"
        echo "FINAL SCORE"
        echo "============================================================"

        if [[ ! -f "chr~{chromosome}.sscore" ]]; then
            echo "ERROR: .sscore file was not created."
            exit 1
        fi


        echo "Created:"
        ls -lh "chr~{chromosome}.sscore"

        echo ""
        echo "Score file line count:"
        wc -l "chr~{chromosome}.sscore"

        echo ""
        echo "First 10 lines:"
        head -10 "chr~{chromosome}.sscore"

        echo ""
        echo "============================================================"
        echo "CHROMOSOME ~{chromosome} COMPLETE"
        echo "============================================================"
        echo "End: $(date -u)"

    >>>


    # ==============================================================
    # OUTPUTS
    # ==============================================================

    output {

        File sscore =
            "chr" + chromosome + ".sscore"

        File score_file =
            "chr" + chromosome + "_aou_prs_score.tsv"

        File extract_file =
            "chr" + chromosome + "_prs_extract.txt"

        File match_table =
            "chr" + chromosome + "_prs_aou_matches.tsv"

        File match_summary =
            "chr" + chromosome + "_prs_match_summary.tsv"

        File log =
            "chr" + chromosome + ".log"

        File task_log =
            "chr" + chromosome + "_task.log"
    }


    # ==============================================================
    # RUNTIME
    # ==============================================================

    runtime {

        docker:
            "ubuntu:22.04"

        cpu:
            cpu

        memory:
            mem + " GB"

        disks:
            "local-disk 750 SSD"
    }
}