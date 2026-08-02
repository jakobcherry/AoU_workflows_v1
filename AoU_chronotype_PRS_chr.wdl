version 1.0


workflow AoU_chronotype_PRS_chr {


    input {

        String chromosome

        String genetic_mt =
            "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/multiMT/hail.mt"

        File prs_file =
            "gs://full-cohort-for-prs-wb-happy-almond-4027/chronotype_AoUPRS_input.csv.bgz"

        String output_dir =
            "gs://full-cohort-for-prs-wb-happy-almond-4027/chronometa_chr_prs/"


    }



    call RunChromosomePRS {

        input:

            chromosome = chromosome,
            genetic_mt = genetic_mt,
            prs_file = prs_file,
            output_dir = output_dir

    }



    output {


        File prs_score =
            RunChromosomePRS.prs_score


        File bed =
            RunChromosomePRS.bed


        File bim =
            RunChromosomePRS.bim


        File fam =
            RunChromosomePRS.fam


    }

}




task RunChromosomePRS {


    input {


        String chromosome

        String genetic_mt

        File prs_file

        String output_dir

    }



    runtime {

        docker:
            "us.gcr.io/broad-gatk/hail:0.2.135"

        cpu:
            32

        memory:
            "120 GB"

        disks:
            "local-disk 750 SSD"

    }



    command <<<


set -euo pipefail


echo "Running chromosome ~{chromosome}"


mkdir -p work



python3 <<'PY'


import hail as hl


chromosome = "~{chromosome}"

genetic_mt = "~{genetic_mt}"

prs_file = "~{prs_file}"



hl.init(
    app_name=f"AoU_chr{chromosome}_PRS",
    tmp_dir="work/hail_tmp"
)



print("Loading WGS")


mt = hl.read_matrix_table(
    genetic_mt
)



print("Filtering chromosome")


mt = mt.filter_rows(
    mt.locus.contig ==
    f"chr{chromosome}"
)



print(
    "Chromosome variants:",
    mt.count_rows()
)



print("Importing PRS")



prs = hl.import_table(
    prs_file,
    delimiter=",",
    force_bgz=True,
    impute=True
)



prs = prs.filter(
    prs.chr ==
    int(chromosome)
)



print(
    "PRS variants:",
    prs.count()
)



prs = prs.annotate(

    locus =
        hl.locus(
            hl.format(
                "chr%s",
                prs.chr
            ),
            prs.bp,
            reference_genome="GRCh38"
        ),

    alleles =
        [
            prs.noneffect_allele,
            prs.effect_allele
        ]

)



prs = prs.key_by(
    "locus",
    "alleles"
)



print("Matching")



mt = mt.filter_rows(
    hl.is_defined(
        prs[mt.row_key]
    )
)



print(
    "Matched variants:",
    mt.count_rows()
)



mt = mt.annotate_rows(

    rsid =
        prs[mt.row_key].rs_number,

    weight =
        prs[mt.row_key].weight

)



print("Exporting PLINK")



hl.export_plink(
    mt,
    output=f"work/chr{chromosome}"
)



print("Creating scoring file")



prs.select(
    SNP = prs.rs_number,
    A1 = prs.effect_allele,
    BETA = prs.weight

).export(
    f"work/chr{chromosome}_score.txt"
)



print("Finished Hail")



PY



echo "Running PLINK2 score"



plink2 \
    --bfile work/chr~{chromosome} \
    --score work/chr~{chromosome}_score.txt \
    1 2 3 header-read \
    cols=+scoresums \
    --out work/chr~{chromosome}_PRS



echo "Copying outputs"



gsutil cp \
    work/chr~{chromosome}_PRS.sscore \
    ~{output_dir}/chronotype_chr~{chromosome}_PRS.sscore



gsutil cp \
    work/chr~{chromosome}.bed \
    ~{output_dir}/chronotype_chr~{chromosome}.bed



gsutil cp \
    work/chr~{chromosome}.bim \
    ~{output_dir}/chronotype_chr~{chromosome}.bim



gsutil cp \
    work/chr~{chromosome}.fam \
    ~{output_dir}/chronotype_chr~{chromosome}.fam



    >>>


    output {


        File prs_score =
            "~{output_dir}/chronotype_chr~{chromosome}_PRS.sscore"



        File bed =
            "~{output_dir}/chronotype_chr~{chromosome}.bed"



        File bim =
            "~{output_dir}/chronotype_chr~{chromosome}.bim"



        File fam =
            "~{output_dir}/chronotype_chr~{chromosome}.fam"



    }

}