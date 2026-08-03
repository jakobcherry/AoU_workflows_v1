version 1.0


workflow AoU_chronotype_PRS_genomewide {


    input {

        String genetic_mt

        File prs_file

        String output_prefix

        Int cpu = 16

        Int memory_gb = 160

    }


    call CalculatePRS {


        input:

            genetic_mt = genetic_mt,

            prs_file = prs_file,

            output_prefix = output_prefix,

            cpu = cpu,

            memory_gb = memory_gb

    }


    output {

        Array[File] prs_output = CalculatePRS.prs_output

    }

}





task CalculatePRS {


    input {


        String genetic_mt

        File prs_file

        String output_prefix

        Int cpu

        Int memory_gb

    }



    command <<<


set -euo pipefail


python3 <<PYTHON


import os

import hail as hl



PROJECT = "wb-happy-almond-4027"



############################################################
# Spark configuration
############################################################


os.environ["PYSPARK_SUBMIT_ARGS"] = (

    "--conf spark.driver.memory=140g "

    "--conf spark.executor.memory=140g "

    "pyspark-shell"

)



############################################################
# Initialize Hail
############################################################


print("Starting Hail")


hl.init(

    app_name="AoU_chronotype_PRS_genomewide",

    backend="spark",

    tmp_dir="gs://dataproc-temp-wb-happy-almond-4027/hail_tmp",

    gcs_requester_pays_configuration=PROJECT

)


print("Hail version:", hl.version())





############################################################
# Load AoU MatrixTable
############################################################


print("Loading AoU WGS MatrixTable")


mt = hl.read_matrix_table(

    "~{genetic_mt}"

)


print("Total variants:")

print(mt.count_rows())


print("Total samples:")

print(mt.count_cols())





############################################################
# Load PRS weights
############################################################


print("Loading PRS weights")


prs = hl.import_table(

    "~{prs_file}",

    delimiter=",",

    force_bgz=True,

    impute=True

)


print("PRS variants:")

print(prs.count())





############################################################
# Create PRS keys
############################################################


print("Creating PRS keys")


prs = prs.annotate(

    locus = hl.locus(

        hl.format(

            "chr%s",

            prs.chr

        ),

        prs.bp,

        reference_genome="GRCh38"

    ),

    alleles = [

        prs.noneffect_allele,

        prs.effect_allele

    ]

)



prs = prs.select(

    "locus",

    "alleles",

    "effect_allele",

    "weight"

)



prs = prs.key_by(

    "locus",

    "alleles"

)





############################################################
# Match AoU variants
############################################################


print("Filtering AoU variants to PRS variants")


mt = mt.filter_rows(

    hl.is_defined(

        prs[mt.row_key]

    )

)


print("Matched variants:")

print(mt.count_rows())





############################################################
# Variant QC
############################################################


print("Running variant QC")


mt = hl.variant_qc(mt)





############################################################
# Add PRS information
############################################################


print("Adding PRS annotations")


mt = mt.annotate_rows(

    prs_weight = prs[mt.row_key].weight,

    prs_effect_allele = prs[mt.row_key].effect_allele

)





############################################################
# Determine allele orientation
############################################################


print("Determining allele orientation")


flip = (

    hl.case()

    .when(

        mt.prs_effect_allele == mt.alleles[0],

        True

    )

    .when(

        mt.prs_effect_allele == mt.alleles[1],

        False

    )

    .or_missing()

)



mt = mt.annotate_rows(

    flip = flip

)





print("REF matches:")

print(

    mt.aggregate_rows(

        hl.agg.count_where(mt.flip == True)

    )

)



print("ALT matches:")

print(

    mt.aggregate_rows(

        hl.agg.count_where(mt.flip == False)

    )

)



print("Unresolved:")

print(

    mt.aggregate_rows(

        hl.agg.count_where(hl.is_missing(mt.flip))

    )

)





############################################################
# Remove unresolved variants
############################################################


print("Removing unresolved variants")


mt = mt.filter_rows(

    hl.is_defined(mt.flip)

)





############################################################
# Calculate effect allele dosage
############################################################


print("Calculating dosage")


dosage = hl.if_else(

    mt.flip,

    2 - mt.GT.n_alt_alleles(),

    mt.GT.n_alt_alleles()

)





############################################################
# Calculate PRS
############################################################


print("Calculating genome-wide PRS")


mt = mt.annotate_cols(

    chronotype_PRS = hl.agg.sum(

        mt.prs_weight * dosage

    )

)





############################################################
# Export
############################################################


print("Exporting PRS")


prs_results = mt.cols().select(

    "chronotype_PRS"

)


prs_results.export(

    "~{output_prefix}.bgz"

)



print("Finished successfully")


PYTHON


    >>>


    output {

        Array[File] prs_output = glob("~{output_prefix}.bgz*")

    }



    runtime {


        docker: "hailgenetics/hail:0.2.135"


        cpu: cpu


        memory: "~{memory_gb} GB"


        disks: "local-disk 2000 SSD"

    }


}