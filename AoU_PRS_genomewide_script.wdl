version 1.0


workflow AoU_chronotype_PRS_genomewide {

    input {

        String genetic_mt

        File prs_file

        String output_prefix

        Int cpu = 16

        Int mem = 160

    }


    call CalculatePRS {

        input:

            genetic_mt = genetic_mt,

            prs_file = prs_file,

            output_prefix = output_prefix,

            cpu = cpu,

            mem = mem

    }


    output {

        File prs_output = CalculatePRS.prs_output

    }

}



task CalculatePRS {


    input {

        String genetic_mt

        File prs_file

        String output_prefix

        Int cpu

        Int mem

    }


    command <<<'PYTHON'

set -euo pipefail


python3 <<PYTHON


import hail as hl



PROJECT = "wb-happy-almond-4027"



############################################################
# Initialize Hail
############################################################


hl.init(

    app_name="AoU_chronotype_genomewide_PRS",

    backend="spark",

    tmp_dir=f"gs://dataproc-temp-{PROJECT}/hail_tmp",

    gcs_requester_pays_configuration=PROJECT

)



print("Hail version:", hl.version())



############################################################
# Load AoU WGS MatrixTable
############################################################


print("Loading AoU MatrixTable")



mt = hl.read_matrix_table(

    "~{genetic_mt}"

)



print("Total variants:")

print(mt.count_rows())



print("Total samples:")

print(mt.count_cols())



############################################################
# Variant QC for allele frequency fallback
############################################################


print("Running variant QC")



mt = hl.variant_qc(mt)



############################################################
# Import PRS weights
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


print("Creating PRS locus keys")



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

    "weight",

    "effect_allele"

)



prs = prs.key_by(

    "locus",

    "alleles"

)



print("PRS keyed variants:")

print(prs.count())



############################################################
# Keep only PRS variants
############################################################


print("Matching AoU variants")



mt = mt.filter_rows(

    hl.is_defined(

        prs[mt.row_key]

    )

)



matched = mt.count_rows()


print("Matched variants:")

print(matched)



############################################################
# Add PRS annotations
############################################################


mt = mt.annotate_rows(

    prs_weight =

        prs[mt.row_key].weight,


    prs_effect_allele =

        prs[mt.row_key].effect_allele

)



############################################################
# Determine allele orientation
############################################################


print("Determining allele orientation")



flip = hl.case() \

    .when(

        mt.prs_effect_allele == mt.alleles[0],

        True

    ) \

    .when(

        mt.prs_effect_allele == mt.alleles[1],

        False

    ) \

    .or_missing()



mt = mt.annotate_rows(

    flip = flip

)



############################################################
# Count flip direction
############################################################


print("REF effect allele variants:")

print(

    mt.aggregate_rows(

        hl.agg.count_where(mt.flip == True)

    )

)



print("ALT effect allele variants:")

print(

    mt.aggregate_rows(

        hl.agg.count_where(mt.flip == False)

    )

)



############################################################
# Calculate dosage
############################################################


dosage = hl.coalesce(

    hl.if_else(

        mt.flip,

        2 - mt.GT.n_alt_alleles(),

        mt.GT.n_alt_alleles()

    ),


    2 * hl.if_else(

        mt.flip,

        mt.variant_qc.AF[0],

        mt.variant_qc.AF[1]

    )

)



############################################################
# Calculate participant PRS
############################################################


print("Calculating PRS")



mt = mt.annotate_cols(

    chronotype_PRS =

        hl.agg.sum(

            mt.prs_weight * dosage

        )

)



############################################################
# Export
############################################################


print("Exporting PRS")



mt.cols() \

    .select(

        "chronotype_PRS"

    ) \

    .export(

        "~{output_prefix}"

    )



print("Finished")



PYTHON


    PYTHON



    output {

        File prs_output = output_prefix

    }


    runtime {

        docker: "hailgenetics/hail:0.2.135"

        cpu: cpu

        memory: "~{mem} GB"

        disks: "local-disk 1500 SSD"

    }

}