version 1.0


workflow AoU_chronotype_PRS_chr {

    input {

        File genetic_mt

        File prs_file

        String chromosome

        String output_prefix

        Int cpu = 16

        Int mem = 120

        Int threads = 16
    }


    call RunChromosomePRS {

        input:

            genetic_mt = genetic_mt,

            prs_file = prs_file,

            chromosome = chromosome,

            output_prefix = output_prefix,

            cpu = cpu,

            mem = mem,

            threads = threads
    }


    output {

        File bed = RunChromosomePRS.bed

        File bim = RunChromosomePRS.bim

        File fam = RunChromosomePRS.fam
    }

}



task RunChromosomePRS {


    input {

        File genetic_mt

        File prs_file

        String chromosome

        String output_prefix

        Int cpu

        Int mem

        Int threads
    }



    command <<<

        set -euo pipefail


        echo "=============================="
        echo "AoU Chromosome PRS"
        echo "Chromosome: ~{chromosome}"
        echo "=============================="


        echo "Files:"
        ls -lh


        python3 <<PYTHON


import hail as hl


print("Starting Hail")


hl.init(
    app_name="AoU_PRs_chr~{chromosome}"
)


print(
    "Hail version:",
    hl.version()
)



print("Loading WGS MT")


mt = hl.read_matrix_table(
    "~{genetic_mt}"
)


print(
    "Total variants:",
    mt.count_rows()
)



print("Loading PRS")


prs = hl.import_table(
    "~{prs_file}",
    delimiter=",",
    force_bgz=True,
    impute=True
)



print(
    "PRS variants:",
    prs.count()
)



print("Formatting PRS variants")


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



prs = prs.key_by(

    "locus",

    "alleles"

)



print("Filtering chromosome")


mt_chr = mt.filter_rows(

    mt.locus.contig ==

    "chr~{chromosome}"

)



print(
    "Chromosome variants:",
    mt_chr.count_rows()
)



print("Matching PRS variants")



mt_prs = mt_chr.filter_rows(

    hl.is_defined(

        prs[mt_chr.row_key]

    )

)



print(

    "Matched variants:",

    mt_prs.count_rows()

)



print("Annotating weights")



mt_prs = mt_prs.annotate_rows(

    prs_weight =

        prs[mt_prs.row_key].weight

)



print("Exporting PLINK")



hl.export_plink(

    mt_prs,

    output="~{output_prefix}"

)



print("DONE")



PYTHON


    >>>



    output {

        File bed = output_prefix + ".bed"

        File bim = output_prefix + ".bim"

        File fam = output_prefix + ".fam"

    }



    runtime {

        docker: "hailgenetics/hail:0.2.135"

        cpu: "~{cpu}"

        memory: "~{mem} GB"

        disks: "local-disk 500 SSD"

    }

}