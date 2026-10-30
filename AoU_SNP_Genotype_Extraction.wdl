
version 1.0

# ============================================================
# AoU v9 CHRONOTYPE VARIANT GENOTYPE EXTRACTION
#
# Input: hg38 TSV produced by the R liftover script.
#
# Required columns:
#   chr
#   bp
#   rs_number
#   allele1
#   allele2
#
# Optional columns:
#   source
#   original_chr
#   original_bp_hg37
#   original_allele1
#   original_allele2
#   liftover_strand
#   allele_orientation_status
#   input_row
#
# No weights, PRS calculation, or liftover.
#
# Final genotype dosage:
#   Dosage of the input allele1.
#
# Allele orientation:
#   Direct and reversed BIM allele order are both accepted.
#   Dosage is counted from the actual VCF allele sequence,
#   not from the BIM A1/A2 ordering.
#
# D/I:
#   Reported separately for later indel harmonization.
#   They do not stop extraction of other variants.
#
# chromosome_selector:
#   0  = chromosomes 1-22 and X
#   1-22 = specified autosome
#   23 = X
#
# Resources:
#   CPU: 16
#   RAM: 120 GB
#   Disk: 1875 GB SSD
# ============================================================


workflow AoU_SNP_Genotype_Extraction {

  input {
    File variant_file

    File plink2_file = "gs://cloned-all-circadian-patients-for-fitbit-and-prs-wb-perky-onion/plink2"

    String plink_bucket = "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/plink_bed/acaf_threshold"

    Int cpu = 16
    Int mem = 120

    # 0 = all chromosomes
    # 1-22 = autosome
    # 23 = X
    Int chromosome_selector = 0
  }

  Array[String] chromosome_labels = [
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
    "12", "13", "14", "15", "16", "17", "18", "19", "20", "21",
    "22", "X"
  ]

  Array[Int] selected_indices =
    if (chromosome_selector == 0) {
      range(23)
    } else {
      [chromosome_selector - 1]
    }

  scatter (idx in selected_indices) {

    File bedfile = "~{plink_bucket}.chr~{chromosome_labels[idx]}.bed"
    File bimfile = "~{plink_bucket}.chr~{chromosome_labels[idx]}.bim"
    File famfile = "~{plink_bucket}.chr~{chromosome_labels[idx]}.fam"

    call ExtractChromosomeGenotypes {
      input:
        variant_file = variant_file,
        plink2_file = plink2_file,
        bedfile = bedfile,
        bimfile = bimfile,
        famfile = famfile,
        chromosome_label = chromosome_labels[idx],
        cpu = cpu,
        mem = mem
    }
  }

  output {
    Array[File] genotype_matrices =
      ExtractChromosomeGenotypes.genotype_matrix

    Array[File] variant_matches =
      ExtractChromosomeGenotypes.variant_match

    Array[File] variant_qc =
      ExtractChromosomeGenotypes.variant_qc

    Array[File] indel_review =
      ExtractChromosomeGenotypes.indel_review

    Array[File] extraction_summaries =
      ExtractChromosomeGenotypes.extraction_summary

    Array[File] plink_logs =
      ExtractChromosomeGenotypes.plink_log
  }
}


task ExtractChromosomeGenotypes {

  input {
    File variant_file
    File plink2_file
    File bedfile
    File bimfile
    File famfile

    String chromosome_label

    Int cpu
    Int mem
  }

  command <<<
    set -euo pipefail

    # --------------------------------------------------------
    # 1. LOCALIZE AND CHECK PLINK FILES
    # --------------------------------------------------------

    PLINK2="./$(basename "~{plink2_file}")"

    chmod 750 "$PLINK2"

    BFILE_PREFIX="$(basename "~{bedfile}" .bed)"

    test -f "${BFILE_PREFIX}.bed"
    test -f "${BFILE_PREFIX}.bim"
    test -f "${BFILE_PREFIX}.fam"

    test -s "$PLINK2"

    # --------------------------------------------------------
    # 2. PYTHON VARIANT MATCHING AND DOSAGE EXTRACTION
    # --------------------------------------------------------

    cat > extract_genotypes.py <<'PYTHON'
import csv
import sys
import os
import re
from collections import defaultdict


CHROMOSOME = sys.argv[1]

MATCH_FIELDS = [
    "input_row",
    "chr",
    "bp",
    "rs_number",
    "allele1",
    "allele2",
    "source",
    "original_chr",
    "original_bp_hg37",
    "original_allele1",
    "original_allele2",
    "liftover_strand",
    "allele_orientation_status",
    "variant_column",
    "matched_variant_id",
    "bim_a1",
    "bim_a2",
    "allele1_vs_bim",
    "match_status",
    "match_detail"
]

QC_FIELDS = [
    "input_row",
    "chr",
    "bp",
    "rs_number",
    "variant_column",
    "allele1",
    "allele2",
    "bim_a1",
    "bim_a2",
    "allele1_vs_bim",
    "match_status",
    "n_samples",
    "n_called",
    "n_missing",
    "missing_rate",
    "dosage_0",
    "dosage_1",
    "dosage_2",
    "n_haploid_calls",
    "n_diploid_calls"
]

SUMMARY_FIELDS = [
    "chromosome",
    "n_input_variants",
    "n_matched",
    "n_not_found",
    "n_allele_mismatch",
    "n_ambiguous",
    "n_indel_review",
    "n_invalid",
    "n_genotype_variants",
    "n_samples",
    "n_genotype_rows"
]


def norm_chr(value):
    value = str(value).strip().upper()

    if value.startswith("CHR"):
        value = value[3:]

    if value == "23":
        value = "X"

    if value == "24":
        value = "Y"

    if value == "M":
        value = "MT"

    return value


def norm_allele(value):
    return str(value).strip().upper()


def write_tsv(path, fields, rows):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n"
        )

        writer.writeheader()
        writer.writerows(rows)


def read_fam(path):
    samples = []

    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.strip().split()

            if len(fields) < 2:
                continue

            fid = fields[0]
            iid = fields[1]

            # AoU person_id is expected to be the IID.
            samples.append(iid)

    if len(samples) != len(set(samples)):
        raise RuntimeError(
            "Duplicate IID values in FAM file. "
            "Cannot safely construct a unique person_id matrix."
        )

    return samples


def read_variants(path):
    variants = []

    with open(path, "r", newline="", encoding="utf-8-sig") as handle:

        reader = csv.DictReader(handle, delimiter="\t")

        if not reader.fieldnames:
            raise RuntimeError("Variant TSV has no header.")

        columns = {
            name.strip().lower(): name
            for name in reader.fieldnames
            if name is not None
        }

        required = [
            "chr",
            "bp",
            "rs_number",
            "allele1",
            "allele2"
        ]

        missing = [x for x in required if x not in columns]

        if missing:
            raise RuntimeError(
                "Missing required TSV columns: "
                + ", ".join(missing)
            )

        optional = [
            "source",
            "original_chr",
            "original_bp_hg37",
            "original_allele1",
            "original_allele2",
            "liftover_strand",
            "allele_orientation_status",
            "input_row"
        ]

        for row_number, row in enumerate(reader, start=2):

            chromosome = norm_chr(row.get(columns["chr"], ""))

            if chromosome != norm_chr(CHROMOSOME):
                continue

            item = {
                "input_row": row.get(
                    columns.get("input_row", ""),
                    str(row_number)
                ) or str(row_number),

                "chr": chromosome,

                "bp": str(
                    row.get(columns["bp"], "")
                ).strip(),

                "rs_number": str(
                    row.get(columns["rs_number"], "")
                ).strip(),

                "allele1": norm_allele(
                    row.get(columns["allele1"], "")
                ),

                "allele2": norm_allele(
                    row.get(columns["allele2"], "")
                )
            }

            for field in optional:
                if field in columns:
                    item[field] = str(
                        row.get(columns[field], "")
                    ).strip()
                else:
                    item[field] = ""

            # A unique output column name for every input row.
            safe_rsid = re.sub(
                r"[^A-Za-z0-9_.-]",
                "_",
                item["rs_number"] or "variant"
            )

            item["variant_column"] = (
                safe_rsid
                + "_chr"
                + chromosome
                + "_bp"
                + item["bp"]
                + "_row"
                + str(row_number)
            )

            variants.append(item)

    return variants


def read_bim(path):
    records = []

    with open(path, "r", encoding="utf-8") as handle:

        for line in handle:

            fields = line.strip().split()

            if len(fields) < 6:
                continue

            try:
                bp = int(fields[3])
            except ValueError:
                continue

            records.append({
                "chr": norm_chr(fields[0]),
                "id": fields[1],
                "bp": bp,
                "a1": norm_allele(fields[4]),
                "a2": norm_allele(fields[5])
            })

    return records


def is_symbolic_indel(a1, a2):
    return a1 in ("D", "I") or a2 in ("D", "I")


def prepare():

    variants = read_variants(sys.argv[2])
    bim = read_bim(sys.argv[3])

    by_position = defaultdict(list)

    for record in bim:
        by_position[
            (record["chr"], record["bp"])
        ].append(record)

    match_rows = []
    matched_ids = []

    for variant in variants:

        row = {
            **variant,
            "matched_variant_id": "",
            "bim_a1": "",
            "bim_a2": "",
            "allele1_vs_bim": "",
            "match_status": "",
            "match_detail": ""
        }

        try:
            bp = int(variant["bp"])

            if bp <= 0:
                raise ValueError()

        except ValueError:
            row["match_status"] = "INVALID_POSITION"
            row["match_detail"] = "bp is not a positive integer"
            match_rows.append(row)
            continue

        a1 = variant["allele1"]
        a2 = variant["allele2"]

        if not a1 or not a2 or a1 == a2:
            row["match_status"] = "INVALID_ALLELES"
            row["match_detail"] = (
                "Alleles are missing or identical"
            )
            match_rows.append(row)
            continue

        candidates = by_position.get(
            (variant["chr"], bp),
            []
        )

        if not candidates:
            row["match_status"] = "NOT_FOUND"
            row["match_detail"] = (
                "No BIM record at this chromosome and position"
            )
            match_rows.append(row)
            continue

        # D/I labels are not nucleotide sequences.
        # Keep them in the report without attempting an unsafe match.
        if is_symbolic_indel(a1, a2):

            row["match_status"] = "INDEL_REVIEW_POSITION_FOUND"

            row["match_detail"] = (
                "D/I allele label found. BIM position exists, "
                "but the sequence representation must be "
                "harmonized before genotype extraction. "
                "BIM candidates: "
                + ";".join(
                    r["id"] + ":" + r["a1"] + "/" + r["a2"]
                    for r in candidates
                )
            )

            if len(candidates) == 1:
                row["matched_variant_id"] = candidates[0]["id"]
                row["bim_a1"] = candidates[0]["a1"]
                row["bim_a2"] = candidates[0]["a2"]

            match_rows.append(row)
            continue

        wanted = frozenset([a1, a2])

        exact = [
            record for record in candidates
            if frozenset([record["a1"], record["a2"]]) == wanted
        ]

        if not exact:

            observed = sorted(set(
                record["a1"] + "/" + record["a2"]
                for record in candidates
            ))

            row["match_status"] = "ALLELE_MISMATCH"

            row["match_detail"] = (
                "Position found, but allele pair did not match. "
                "BIM allele pairs: "
                + ",".join(observed)
            )

            match_rows.append(row)
            continue

        # If multiple identical BIM records exist, use rsID
        # only when it uniquely resolves the match.
        if len(exact) > 1 and variant["rs_number"]:

            rs_matches = [
                record for record in exact
                if record["id"] == variant["rs_number"]
            ]

            if len(rs_matches) == 1:
                exact = rs_matches

        if len(exact) != 1:

            row["match_status"] = "AMBIGUOUS"

            row["match_detail"] = (
                "Multiple BIM records match position and alleles"
            )

            match_rows.append(row)
            continue

        record = exact[0]

        if not record["id"] or record["id"] == ".":
            row["match_status"] = "AMBIGUOUS"
            row["match_detail"] = (
                "Matched BIM record has no usable variant ID"
            )
            match_rows.append(row)
            continue

        row["matched_variant_id"] = record["id"]
        row["bim_a1"] = record["a1"]
        row["bim_a2"] = record["a2"]

        if a1 == record["a1"]:
            row["allele1_vs_bim"] = "BIM_A1"
            row["match_detail"] = (
                "Direct allele order: input allele1 equals BIM A1"
            )

        elif a1 == record["a2"]:
            row["allele1_vs_bim"] = "BIM_A2"
            row["match_detail"] = (
                "Reversed allele order: input allele1 equals BIM A2"
            )

        else:
            row["match_status"] = "ALLELE_MISMATCH"
            row["match_detail"] = "Unexpected allele orientation"
            match_rows.append(row)
            continue

        row["match_status"] = "MATCHED"

        matched_ids.append(record["id"])
        match_rows.append(row)

    write_tsv(
        "variant_match.tsv",
        MATCH_FIELDS,
        match_rows
    )

    # PLINK --extract uses variant IDs.
    # Duplicate IDs are harmless for extraction; the VCF is
    # subsequently checked by chromosome, position and alleles.
    with open("extract_ids.txt", "w", encoding="utf-8") as handle:

        for variant_id in sorted(set(matched_ids)):
            handle.write(variant_id + "\n")

    # Separate indel report.
    indel_rows = [
        row for row in match_rows
        if row["match_status"].startswith("INDEL_REVIEW")
    ]

    write_tsv(
        "indel_review.tsv",
        MATCH_FIELDS,
        indel_rows
    )


def parse_vcf():

    with open("variant_match.tsv", "r", newline="", encoding="utf-8") as handle:
        match_rows = list(
            csv.DictReader(handle, delimiter="\t")
        )

    fam_samples = read_fam(sys.argv[4])

    matched = [
        row for row in match_rows
        if row["match_status"] == "MATCHED"
    ]

    # Key by chromosome, position and unordered allele pair.
    # This avoids relying on variant ID alone.
    by_key = defaultdict(list)

    for row in matched:

        key = (
            norm_chr(row["chr"]),
            int(row["bp"]),
            frozenset([
                row["allele1"],
                row["allele2"]
            ])
        )

        by_key[key].append(row)

    # Initialize an output matrix with every FAM sample.
    matrix = {
        sample: {}
        for sample in fam_samples
    }

    qc = {}

    for row in matched:

        qc[row["input_row"]] = {
            "input_row": row["input_row"],
            "chr": row["chr"],
            "bp": row["bp"],
            "rs_number": row["rs_number"],
            "variant_column": row["variant_column"],
            "allele1": row["allele1"],
            "allele2": row["allele2"],
            "bim_a1": row["bim_a1"],
            "bim_a2": row["bim_a2"],
            "allele1_vs_bim": row["allele1_vs_bim"],
            "match_status": "MATCHED_NOT_EXPORTED",
            "n_samples": 0,
            "n_called": 0,
            "n_missing": 0,
            "missing_rate": "",
            "dosage_0": 0,
            "dosage_1": 0,
            "dosage_2": 0,
            "n_haploid_calls": 0,
            "n_diploid_calls": 0
        }

        for sample in fam_samples:
            matrix[sample][row["variant_column"]] = "NA"

    vcf_samples = []

    with open("plink_extract.vcf", "r", encoding="utf-8") as handle:

        for line in handle:

            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):

                header = line.rstrip("\n").split("\t")
                vcf_samples = header[9:]

                if vcf_samples != fam_samples:
                    raise RuntimeError(
                        "VCF sample IDs/order do not match FAM IID order. "
                        "Cannot safely assign genotypes to person_id. "
                        "Check PLINK VCF sample-ID behavior."
                    )

                continue

            if line.startswith("#") or not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 9:
                continue

            chromosome = norm_chr(fields[0])

            try:
                bp = int(fields[1])
            except ValueError:
                continue

            ref = norm_allele(fields[3])
            alt = norm_allele(fields[4])

            # Only biallelic records are accepted here.
            if "," in alt:
                continue

            key = (
                chromosome,
                bp,
                frozenset([ref, alt])
            )

            candidates = by_key.get(key, [])

            if len(candidates) != 1:
                continue

            row = candidates[0]

            # The input allele1 must be exactly REF or ALT.
            if row["allele1"] == ref:
                target_vcf_allele = "0"
            elif row["allele1"] == alt:
                target_vcf_allele = "1"
            else:
                qc[row["input_row"]]["match_status"] = (
                    "VCF_ALLELE_MISMATCH"
                )
                continue

            format_fields = fields[8].split(":")

            if "GT" not in format_fields:
                qc[row["input_row"]]["match_status"] = (
                    "VCF_GT_FIELD_MISSING"
                )
                continue

            gt_index = format_fields.index("GT")

            q = qc[row["input_row"]]
            q["match_status"] = "EXTRACTED"

            for sample, sample_field in zip(
                vcf_samples,
                fields[9:]
            ):

                sample_parts = sample_field.split(":")

                gt = (
                    sample_parts[gt_index]
                    if gt_index < len(sample_parts)
                    else "."
                )

                alleles = gt.replace("|", "/").split("/")
                ploidy = len(alleles)

                dosage = None

                # Count copies of input allele1.
                # Haploid calls naturally produce dosage 0 or 1.
                # Diploid calls produce dosage 0, 1 or 2.
                if (
                    ploidy > 0
                    and all(a in ("0", "1") for a in alleles)
                ):

                    dosage = sum(
                        1
                        for allele in alleles
                        if allele == target_vcf_allele
                    )

                if dosage is None:
                    matrix[sample][row["variant_column"]] = "NA"
                    q["n_missing"] += 1

                else:
                    matrix[sample][row["variant_column"]] = dosage

                    q["n_called"] += 1
                    q["dosage_" + str(dosage)] += 1

                    if ploidy == 1:
                        q["n_haploid_calls"] += 1
                    elif ploidy == 2:
                        q["n_diploid_calls"] += 1

                q["n_samples"] += 1

    # Calculate missingness.
    for q in qc.values():

        if q["n_samples"] > 0:
            q["missing_rate"] = (
                q["n_missing"] / q["n_samples"]
            )

    # --------------------------------------------------------
    # WRITE WIDE PERSON-LEVEL MATRIX
    # --------------------------------------------------------

    variant_columns = [
        row["variant_column"]
        for row in matched
    ]

    with open(
        "genotype_matrix.tsv",
        "w",
        newline="",
        encoding="utf-8"
    ) as handle:

        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n"
        )

        writer.writerow(["person_id"] + variant_columns)

        for sample in fam_samples:

            writer.writerow(
                [sample]
                + [
                    matrix[sample].get(column, "NA")
                    for column in variant_columns
                ]
            )

    write_tsv(
        "variant_qc.tsv",
        QC_FIELDS,
        list(qc.values())
    )

    write_summary(
        match_rows,
        len(fam_samples),
        len(variant_columns),
        sum(
            1 for q in qc.values()
            if q["match_status"] == "EXTRACTED"
        )
    )


def empty_outputs():

    with open("variant_match.tsv", "r", newline="", encoding="utf-8") as handle:
        match_rows = list(
            csv.DictReader(handle, delimiter="\t")
        )

    fam_samples = read_fam(sys.argv[4])

    with open(
        "genotype_matrix.tsv",
        "w",
        newline="",
        encoding="utf-8"
    ) as handle:

        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n"
        )

        writer.writerow(["person_id"])

        for sample in fam_samples:
            writer.writerow([sample])

    qc_rows = []

    for row in match_rows:

        if row["match_status"] == "MATCHED":
            status = "MATCHED_NO_GENOTYPES_EXPORTED"
        else:
            status = row["match_status"]

        qc_rows.append({
            "input_row": row["input_row"],
            "chr": row["chr"],
            "bp": row["bp"],
            "rs_number": row["rs_number"],
            "variant_column": row["variant_column"],
            "allele1": row["allele1"],
            "allele2": row["allele2"],
            "bim_a1": row["bim_a1"],
            "bim_a2": row["bim_a2"],
            "allele1_vs_bim": row["allele1_vs_bim"],
            "match_status": status,
            "n_samples": 0,
            "n_called": 0,
            "n_missing": 0,
            "missing_rate": "",
            "dosage_0": 0,
            "dosage_1": 0,
            "dosage_2": 0,
            "n_haploid_calls": 0,
            "n_diploid_calls": 0
        })

    write_tsv("variant_qc.tsv", QC_FIELDS, qc_rows)

    write_summary(
        match_rows,
        len(fam_samples),
        0,
        0
    )


def write_summary(match_rows, n_samples, n_genotype_variants, n_extracted):

    counts = defaultdict(int)

    for row in match_rows:
        counts[row["match_status"]] += 1

    n_genotype_rows = n_samples * n_genotype_variants

    summary = {
        "chromosome": CHROMOSOME,
        "n_input_variants": len(match_rows),
        "n_matched": counts["MATCHED"],
        "n_not_found": counts["NOT_FOUND"],
        "n_allele_mismatch": counts["ALLELE_MISMATCH"],
        "n_ambiguous": counts["AMBIGUOUS"],
        "n_indel_review": (
            counts["INDEL_REVIEW_POSITION_FOUND"]
            + counts["INDEL_REVIEW_POSITION_NOT_FOUND"]
        ),
        "n_invalid": (
            counts["INVALID_POSITION"]
            + counts["INVALID_ALLELES"]
        ),
        "n_genotype_variants": n_extracted,
        "n_samples": n_samples,
        "n_genotype_rows": n_genotype_rows
    }

    write_tsv(
        "extraction_summary.tsv",
        SUMMARY_FIELDS,
        [summary]
    )


if __name__ == "__main__":

    mode = sys.argv[2]

    if mode == "prepare":
        prepare()

    elif mode == "parse":
        parse_vcf()

    elif mode == "empty":
        empty_outputs()

    else:
        raise RuntimeError("Unknown mode: " + mode)

PYTHON

    # --------------------------------------------------------
    # 3. MATCH INPUT VARIANTS AGAINST BIM
    # --------------------------------------------------------

    python3 extract_genotypes.py \
      "~{chromosome_label}" \
      prepare \
      "~{variant_file}" \
      "${BFILE_PREFIX}.bim"

    # --------------------------------------------------------
    # 4. EXTRACT MATCHED VARIANTS
    # --------------------------------------------------------

    if [ -s extract_ids.txt ]; then

      "${PLINK2}" \
        --bfile "${BFILE_PREFIX}" \
        --extract extract_ids.txt \
        --export vcf id-paste=iid \
        --threads ~{cpu} \
        --memory ~{mem}000 \
        --out plink_extract

      test -s plink_extract.vcf

      python3 extract_genotypes.py \
        "~{chromosome_label}" \
        parse \
        plink_extract.vcf \
        "${BFILE_PREFIX}.fam"

    else

      python3 extract_genotypes.py \
        "~{chromosome_label}" \
        empty \
        unused \
        "${BFILE_PREFIX}.fam"

      echo "No extractable SNP variants on chromosome ~{chromosome_label}." \
        > plink_extract.log

    fi

    # Ensure an indel report exists even when no indels occur.
    test -f indel_review.tsv

    # Ensure a log exists if PLINK did not create one.
    if [ ! -f plink_extract.log ]; then
      touch plink_extract.log
    fi

  >>>

  output {
    File genotype_matrix = "genotype_matrix.tsv"

    File variant_match = "variant_match.tsv"

    File variant_qc = "variant_qc.tsv"

    File indel_review = "indel_review.tsv"

    File extraction_summary = "extraction_summary.tsv"

    File plink_log = "plink_extract.log"
  }

  runtime {
    docker: "python:3.11-slim"
    cpu: cpu
    memory: "~{mem} GB"
    disks: "local-disk 1875 SSD"
  }
}