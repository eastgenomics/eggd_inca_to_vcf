# eggd_inca_to_vcf

## What does this app do?

Generates a VCF and aggregates variant information from the Inca database export CSV.

## What data are required for this app to run?

**Packages**

* Python packages as seen in resources/home/dnanexus/packages.

**Inputs**

* `input_file`: Export CSV file from the Inca database
* `output_filename`: Output filename for annotated VCF
* `genome_build`: Genome build to populate the VCF header with contigs. Options are: GRCh37, GRCh38
* `probeset`: Probeset or allele origin to filter the database, comma-separated if more than one option used. Options are: germline, somatic, 99347387, 96527893

## How to run

```bash
# in DNAnexus
dx run ${app_id} \
-iinput_file= \
-igenome_build= \
-iprobeset= \
-ioutput_filename= \
-y

# locally
python3 generate_inca_vcf.py \
--input_file ${input_file} \
--output_filename ${output_filename} \
--genome_build ${genome_build} \
--probeset ${probeset}
```

## What does this app output?

This app outputs an annotated VCF.
