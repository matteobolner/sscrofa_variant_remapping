rule all:
    input:
        #"data/complete_coords/complete_coords_v10.csv" #output of the first part of the pipeline
        "data/var_coords/var_coords_v11_report_duplicates.txt",
        "data/final_files/vars_sscrofa11.vcf" #output of the second part of the snakemake pipeline

rule get_coords:
    input:
        "data/starting_files/starting_data.csv"
    output:
        vars_file = "data/var_coords/var_coords_v10.csv",
        complete_coords = "data/complete_coords/complete_coords_v10.csv",
        genes_info = "data/genes_info/genes_info_v10.csv"
    script:
        "scripts/get_coords.py"

#must fix problem with perl in snakemake env first
rule remap_vars:
    input:
        "data/var_coords/var_coords_v10.csv"
    output:
        annot = "data/var_coords/var_coords_v11.csv",
        report = "data/var_coords/var_coords_v11.report"
    shell:
        "/usr/bin/perl5.30.3 /home/pelmo/data_and_pipelines/sscrofa_variant_remapping/scripts/remap_api.pl --mode asm-asm --annotation {input} --from GCF_000003025.5 --dest GCF_000003025.6 --annot_out {output.annot} --report_out {output.report}"

#the report file has more lines than variant positions due to multiple remappings; after finding the duplicates, the file was inspected and all lines with Contig630 (14) were removed to avoid problems
#UPDATE:appartently now it works

rule report_duplicates:
    input:
        "data/var_coords/var_coords_v11.report"
    output:
        "data/var_coords/var_coords_v11_report_duplicates.txt"
    shell:
        "cut -f 1 {input} | uniq -D > {output}"

rule remove_contig630:
    input:
        "data/var_coords/var_coords_v11.report"
    output:
        "data/var_coords/var_coords_v11_nocontig630.report"
    shell:
        "grep -v 'Contig630' {input} > {output}"


rule merge_results:
    input:
        "data/var_coords/var_coords_v11_nocontig630.report",
        "data/complete_coords/complete_coords_v10.csv"
    output:
        "data/complete_coords/merged_coords.csv",
        "data/final_files/unmapped_vars.csv"

    script:
        "scripts/merge_results.py"

#after manual curation of ids
rule update_ensembl_ids:
    input:
        "data/starting_files/remapped_ids.csv",
        "data/complete_coords/merged_coords.csv"
    output:
        "data/genes_info/genes_info_v11.csv",
        "data/complete_coords/merged_updated_coords_v11.csv"
    script:
        "scripts/update_ids.py"


rule verify_vars:
    input:
        "data/complete_coords/merged_updated_coords_v11.csv",
        "data/genes_info/genes_info_v10.csv",
        "data/genes_info/genes_info_v11.csv",
    output:
        "data/final_files/remapped_and_verified_vars.csv",
        "data/final_files/unsure_remappings.csv"
    script:
        "scripts/verify_vars.py"

rule update_vcf:
    input:
        "data/vcf_files/whole_POP_all.vcf",
        "data/final_files/remapped_and_verified_vars.csv"
    output:
        temp("data/vcf_files/remapped_vcf_no_header.vcf")
    script:
        "scripts/update_vcf.py"

rule update_vcf_header: #needs a manually curated header
    input:
        vcf = "data/vcf_files/remapped_vcf_no_header.vcf",
        header = "data/vcf_files/headers/updated_header.txt"
    output:
        "data/final_files/vars_sscrofa11.vcf"
    shell:
        "cat {input.header} {input.vcf} > {output}"
