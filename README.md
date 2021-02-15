# sscrofa_variant_remapping

Snakemake pipeline to remap a series of variant positions from the Sscrofa10.2 genome assembly to Sscrofa11.1.  
The method used for the remapping is the NCBI remap tool (https://www.ncbi.nlm.nih.gov/genome/tools/remap), which has a PERL API available (https://www.ncbi.nlm.nih.gov/genome/tools/remap/docs/api).  
Since the variant positions used in input are taken from Ensembl, there are some issues due to inconsistencies between Ensembl and NCBI; however, most variants can be correctly remapped.  

Format for the input variants file :  
gene_name|ensembl_id|chr|gene_start|gene_end|var_pos

Format for the input id update file:  
id_10_2,id_11_1,updated_id



TO DO:
-pipeline in docker container for better reproducibility
-push pipeline on workflowhub for better accessibility
-restructure the pipeline to take as input the full .VCF file and not the single #CHROM and POS columns as was done in this project to make the pipeline easier to implement on other VCF files
