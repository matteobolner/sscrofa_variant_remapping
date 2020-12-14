# sscrofa_variant_remapping

Snakemake pipeline to remap a series of variant positions from the Sscrofa10.2 genome assembly to Sscrofa11.1.  
The method used for the remapping is the NCBI remap tool (https://www.ncbi.nlm.nih.gov/genome/tools/remap), which has a PERL API available (https://www.ncbi.nlm.nih.gov/genome/tools/remap/docs/api).  
Since the variant positions used in input are taken from Ensembl, there are some issues with due to inconsistencies between Ensembl and NCBI; however, most variants can be correctly remapped.  
