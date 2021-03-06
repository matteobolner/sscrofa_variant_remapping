# Sus scrofa variant remapping

Snakemake pipeline to remap a series of variant positions from the Sscrofa10.2 genome assembly to Sscrofa11.1.
The method used for the remapping is the NCBI remap tool (https://www.ncbi.nlm.nih.gov/genome/tools/remap), which has a PERL API available (https://www.ncbi.nlm.nih.gov/genome/tools/remap/docs/api).  

Since the variant positions used in input are taken from Ensembl, there are some issues due to inconsistencies between Ensembl and NCBI; however, most variants can be correctly remapped.  

Format for the input variants file :  
gene_name|ensembl_id|chr|gene_start|gene_end|var_pos

Format for the input id update file:  
id_10_2,id_11_1,updated_id

All the generic libraries required by the pipeline are specified in the environment.yml file in workflow/environments and should be automatically fetched by snakemake.   

Since the PERL API requires specific versions and libraries to run, a Docker container encapsulating it was built starting from the PERL Dockerfile (see https://github.com/matteobolner/ncbi_remap_api_docker). The container runs a PERL environment from which the script can be succesfully executed.

**WARNING:**  
The resources and results folders present on this repository contain only the variants for one gene (DEFB114, 365 variants) to provide a showcase of the pipeline; in the thesis project the pipeline was applied to the whole 102948 variants.  

**REQUIREMENTS TO RUN THE PIPELINE:**  
- Snakemake  
- Docker and singularity  


**Graph of the pipeline:**  
![alt text](https://raw.githubusercontent.com/matteobolner/sscrofa_variant_remapping/247bd75af598d639d22bafbe92613706c290a65b/workflow/report/dag.svg)

**Workflow report**  
An HTML report of the final results is available in workflow/reports  




**TO DO:**  
-	~~remap api in docker container for better reproducibility~~  
- push pipeline on workflowhub for better accessibility  
- restructure the pipeline to take as input the full .VCF file and not the single #CHROM and POS columns as was done in this project to make the pipeline easier to implement on other VCF files  
