This github repository houses a wrapper program (`microSLURM_16S.sh`) for processing 16S amplicon sequences (derived from Illumina paired-end sequencing) to detect amplicon sequence variant (ASV) profiles (ASV abundance counts per sample) on a high performance computing cluster using a SLURM scheduling system. The overall pipeline includes performing an intitial quality assessment on raw sequences using FastQC [https://www.bioinformatics.babraham.ac.uk/projects/fastqc/] and MultiQC [https://multiqc.info/], removal of PCR primers using cutadapt [https://cutadapt.readthedocs.io/en/stable/], quality trimming/filtering, ASV detection, and taxonomic profiling using DADA2 R package [https://benjjneb.github.io/dada2/], and phylogenetic tree construction of ASVs using DECIPHER [http://www2.decipher.codes/] and phangorn [https://github.com/KlausVigo/phangorn] R packages. Alternatively, RAxML-NG [https://github.com/amkozlov/raxml-ng] can be used for phylogenetic tree construction. All generated data is packaged nicely into a phyloseq object containing ASV abundances, taxonomic assignments, and phylogenetic tree using phyloseq R package [https://joey711.github.io/phyloseq/], that can be used in further statistical analyses.

The following gives an overview of the overall structure of the repository:

## Directory tree for repository
```
microSLURM_16S
|
|-- Environment -- Contains a .yml file that can be used to build a conda environment
|                  that contains all the necessary programs to run the pipeline.
|
|-- Reference_Files -- Directory that contains the taxonomic reference files 
|                      used in the pipeline.
|
|-- microSLURM_16S.job -- An example sbatch job script for submitting the 
|                         microSLURM_16S.sh to a SLURM scheduling 
|                         system on a high performance computing cluster.
|
|-- microSLURM_16S.sh -- The wrapper program that runs the 16S amplicon pipeline.

```
## Important notes about the pipeline program

### Submission of internal sbatch jobs
The `microSLURM_16S.sh` script will internally submit jobs for each step of the pipeline. For each job step, a `run.sh` file is created in the current directory that contains the code for the currently running step, and is deleted once the step completes. Some steps will also produce a `run.R` script that will also be deleted onece the step finishes. Partitions, time limits, number of cores (cpus), and memory per cpu can be requested through the wrapper scripts to meet the demands of a particular SLURM scheduler, dataset, etc.

### Required programs/databases and parameter descriptions
For descriptions of required programs/databases and parameters for `microSLURM_16S.sh` script, run the script with parameter `-h`. As stated in the directory tree above, the directory `Environment` contains a `.yml` file that can be used to build a conda environment with all the necessary programs for running the `microSLURM_16S.sh` script.

### Removing intermediate sequence files generated during pipeline
During running of the `microSLURM_16S.sh` script, sequence files are generated at the primer removal step and after quality trimming/filtering step. This may take up a lot of storage space depending on the dataset size, so if you would like intermediate files removed at the end of the pipeline, keeping only the processed, quality controlled sequences (that come after quality trimming/filtering step), use the flag `-d` to remove intermediate sequence files and the corresponding directory. A combined log file for each sample will be created joining all logs from both steps. As the numbering of directories no longer makes sense after removing the intermediate sequence directory, the pipeline output directory is reorganized to include a directory for initial FastQC reports, the processed sequences, final FastQC reports, and results for taxonomic profiling and phylogenetic tree construction all with new numbering.

```
./microSLURM_16S.sh -h

##############################################################
# microSLURM_16S                                             #
# 16S rRNA Gene Amplicon Sequencing Processing Pipeline      #
# Last updated: 10 Nov 2022                                  #
##############################################################
 
 Description: This is a wrapper program that wraps various  
 programs to process raw paired-end 16S amplicon sequences. 
 The end product of this pipeline are amplicon sequence     
 variant (ASV) abundances per sample along with             
 taxonomical classifications and phylogenetic tree of ASV   
 sequences. ASV tables produced by this pipeline should be  
 ready for further statistical analyses after performing    
 the appropriate data transformations.                      
                                                            
 Required programs and databases:                           
    SLURM:      Program is designed to work with a SLURM    
                high performance computing cluster          
                scheduling system.                          
    R base:     For performing various pipeline functions.  
    FastQC:     For performing initial quality reports.     
    MultiQC:    For summarizing multiple FastQC reports.    
    Cutadapt:   For removing PCR primers used to amplify    
                16S rRNA region.                            
    DADA2:      For quality trimming/filtering of reads,    
                ASV detection, and assignment of            
                ASV taxonomy. Taxonomy assignment requires  
                a DADA2 formatted taxonomy reference fasta  
                to be given.                                
    DECIPHER:   For generating a phylogenetic tree.         
                Performs the sequence alignment step of     
                creating the tree.                          
    phangorn:   For generating a phylogenetic tree.         
                Constructs the tree from multiple sequence  
                alignment that's been performed by DECIPHER.
                Required even when using RAxML-NG for       
                phylogenetic tree construction.             
    RAxML-NG:   Alternative method for generating a         
                phylogenetic tree. Constructs the tree from 
                multiple sequence alignment that's been     
                performed by DECIPHER, and starting tree    
                produced by NJ() function of phangorn.      
    Phyloseq:   For packaging together data that were       
                produced by the pipeline.                   
                                                            
 Usage:                                                     
 microSLURM_16S.sh -i input_seqs_dir \         
                    -o output_dir \                         
                    -p 'commands; to; load; programs' \     
                    -f notificationEmail@forFailures.edu \  
                    -c FORWARD_PRIMER,REVERSE_PRIMER        
                    -r path/to/taxonomy/reference.fa.gz \   
                    -n list,of,node,partitions \            
                    -t list,of,time,requests \              
                    -k list,of,cpu,number,requests \        
                    -m list,of,memory,per,cpu,requests \    
                    [additional options]                    
                                                            
 Parameters:                                                
     -h    Print the parameter list below then exit.        
                                                            
 Required analysis parameters                               
     -i    (Required) Directory that contains the raw       
           fastq files to be processed. Sequences must have 
           file extensions .fastq OR .fq,                   
           and can be gzipped or not.                       
     -o    (Required) Directory to put output of pipeline   
           into. NOTE: make sure output directory is in an  
           area that has plenty of data storage space       
           available if processing large datasets.          
     -p    (Required) Single quoted string that contains    
           commands to load all the necessary programs      
           needed to run pipeline steps (e.g. activating    
           conda environments, loading modules, adding to   
           PATH, etc.).                                     
     -f    (Required) E-mail to send notifications to upon  
           failure of any jobs.                             
     -c    (Required) Primer sequences used to amplify      
           the 16S rRNA region. Need to provide both        
           forward and reverse primers in a comma           
           separated list (no spaces). Make sure it is      
           in ACGT format. Wildcard nucleotides are         
           acceptable.                                      
     -r    (Required) Taxonomy reference file to use for    
           assigning taxonomic classifications. Provide     
           the full path to the file. Should be formatted   
           for use with DADA2. 
                                                            
 Required SLURM parameters                                  
 Note: For each parameter, a comma separated list of 8      
 entries must be given (one for each pipeline step).        
     -n    (Required) Names of the partitions to request for
           submitting jobs.                                 
     -t    (Required) Time requests for running jobs.       
           Specify in hours (e.g. 12:00:00 for 12 hours).   
     -k    (Required) Number of cores wanted for jobs.      
     -m    (Required) Amount of memory requested for each   
           requested core. Specify in Megabytes.            
     -f    (Required) E-mail to send notifications to upon  
           failure of any jobs.                             
                                                            
 Optional pipeline parameters                               
     -x    (Optional) How many nucleotide bases to trim off 
           the 5' end of reads during quality trimming/     
           filtering (as quality usually crashes towards the
           5' end). Default is 10.                          
     -a    (Optional) Add species taxonomy information      
           through exact matching (see DADA2 documentation  
           for how this is done). Provide full path to the  
           reference file for adding species taxonomy. If   
           none provided, species will not be added to      
           genus classifications. Requires the -t parameter 
           to be specified with a reference file that only  
           goes to genus level.                             
     -R    (Optional) Use RAxML-NG for phylogenetic tree    
           construction. Default is to use phangorn R       
           package.                                         
     -s    (Optional) Skip certain steps in the pipeline if 
           need be. Provide a comma separated list of steps 
           that you wish to skip in the pipeline. List may  
           have the values: fastqc_initial, cutadapt,       
           trim_filter, fastqc_final, taxonomic_profiling,  
           phylogenetic_tree.                               
                             
```
