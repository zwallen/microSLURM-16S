#!/bin/bash
set -e

##############################################################
# microSLURM_16S                                             #
# 16S rRNA Gene Amplicon Sequencing Processing Pipeline      #
# Last updated: 10 Nov 2022                                  #
#                                                            #
# Description: This is a wrapper program that wraps various  #
# programs to process raw paired-end 16S amplicon sequences. #
# The end product of this pipeline are amplicon sequence     #
# variant (ASV) abundances per sample along with             #
# taxonomic classifications and phylogenetic tree of ASV     #
# sequences. ASV tables produced by this pipeline should be  #
# ready for further statistical analyses after performing    #
# the appropriate data transformations.                      #
#                                                            #
# Required programs and databases:                           #
#    SLURM:      Program is designed to work with a SLURM    #
#                high performance computing cluster          #
#                scheduling system.                          #
#    R base:     For performing various pipeline functions.  #
#    FastQC:     For performing initial quality reports.     #
#    MultiQC:    For summarizing multiple FastQC reports.    #
#    Cutadapt:   For removing PCR primers used to amplify    #
#                16S rRNA region.                            #
#    DADA2:      For quality trimming/filtering of reads,    #
#                ASV detection, and assignment of            #
#                ASV taxonomy. Taxonomy assignment requires  #
#                a DADA2 formatted taxonomy reference fasta  #
#                to be given.                                #
#    DECIPHER:   For generating a phylogenetic tree.         #
#                Performs the sequence alignment step of     #
#                creating the tree.                          #
#    phangorn:   For generating a phylogenetic tree.         #
#                Constructs the tree from multiple sequence  #
#                alignment that's been performed by DECIPHER.#
#                Required even when using RAxML-NG for       #
#                phylogenetic tree construction.             #
#    RAxML-NG:   Alternative method for generating a         #
#                phylogenetic tree. Constructs the tree from #
#                multiple sequence alignment that's been     #
#                performed by DECIPHER, and starting tree    #
#                produced by NJ() function of phangorn.      #
#    Phyloseq:   For packaging together data that were       #
#                produced by the pipeline.                   #
#                                                            #
# Usage:                                                     #
# microSLURM_16S.sh -i input_seqs_dir \                      #
#                    -o output_dir \                         #
#                    -p 'commands; to; load; programs' \     #
#                    -f notificationEmail@forFailures.edu \  #
#                    -c FORWARD_PRIMER,REVERSE_PRIMER \      #
#                    -r path/to/taxonomy/reference.fa.gz \   #
#                    -n list,of,node,partitions \            #
#                    -t list,of,time,requests \              #
#                    -k list,of,cpu,number,requests \        #
#                    -m list,of,memory,per,cpu,requests \    #
#                    [additional options]                    #
#                                                            #
# Parameters:                                                #
#     -h    Print the parameter list below then exit.        #
#                                                            #
# Required analysis parameters                               #
#     -i    (Required) Directory that contains the raw       #
#           fastq files to be processed. Sequences must have #
#           file extensions .fastq OR .fq,                   #
#           and can be gzipped or not.                       #
#     -o    (Required) Directory to put output of pipeline   #
#           into. NOTE: make sure output directory is in an  #
#           area that has plenty of data storage space       #
#           available if processing large datasets.          #
#     -p    (Required) Single quoted string that contains    #
#           commands to load all the necessary programs      #
#           needed to run pipeline steps (e.g. activating    #
#           conda environments, loading modules, adding to   #
#           PATH, etc.).                                     #
#     -c    (Required) Primer sequences used to amplify      #
#           the 16S rRNA gene region. Need to provide both   #
#           forward and reverse primers in a comma           #
#           separated list (no spaces). Make sure it is      #
#           in ACGT format. Wildcard nucleotides are         #
#           acceptable. List forward then reverse primer.    #
#     -r    (Required) Taxonomy reference file to use for    #
#           assigning taxonomic classifications. Provide     #
#           the full path to the file. Should be formatted   #
#           for use with DADA2.                              #
#                                                            #
# Required SLURM parameters                                  #
# Note: For each parameter, a comma separated list of 8      #
# entries must be given (one for each pipeline step).        #
#     -n    (Required) Names of the partitions to request for#
#           submitting jobs.                                 #
#     -t    (Required) Time requests for running jobs.       #
#           Specify in hours (e.g. 12:00:00 for 12 hours).   #
#     -k    (Required) Number of cores wanted for jobs.      #
#     -m    (Required) Amount of memory requested for each   #
#           requested core. Specify in Megabytes.            #
#     -f    (Required) E-mail to send notifications to upon  #
#           failure of any jobs.                             #
#                                                            #
# Optional pipeline parameters                               #
#     -x    (Optional) How many nucleotide bases to trim off #
#           the 5' end of reads during quality trimming/     #
#           filtering (as quality usually crashes towards the#
#           5' end). Default is 10.                          #
#     -a    (Optional) Add species taxonomy information      #
#           through exact matching (see DADA2 documentation  #
#           for how this is done). Provide full path to the  #
#           reference file for adding species taxonomy. If   #
#           none provided, species will not be added to      #
#           genus classifications. Requires the -t parameter #
#           to be specified with a reference file that only  #
#           goes to genus level.                             #
#     -R    (Optional) Use RAxML-NG for phylogenetic tree    #
#           construction. Default is to use phangorn R       #
#           package.                                         #
#     -s    (Optional) Skip certain steps in the pipeline if #
#           need be. Provide a comma separated list of steps #
#           that you wish to skip in the pipeline. List may  #
#           have the values: fastqc_initial, cutadapt,       #
#           trim_filter, fastqc_final, taxonomic_profiling,  #
#           phylogenetic_tree.                               #
#     -d    (Optional) Delete sequence files generated during#
#           intermediate pipeline steps (i.e. cutadapt).     #
#           Will keep log files and reorganize output folder #
#           if flag is give since default numbering of folder#
#           will no longer make sense.                       #
##############################################################

echo " "
echo "##############################################################"
echo "# microSLURM_16S                                             #"
echo "# 16S rRNA Gene Amplicon Sequencing Processing Pipeline      #"
echo "# Last updated: 10 Nov 2022                                  #"
echo "##############################################################"
echo " "

# Argument parsing
while getopts ":hi:o:p:c:x:r:n:t:k:m:f:a:Rs:d" opt; do
  case $opt in
    h)
    echo " Description: This is a wrapper program that wraps various  "
    echo " programs to process raw paired-end 16S amplicon sequences. "
    echo " The end product of this pipeline are amplicon sequence     "
    echo " variant (ASV) abundances per sample along with             "
    echo " taxonomical classifications and phylogenetic tree of ASV   "
    echo " sequences. ASV tables produced by this pipeline should be  "
    echo " ready for further statistical analyses after performing    "
    echo " the appropriate data transformations.                      "
    echo "                                                            "
    echo " Required programs and databases:                           "
    echo "    SLURM:      Program is designed to work with a SLURM    "
    echo "                high performance computing cluster          "
    echo "                scheduling system.                          "
    echo "    R base:     For performing various pipeline functions.  "
    echo "    FastQC:     For performing initial quality reports.     "
    echo "    MultiQC:    For summarizing multiple FastQC reports.    "
    echo "    Cutadapt:   For removing PCR primers used to amplify    "
    echo "                16S rRNA region.                            "
    echo "    DADA2:      For quality trimming/filtering of reads,    "
    echo "                ASV detection, and assignment of            "
    echo "                ASV taxonomy. Taxonomy assignment requires  "
    echo "                a DADA2 formatted taxonomy reference fasta  "
    echo "                to be given.                                "
    echo "    DECIPHER:   For generating a phylogenetic tree.         "
    echo "                Performs the sequence alignment step of     "
    echo "                creating the tree.                          "
    echo "    phangorn:   For generating a phylogenetic tree.         "
    echo "                Constructs the tree from multiple sequence  "
    echo "                alignment that's been performed by DECIPHER."
    echo "                Required even when using RAxML-NG for       "
    echo "                phylogenetic tree construction.             "
    echo "    RAxML-NG:   Alternative method for generating a         "
    echo "                phylogenetic tree. Constructs the tree from "
    echo "                multiple sequence alignment that's been     "
    echo "                performed by DECIPHER, and starting tree    "
    echo "                produced by NJ() function of phangorn.      "
    echo "    Phyloseq:   For packaging together data that were       "
    echo "                produced by the pipeline.                   "
    echo "                                                            "
    echo " Usage:                                                     "
    echo " microSLURM_16S.sh -i input_seqs_dir \                      "
    echo "                    -o output_dir \                         "
    echo "                    -p 'commands; to; load; programs' \     "
    echo "                    -f notificationEmail@forFailures.edu \  "
    echo "                    -c FORWARD_PRIMER,REVERSE_PRIMER        "
    echo "                    -r path/to/taxonomy/reference.fa.gz \   "
    echo "                    -n list,of,node,partitions \            "
    echo "                    -t list,of,time,requests \              "
    echo "                    -k list,of,cpu,number,requests \        "
    echo "                    -m list,of,memory,per,cpu,requests \    "
    echo "                    [additional options]                    "
    echo "                                                            "
    echo " Parameters:                                                "
    echo "     -h    Print the parameter list below then exit.        "
    echo "                                                            "
    echo " Required analysis parameters                               "
    echo "     -i    (Required) Directory that contains the raw       "
    echo "           fastq files to be processed. Sequences must have "
    echo "           file extensions .fastq OR .fq,                   "
    echo "           and can be gzipped or not.                       "
    echo "     -o    (Required) Directory to put output of pipeline   "
    echo "           into. NOTE: make sure output directory is in an  "
    echo "           area that has plenty of data storage space       "
    echo "           available if processing large datasets.          "
    echo "     -p    (Required) Single quoted string that contains    "
    echo "           commands to load all the necessary programs      "
    echo "           needed to run pipeline steps (e.g. activating    "
    echo "           conda environments, loading modules, adding to   "
    echo "           PATH, etc.).                                     "
    echo "     -f    (Required) E-mail to send notifications to upon  "
    echo "           failure of any jobs.                             "
    echo "     -c    (Required) Primer sequences used to amplify      "
    echo "           the 16S rRNA region. Need to provide both        "
    echo "           forward and reverse primers in a comma           "
    echo "           separated list (no spaces). Make sure it is      "
    echo "           in ACGT format. Wildcard nucleotides are         "
    echo "           acceptable.                                      "
    echo "     -r    (Required) Taxonomy reference file to use for    "
    echo "           assigning taxonomic classifications. Provide     "
    echo "           the full path to the file. Should be formatted   "
    echo "           for use with DADA2.                              "
    echo "                                                            "
    echo " Required SLURM parameters                                  "
    echo " Note: For each parameter, a comma separated list of 8      "
    echo " entries must be given (one for each pipeline step).        "
    echo "     -n    (Required) Names of the partitions to request for"
    echo "           submitting jobs.                                 "
    echo "     -t    (Required) Time requests for running jobs.       "
    echo "           Specify in hours (e.g. 12:00:00 for 12 hours).   "
    echo "     -k    (Required) Number of cores wanted for jobs.      "
    echo "     -m    (Required) Amount of memory requested for each   "
    echo "           requested core. Specify in Megabytes.            "
    echo "     -f    (Required) E-mail to send notifications to upon  "
    echo "           failure of any jobs.                             "
    echo "                                                            "
    echo " Optional pipeline parameters                               "
    echo "     -x    (Optional) How many nucleotide bases to trim off "
    echo "           the 5' end of reads during quality trimming/     "
    echo "           filtering (as quality usually crashes towards the"
    echo "           5' end). Default is 10.                          "
    echo "     -a    (Optional) Add species taxonomy information      "
    echo "           through exact matching (see DADA2 documentation  "
    echo "           for how this is done). Provide full path to the  "
    echo "           reference file for adding species taxonomy. If   "
    echo "           none provided, species will not be added to      "
    echo "           genus classifications. Requires the -t parameter "
    echo "           to be specified with a reference file that only  "
    echo "           goes to genus level.                             "
    echo "     -R    (Optional) Use RAxML-NG for phylogenetic tree    "
    echo "           construction. Default is to use phangorn R       "
    echo "           package.                                         "
    echo "     -s    (Optional) Skip certain steps in the pipeline if "
    echo "           need be. Provide a comma separated list of steps "
    echo "           that you wish to skip in the pipeline. List may  "
    echo "           have the values: fastqc_initial, cutadapt,       "
    echo "           trim_filter, fastqc_final, taxonomic_profiling,  "
    echo "           phylogenetic_tree.                               "
    echo "     -d    (Optional) Delete sequence files generated during"
    echo "           intermediate pipeline steps (i.e. cutadapt).     "
    echo "           Will keep log files and reorganize output folder "
    echo "           if flag is give since default numbering of folder"
    echo "           will no longer make sense.                       "
    echo " "
    exit 0
    ;;
    i) SEQ_DIR=$(echo $OPTARG | sed 's#/$##')
    ;;
    o) OUT_DIR=$(echo $OPTARG | sed 's#/$##')
    ;;
    p) PROG_LOAD="$OPTARG"
    ;;
    c) PRIMERS="$OPTARG"
    ;;
    x) TRIM="$OPTARG"
    ;;
    r) TAXA_REF="$OPTARG"
    ;;
    n) PARTITION="$OPTARG"
    ;;
    t) TIME_REQUEST="$OPTARG"
    ;;
    k) CPU_REQUEST="$OPTARG"
    ;;
    m) MEM_PER_CPU="$OPTARG"
    ;;
    f) FAIL_EMAIL="$OPTARG"
    ;;
    a) SPEC_REF="$OPTARG"
    ;;
    R) RAXML=1
    ;;
    s) SKIP="$OPTARG"
    ;;
    d) DELETE=1
    ;;
    \?) echo "Invalid option: $OPTARG" 1>&2
        exit 1
    ;;
    :) echo "Invalid option: $OPTARG requires an argument" 1>&2
       exit 1
    ;;
  esac
done

# Check that valid arguments were entered

# -i
if [[ -z "$SEQ_DIR" ]]; then
  echo "ERROR: Argument -i is required, please supply a directory with input fastq files"
  exit 1
fi
if [[ ! -d "$SEQ_DIR" ]]; then
  echo "ERROR: Argument -i should be a directory, please supply a directory with input fastq files"
  exit 1
fi
if ls -l $SEQ_DIR | grep -q ".fastq.gz"; then
  SEQ_EXT=fastq.gz
elif ls -l $SEQ_DIR | grep -q ".fastq"; then
  SEQ_EXT=fastq
elif ls -l $SEQ_DIR | grep -q ".fq.gz"; then
  SEQ_EXT=fq.gz
elif ls -l $SEQ_DIR | grep -q ".fq"; then
  SEQ_EXT=fq
else
  echo "ERROR: Sequences in input directory should have file extension of either .fastq[.gz] OR .fq[.gz]"
  exit 1
fi
found_1=$(ls -l $SEQ_DIR | grep -q "_R1_001")
found_2=$(ls -l $SEQ_DIR | grep -q "_R2_001")
if [[ -n "$found_1" ]] && [[ -n "$found_2" ]]; then
  echo "ERROR: Sequences in input directory expected to be paired-end Illumina sequence fastq files whose file names contain the strings '_R1_001' and '_R2_001'"
  exit 1
else
  :
fi

# -o
if [[ -z "$OUT_DIR" ]]; then
  echo "ERROR: Argument -o is required, please supply an output directory"
  exit 1
fi
if [[ ! -d "$OUT_DIR" ]]; then
  echo "ERROR: Argument -o should be a directory, please supply an output directory"
  exit 1
fi

# -p
if [[ -z "$PROG_LOAD" ]]; then
  echo "ERROR: Argument -p is required, please supply a single quoted string of commands needed to load required programs (can be an empty string ' ' if none required)"
  exit 1
fi

# -c
if [[ -z "$PRIMERS" ]]; then
  echo "ERROR: Argument -c is required, please supply a comma separated list specifying forward and reverse PCR primers"
  exit 1
fi
if echo $PRIMERS | grep -q ","; then
  :
else
  echo "ERROR: Invalid input given to -c, should be a comma separated list specifying forward and reverse PCR primers"
  exit 1
fi
if echo $PRIMERS | grep -q " "; then
  echo "ERROR: Invalid input given to -c, should be a comma separated list specifying forward and reverse PCR primers with no spaces"
  exit 1
fi
if [[ $(echo $PRIMERS | awk -F"," '{print NF}') -eq 2 ]]; then
  :
else
  echo "ERROR: Invalid input given to -c, should be a comma separated list specifying only two PCR primers (forward and reverse)"
  exit 1
fi
FWD_PRIMER=$(echo $PRIMERS | awk -F"," '{print $1}')
REV_PRIMER=$(echo $PRIMERS | awk -F"," '{print $2}')

# -x
if [[ ! -z "$TRIM" ]]; then
  number='^[0-9]+$'
  if [[ ! "$TRIM" =~ $number ]]; then
    echo "ERROR: Please supply an integer >= to 0 for parameter -x"
    exit 1
  fi
  if [[ "$TRIM" -lt 0 ]]; then
    echo "ERROR: Please supply an integer >= to 0 for parameter -x"
    exit 1
  fi
fi
if [[ -z "$TRIM" ]]; then
  TRIM=10
fi

# -r
if [[ -z "$TAXA_REF" ]]; then
  echo "ERROR: Argument -t is required, please supply a DADA2 formatted reference file to use for taxonomic assignment of ASVs"
  exit 1
fi
if [[ -d "$TAXA_REF" ]]; then
  echo "ERROR: Argument -t should be the path to a single file, not a directory, please supply path to a taxonomic reference file"
  exit 1
fi

# -n
if [[ -z "$PARTITION" ]]; then
  echo "ERROR: Argument -n is required, please supply a node partition name to send jobs to."
  exit 1
fi
PARTITION_1=$(echo $PARTITION | awk -F',' '{print $1}')
PARTITION_2=$(echo $PARTITION | awk -F',' '{print $2}')
PARTITION_3=$(echo $PARTITION | awk -F',' '{print $3}')
PARTITION_4=$(echo $PARTITION | awk -F',' '{print $4}')
PARTITION_5=$(echo $PARTITION | awk -F',' '{print $5}')

# -t
if [[ -z "$TIME_REQUEST" ]]; then
  echo "ERROR: Argument -t is required, please supply a max length of time to run each job for."
  exit 1
fi
TIME_1=$(echo $TIME_REQUEST | awk -F',' '{print $1}')
TIME_2=$(echo $TIME_REQUEST | awk -F',' '{print $2}')
TIME_3=$(echo $TIME_REQUEST | awk -F',' '{print $3}')
TIME_4=$(echo $TIME_REQUEST | awk -F',' '{print $4}')
TIME_5=$(echo $TIME_REQUEST | awk -F',' '{print $5}')

# -k
if [[ -z "$CPU_REQUEST" ]]; then
  echo "ERROR: Argument -k is required, please supply the number of cores being requested to run each job."
  exit 1
fi
CPU_1=$(echo $CPU_REQUEST | awk -F',' '{print $1}')
CPU_2=$(echo $CPU_REQUEST | awk -F',' '{print $2}')
CPU_3=$(echo $CPU_REQUEST | awk -F',' '{print $3}')
CPU_4=$(echo $CPU_REQUEST | awk -F',' '{print $4}')
CPU_5=$(echo $CPU_REQUEST | awk -F',' '{print $5}')

# -m
if [[ -z "$MEM_PER_CPU" ]]; then
  echo "ERROR: Argument -m is required, please supply a memory request for each core of each job."
  exit 1
fi
MEM_1=$(echo $MEM_PER_CPU | awk -F',' '{print $1}')
MEM_2=$(echo $MEM_PER_CPU | awk -F',' '{print $2}')
MEM_3=$(echo $MEM_PER_CPU | awk -F',' '{print $3}')
MEM_4=$(echo $MEM_PER_CPU | awk -F',' '{print $4}')
MEM_5=$(echo $MEM_PER_CPU | awk -F',' '{print $5}')

# -f
if [[ -z "$FAIL_EMAIL" ]]; then
  echo "ERROR: Argument -f is required, please supply an email that can be notified upon failure of any jobs ran during the pipeline"
  exit 1
elif echo $FAIL_EMAIL | grep -q -v '@'; then
  echo "ERROR: Argument -f requires a valid email, please give an email in the form of xxxx@xxxx.xxx"
fi

# -a
if [[ ! -z "$SPEC_REF" ]]; then
  if [[ -d "$SPEC_REF" ]]; then
    echo "ERROR: Argument -a should be the path to a single file, not a directory, please supply path to a species taxonomic reference file"
    exit 1
  fi
fi

# -s
if [[ ! -z "$SKIP" ]]; then
  if echo $SKIP | grep -q "fastqc_initial"; then
    :
  elif echo $SKIP | grep -q "cutadapt"; then
    :
  elif echo $SKIP | grep -q "trim_filter"; then
    :
  elif echo $SKIP | grep -q "fastqc_final"; then
    :
  elif echo $SKIP | grep -q "taxonomic_profiling"; then
    :
  elif echo $SKIP | grep -q "phylogenetic_tree"; then
    :
  else
    echo "ERROR: Invalid argument given to -s, please specify one or more of: fastqc_initial, cutadapt, trim_filter, fastqc_final, taxonomic_profiling, phylogenetic_tree"
    exit 1
  fi
fi

###### CREATE DIRECTORY FOR PIPELINE OUTPUT #####
DATE=$(date | awk '{print $3"_"$2"_"$6}')
if [ -d "${OUT_DIR}/16S_Amplicon_Pipeline_${DATE}" ]
then
	 :
else
	 mkdir ${OUT_DIR}/16S_Amplicon_Pipeline_${DATE}
fi

echo " "
echo "Directory for 16S amplicon pipeline output: 16S_Amplicon_Pipeline_${DATE}"
echo " "
RESULTS_DIR="${OUT_DIR}/16S_Amplicon_Pipeline_${DATE}"

############# INITIAL FASTQC REPORT #############
if echo $SKIP | grep -q "fastqc_initial"; then
  echo " "
  echo "*** Skipping running of initial FastQC report generation on input fastq files ***"
  echo " "

  #Create directory
  if [ -d "${RESULTS_DIR}/1.FastQC_Initial_Reports" ]
  then
	  :
  else
	  mkdir ${RESULTS_DIR}/1.FastQC_Initial_Reports
	  mkdir ${RESULTS_DIR}/1.FastQC_Initial_Reports/0.ErrorOut
	  mkdir ${RESULTS_DIR}/1.FastQC_Initial_Reports/0.Output
  fi
else
  SECONDS=0
  echo " "
  echo "*** Running FastQC on input fastq files ***"
  echo " "

  #Create directory for output
  if [ -d "${RESULTS_DIR}/1.FastQC_Initial_Reports" ]
  then
	  :
  else
	  mkdir ${RESULTS_DIR}/1.FastQC_Initial_Reports
	  mkdir ${RESULTS_DIR}/1.FastQC_Initial_Reports/0.ErrorOut
	  mkdir ${RESULTS_DIR}/1.FastQC_Initial_Reports/0.Output
  fi

  ##### Run FastQC #####
  #Create script for running program and submit
  echo '#!/bin/bash' > run.sh
  echo "#SBATCH --partition=$PARTITION_1" >> run.sh
  echo "#SBATCH --job-name=FastQC" >> run.sh
  echo "#SBATCH --error=${RESULTS_DIR}/1.FastQC_Initial_Reports/0.ErrorOut/FastQC_%A_%a.err" >> run.sh
  echo "#SBATCH --output=${RESULTS_DIR}/1.FastQC_Initial_Reports/0.Output/FastQC_%A_%a.out" >> run.sh
  echo "#SBATCH --time=$TIME_1" >> run.sh
  echo "#SBATCH --ntasks=1" >> run.sh
  echo "#SBATCH --cpus-per-task=$CPU_1" >> run.sh
  echo "#SBATCH --mem-per-cpu=$MEM_1" >> run.sh
  echo "#SBATCH --mail-type=FAIL" >> run.sh
  echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> run.sh
  echo "#SBATCH --array=1-$(ls -l ${SEQ_DIR}/*R1_001.${SEQ_EXT} | wc -l)" >> run.sh
  echo "#SBATCH --wait" >> run.sh
  echo "$PROG_LOAD" >> run.sh
  echo "FILE1=\$(ls ${SEQ_DIR}/*R1_001.${SEQ_EXT} | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
  echo "FILE2=\$(ls ${SEQ_DIR}/*R2_001.${SEQ_EXT} | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
  echo "FILE_NAME=\$(echo \$FILE1 | awk -F '/' '{print \$NF}' | awk -F '_R1_001' '{print \$1}')" >> run.sh
  echo "fastqc \$FILE1 \$FILE2 -d ${RESULTS_DIR}/1.FastQC_Initial_Reports -o ${RESULTS_DIR}/1.FastQC_Initial_Reports \\" >> run.sh
  echo "> ${RESULTS_DIR}/1.FastQC_Initial_Reports/\${FILE_NAME}.log 2>&1" >> run.sh
  chmod +x run.sh

  sbatch run.sh > /dev/null

  echo '#!/bin/bash' > run.sh
  echo "#SBATCH --partition=$PARTITION_1" >> run.sh
  echo "#SBATCH --job-name=MultiQC" >> run.sh
  echo "#SBATCH --error=${RESULTS_DIR}/1.FastQC_Initial_Reports/0.ErrorOut/MultiQC.err" >> run.sh
  echo "#SBATCH --output=${RESULTS_DIR}/1.FastQC_Initial_Reports/0.Output/MultiQC.out" >> run.sh
  echo "#SBATCH --time=$TIME_1" >> run.sh
  echo "#SBATCH --ntasks=1" >> run.sh
  echo "#SBATCH --cpus-per-task=$CPU_1" >> run.sh
  echo "#SBATCH --mem-per-cpu=$MEM_1" >> run.sh
  echo "#SBATCH --mail-type=FAIL" >> run.sh
  echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> run.sh
  echo "#SBATCH --wait" >> run.sh
  echo "$PROG_LOAD" >> run.sh
  echo "multiqc ${RESULTS_DIR}/1.FastQC_Initial_Reports -o ${RESULTS_DIR}/1.FastQC_Initial_Reports \\" >> run.sh
  echo "> ${RESULTS_DIR}/1.FastQC_Initial_Reports/multiqc.log 2>&1" >> run.sh
  chmod +x run.sh

  sbatch run.sh > /dev/null

  rm run.sh

  #Signal jobs have ended
  echo "Initial FastQC reports complete"
  echo "Elapsed time: $(($SECONDS / 3600)) hr : $(($(($SECONDS % 3600)) / 60)) min : $(($SECONDS % 60)) sec"
  echo " "
fi
#################################################

######### REMOVE PRIMERS WITH CUTADAPT ##########
if echo $SKIP | grep -q "cutadapt"; then
  echo " "
  echo "*** Skipping running of cutadapt for primer removal ***"
  echo " "
  
  #Create directory
  if [ -d "${RESULTS_DIR}/2.Primer_Trimmed_Sequences" ]
  then
    :
  else
    mkdir ${RESULTS_DIR}/2.Primer_Trimmed_Sequences
    mkdir ${RESULTS_DIR}/2.Primer_Trimmed_Sequences/0.ErrorOut
    mkdir ${RESULTS_DIR}/2.Primer_Trimmed_Sequences/0.Output
  fi
else
  SECONDS=0
  echo " "
  echo "*** Running cutadapt on input fastq files to remove primer sequences ***"
  echo " "
  
  #Create directory for output
  if [ -d "${RESULTS_DIR}/2.Primer_Trimmed_Sequences" ]
  then
    :
  else
    mkdir ${RESULTS_DIR}/2.Primer_Trimmed_Sequences
    mkdir ${RESULTS_DIR}/2.Primer_Trimmed_Sequences/0.ErrorOut
    mkdir ${RESULTS_DIR}/2.Primer_Trimmed_Sequences/0.Output
  fi
  
  ##### Run cutadapt #####
  #Create script for running program and submit
  echo '#!/bin/bash' > run.sh
  echo "#SBATCH --partition=$PARTITION_2" >> run.sh
  echo "#SBATCH --job-name=cutadapt" >> run.sh
  echo "#SBATCH --error=${RESULTS_DIR}/2.Primer_Trimmed_Sequences/0.ErrorOut/Cutadapt_%A_%a.err" >> run.sh
  echo "#SBATCH --output=${RESULTS_DIR}/2.Primer_Trimmed_Sequences/0.Output/Cutadapt_%A_%a.out" >> run.sh
  echo "#SBATCH --time=$TIME_2" >> run.sh
  echo "#SBATCH --ntasks=1" >> run.sh
  echo "#SBATCH --cpus-per-task=$CPU_2" >> run.sh
  echo "#SBATCH --mem-per-cpu=$MEM_2" >> run.sh
  echo "#SBATCH --mail-type=FAIL" >> run.sh
  echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> run.sh
  echo "#SBATCH --array=1-$(ls -l ${SEQ_DIR}/*R1_001.${SEQ_EXT} | wc -l)" >> run.sh
  echo "#SBATCH --wait" >> run.sh
  echo "$PROG_LOAD" >> run.sh
  echo "FILE1=\$(ls ${SEQ_DIR}/*R1_001.${SEQ_EXT} | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
  echo "FILE2=\$(ls ${SEQ_DIR}/*R2_001.${SEQ_EXT} | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
  echo "FILE_NAME=\$(echo \$FILE1 | awk -F '/' '{print \$NF}' | awk -F '_R1_001' '{print \$1}')" >> run.sh
  echo "cutadapt -b $FWD_PRIMER -B $REV_PRIMER \\" >> run.sh
  echo "-o ${RESULTS_DIR}/2.Primer_Trimmed_Sequences/\${FILE_NAME}_R1_001.fastq.gz \\" >> run.sh
  echo "-p ${RESULTS_DIR}/2.Primer_Trimmed_Sequences/\${FILE_NAME}_R2_001.fastq.gz \\" >> run.sh
  echo "\$FILE1 \$FILE2 \\" >> run.sh
  echo "> ${RESULTS_DIR}/2.Primer_Trimmed_Sequences/\${FILE_NAME}.log 2>&1" >> run.sh
  chmod +x run.sh
  
  sbatch run.sh > /dev/null
  
  rm run.sh
  
  #Signal jobs have ended
  echo "Removal of primer sequences with cutadapt complete"
  echo "Elapsed time: $(($SECONDS / 3600)) hr : $(($(($SECONDS % 3600)) / 60)) min : $(($SECONDS % 60)) sec"
  echo " "
fi
#################################################

##### QUALITY TRIMMING/FILTERING WITH DADA2 #####
if echo $SKIP | grep -q "trim_filter"; then
  echo " "
  echo "*** Skipping quality trimming/filtering of sequences ***"
  echo " "
  
  #Create directory
  if [ -d "${RESULTS_DIR}/3.Quality_Controlled_Sequences" ]
  then
    :
  else
    mkdir ${RESULTS_DIR}/3.Quality_Controlled_Sequences
    mkdir ${RESULTS_DIR}/3.Quality_Controlled_Sequences/0.ErrorOut
    mkdir ${RESULTS_DIR}/3.Quality_Controlled_Sequences/0.Output
  fi
else
  SECONDS=0
  echo " "
  echo "*** Running filterAndTrim DADA2 function for quality trimming/filtering sequences ***"
  echo " "
  
  #Create directory for output
  if [ -d "${RESULTS_DIR}/3.Quality_Controlled_Sequences" ]
  then
    :
  else
    mkdir ${RESULTS_DIR}/3.Quality_Controlled_Sequences
    mkdir ${RESULTS_DIR}/3.Quality_Controlled_Sequences/0.ErrorOut
    mkdir ${RESULTS_DIR}/3.Quality_Controlled_Sequences/0.Output
  fi
  
  ##### Run filterAndTrim #####
  #Create script for running program and submit
  echo '#!/bin/bash' > run.sh
  echo "#SBATCH --partition=$PARTITION_3" >> run.sh
  echo "#SBATCH --job-name=filterAndTrim" >> run.sh
  echo "#SBATCH --error=${RESULTS_DIR}/3.Quality_Controlled_Sequences/0.ErrorOut/filterAndTrim.err" >> run.sh
  echo "#SBATCH --output=${RESULTS_DIR}/3.Quality_Controlled_Sequences/0.Output/filterAndTrim.out" >> run.sh
  echo "#SBATCH --time=$TIME_3" >> run.sh
  echo "#SBATCH --ntasks=1" >> run.sh
  echo "#SBATCH --cpus-per-task=$CPU_3" >> run.sh
  echo "#SBATCH --mem-per-cpu=$MEM_3" >> run.sh
  echo "#SBATCH --mail-type=FAIL" >> run.sh
  echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> run.sh
  echo "#SBATCH --wait" >> run.sh
  echo "$PROG_LOAD" >> run.sh
  echo "echo \"library(dada2)\" > run.R" >> run.sh
  echo "echo \"seqpath <- '${RESULTS_DIR}/2.Primer_Trimmed_Sequences'\" >> run.R" >> run.sh
  echo "echo \"filtpath <- '${RESULTS_DIR}/3.Quality_Controlled_Sequences'\" >> run.R" >> run.sh
  echo "echo \"fastqFs <- sort(list.files(seqpath, pattern='_R1_001.fastq.gz'))\" >> run.R" >> run.sh
  echo "echo \"fastqRs <- sort(list.files(seqpath, pattern='_R2_001.fastq.gz'))\" >> run.R" >> run.sh
  echo "echo \"cat('\n', 'Performing quality control with filterAndTrim function from DADA2:', '\n')\" >> run.R" >> run.sh
  echo "echo \"out <- filterAndTrim(fwd=file.path(seqpath, fastqFs), filt=file.path(filtpath, fastqFs), \" >> run.R" >> run.sh
  echo "echo \"                     rev=file.path(seqpath, fastqRs), filt.rev=file.path(filtpath, fastqRs), \" >> run.R" >> run.sh
  echo "echo \"                     maxEE=2, truncQ=2, maxN=0, minLen=50, trimRight=${TRIM}, rm.phix=T, n=1e5, compress=T, multithread=${CPU_3}, verbose=T)\" >> run.R" >> run.sh
  echo "echo \"out\" >> run.R" >> run.sh
  echo "echo \"cat('\n', 'Average proportion of reads remaining after quality trimming/filtering:',round(mean(out[,2]/out[,1]),2),'\n')\" >> run.R" >> run.sh
  echo "Rscript run.R > ${RESULTS_DIR}/3.Quality_Controlled_Sequences/0.filterAndTrim.log 2>&1" >> run.sh
  echo "rm run.R" >> run.sh
  chmod +x run.sh
  
  sbatch run.sh > /dev/null
  
  rm run.sh
  
  #Signal jobs have ended
  echo "Quality trimming/filtering with DADA2's filterAndTrim complete"
  echo "Elapsed time: $(($SECONDS / 3600)) hr : $(($(($SECONDS % 3600)) / 60)) min : $(($SECONDS % 60)) sec"
  echo " "
fi
#################################################

############# FINAL FASTQC REPORT #############
if echo $SKIP | grep -q "fastqc_final"; then
  echo " "
  echo "*** Skipping running of final FastQC report generation on quality controlled fastq files ***"
  echo " "

  #Create directory
  if [ -d "${RESULTS_DIR}/4.FastQC_Final_Reports" ]
  then
	  :
  else
	  mkdir ${RESULTS_DIR}/4.FastQC_Final_Reports
	  mkdir ${RESULTS_DIR}/4.FastQC_Final_Reports/0.ErrorOut
	  mkdir ${RESULTS_DIR}/4.FastQC_Final_Reports/0.Output
  fi
else
  SECONDS=0
  echo " "
  echo "*** Running final FastQC report on QCed fastq files ***"
  echo " "

  #Create directory for output
  if [ -d "${RESULTS_DIR}/4.FastQC_Final_Reports" ]
  then
	  :
  else
	  mkdir ${RESULTS_DIR}/4.FastQC_Final_Reports
	  mkdir ${RESULTS_DIR}/4.FastQC_Final_Reports/0.ErrorOut
	  mkdir ${RESULTS_DIR}/4.FastQC_Final_Reports/0.Output
  fi

  ##### Run FastQC #####
  #Create script for running program and submit
  echo '#!/bin/bash' > run.sh
  echo "#SBATCH --partition=$PARTITION_4" >> run.sh
  echo "#SBATCH --job-name=FastQC" >> run.sh
  echo "#SBATCH --error=${RESULTS_DIR}/4.FastQC_Final_Reports/0.ErrorOut/FastQC_%A_%a.err" >> run.sh
  echo "#SBATCH --output=${RESULTS_DIR}/4.FastQC_Final_Reports/0.Output/FastQC_%A_%a.out" >> run.sh
  echo "#SBATCH --time=$TIME_4" >> run.sh
  echo "#SBATCH --ntasks=1" >> run.sh
  echo "#SBATCH --cpus-per-task=$CPU_4" >> run.sh
  echo "#SBATCH --mem-per-cpu=$MEM_4" >> run.sh
  echo "#SBATCH --mail-type=FAIL" >> run.sh
  echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> run.sh
  echo "#SBATCH --array=1-$(ls -l ${RESULTS_DIR}/3.Quality_Controlled_Sequences/*R1_001.fastq.gz | wc -l)" >> run.sh
  echo "#SBATCH --wait" >> run.sh
  echo "$PROG_LOAD" >> run.sh
  echo "FILE1=\$(ls ${RESULTS_DIR}/3.Quality_Controlled_Sequences/*R1_001.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
  echo "FILE2=\$(ls ${RESULTS_DIR}/3.Quality_Controlled_Sequences/*R2_001.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
  echo "FILE_NAME=\$(echo \$FILE1 | awk -F '/' '{print \$NF}' | awk -F '_R1_001' '{print \$1}')" >> run.sh
  echo "fastqc \$FILE1 \$FILE2 -d ${RESULTS_DIR}/4.FastQC_Final_Reports -o ${RESULTS_DIR}/4.FastQC_Final_Reports \\" >> run.sh
  echo "> ${RESULTS_DIR}/4.FastQC_Final_Reports/\${FILE_NAME}.log 2>&1" >> run.sh
  chmod +x run.sh

  sbatch run.sh > /dev/null

  echo '#!/bin/bash' > run.sh
  echo "#SBATCH --partition=$PARTITION_4" >> run.sh
  echo "#SBATCH --job-name=MultiQC" >> run.sh
  echo "#SBATCH --error=${RESULTS_DIR}/4.FastQC_Final_Reports/0.ErrorOut/MultiQC.err" >> run.sh
  echo "#SBATCH --output=${RESULTS_DIR}/4.FastQC_Final_Reports/0.Output/MultiQC.out" >> run.sh
  echo "#SBATCH --time=$TIME_4" >> run.sh
  echo "#SBATCH --ntasks=1" >> run.sh
  echo "#SBATCH --cpus-per-task=$CPU_4" >> run.sh
  echo "#SBATCH --mem-per-cpu=$MEM_4" >> run.sh
  echo "#SBATCH --mail-type=FAIL" >> run.sh
  echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> run.sh
  echo "#SBATCH --wait" >> run.sh
  echo "$PROG_LOAD" >> run.sh
  echo "multiqc ${RESULTS_DIR}/4.FastQC_Final_Reports -o ${RESULTS_DIR}/4.FastQC_Final_Reports \\" >> run.sh
  echo "> ${RESULTS_DIR}/4.FastQC_Final_Reports/multiqc.log 2>&1" >> run.sh
  chmod +x run.sh

  sbatch run.sh > /dev/null

  rm run.sh

  #Signal jobs have ended
  echo "Final FastQC reports complete"
  echo "Elapsed time: $(($SECONDS / 3600)) hr : $(($(($SECONDS % 3600)) / 60)) min : $(($SECONDS % 60)) sec"
  echo " "
fi
#################################################

#### ASV DETECTION, TAXONOMY ASSIGNMENT, AND ####
######## PHYLOGENETIC TREE CONSTRUCTION #########
if echo $SKIP | grep -q "taxonomic_profiling"; then
  echo " "
  echo "*** Skipping ASV detection / taxonomy classification ***"
  echo " "
  
  #Create directory
  if [ -d "${RESULTS_DIR}/5.Amplicon_Sequence_Variants" ]
  then
    :
  else
    mkdir ${RESULTS_DIR}/5.Amplicon_Sequence_Variants
    mkdir ${RESULTS_DIR}/5.Amplicon_Sequence_Variants/0.ErrorOut
    mkdir ${RESULTS_DIR}/5.Amplicon_Sequence_Variants/0.Output
  fi
else
  SECONDS=0
  echo " "
  echo "*** Performing ASV detection / taxonomy classification ***"
  echo " "
  
  #Create directory for output
  if [ -d "${RESULTS_DIR}/5.Amplicon_Sequence_Variants" ]
  then
    :
  else
    mkdir ${RESULTS_DIR}/5.Amplicon_Sequence_Variants
    mkdir ${RESULTS_DIR}/5.Amplicon_Sequence_Variants/0.ErrorOut
    mkdir ${RESULTS_DIR}/5.Amplicon_Sequence_Variants/0.Output
  fi
  
  ##### Run main DADA2 workflow #####
  #Create script for running program and submit
  echo '#!/bin/bash' > run.sh
  echo "#SBATCH --partition=$PARTITION_5" >> run.sh
  echo "#SBATCH --job-name=DADA2" >> run.sh
  echo "#SBATCH --error=${RESULTS_DIR}/5.Amplicon_Sequence_Variants/0.ErrorOut/DADA2.err" >> run.sh
  echo "#SBATCH --output=${RESULTS_DIR}/5.Amplicon_Sequence_Variants/0.Output/DADA2.out" >> run.sh
  echo "#SBATCH --time=$TIME_5" >> run.sh
  echo "#SBATCH --ntasks=1" >> run.sh
  echo "#SBATCH --cpus-per-task=$CPU_5" >> run.sh
  echo "#SBATCH --mem-per-cpu=$MEM_5" >> run.sh
  echo "#SBATCH --mail-type=FAIL" >> run.sh
  echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> run.sh
  echo "#SBATCH --wait" >> run.sh
  echo "$PROG_LOAD" >> run.sh
  echo "echo \"suppressMessages(library(dada2))\" > run.R" >> run.sh
  echo "echo \"suppressMessages(library(DECIPHER))\" >> run.R" >> run.sh
  echo "echo \"suppressMessages(library(phangorn))\" >> run.R" >> run.sh
  echo "echo \"suppressMessages(library(phyloseq))\" >> run.R" >> run.sh
  echo "echo \"cat('\n','Package versions:', '\n')\" >> run.R" >> run.sh
  echo "echo \"cat('dada2 -', paste(packageVersion('dada2')), '\n')\" >> run.R" >> run.sh
  echo "echo \"cat('DECIPHER -', paste(packageVersion('DECIPHER')), '\n')\" >> run.R" >> run.sh
  echo "echo \"cat('phangorn -', paste(packageVersion('phangorn')), '\n')\" >> run.R" >> run.sh
  echo "echo \"cat('phyloseq -', paste(packageVersion('phyloseq')), '\n')\" >> run.R" >> run.sh
  echo "echo \"getN <- function(x) sum(getUniques(x))\" >> run.R" >> run.sh
  echo "echo \"filtpath <- '${RESULTS_DIR}/3.Quality_Controlled_Sequences'\" >> run.R" >> run.sh
  echo "echo \"filtFs <- sort(list.files(filtpath, pattern='_R1_001.fastq.gz', full.names=TRUE))\" >> run.R" >> run.sh
  echo "echo \"filtRs <- sort(list.files(filtpath, pattern='_R2_001.fastq.gz', full.names=TRUE))\" >> run.R" >> run.sh
  echo "echo 'sample.names <- sapply(strsplit(basename(filtFs), \"_\"), \`[\`, 1)' >> run.R" >> run.sh
  echo "echo 'sample.namesR <- sapply(strsplit(basename(filtRs), \"_\"), \`[\`, 1)' >> run.R" >> run.sh
  echo "echo \"names(filtFs) <- sample.names\" >> run.R" >> run.sh
  echo "echo \"names(filtRs) <- sample.namesR\" >> run.R" >> run.sh
  echo "echo \"cat('\n', 'Constructing the error model for forward reads:', '\n')\" >> run.R" >> run.sh
  echo "echo \"set.seed(1234)\" >> run.R" >> run.sh
  echo "echo \"errF <- learnErrors(filtFs, nbases=Inf, multithread=${CPU_5})\" >> run.R" >> run.sh
  echo "echo \"cat('\n', 'Constructing the error model for reverse reads:', '\n')\" >> run.R" >> run.sh
  echo "echo \"set.seed(1234)\" >> run.R" >> run.sh
  echo "echo \"errR <- learnErrors(filtRs, nbases=Inf, multithread=${CPU_5})\" >> run.R" >> run.sh
  echo "echo \"cat('\n', 'Plotting error model for forward reads.', '\n')\" >> run.R" >> run.sh
  echo "echo \"pdf('${RESULTS_DIR}/5.Amplicon_Sequence_Variants/ErrorRatePlotF.pdf')\" >> run.R" >> run.sh
  echo "echo \"plotErrors(errF, nominalQ=TRUE)\" >> run.R" >> run.sh
  echo "echo \"trash <- dev.off()\" >> run.R" >> run.sh
  echo "echo \"cat('\n', 'Plotting error model for reverse reads.', '\n')\" >> run.R" >> run.sh
  echo "echo \"pdf('${RESULTS_DIR}/5.Amplicon_Sequence_Variants/ErrorRatePlotR.pdf')\" >> run.R" >> run.sh
  echo "echo \"plotErrors(errR, nominalQ=TRUE)\" >> run.R" >> run.sh
  echo "echo \"trash <- dev.off()\" >> run.R" >> run.sh
  echo "echo \"mergers <- vector('list', length(sample.names))\" >> run.R" >> run.sh
  echo "echo \"names(mergers) <- sample.names\" >> run.R" >> run.sh
  echo "echo \"cat('\n', 'Performing main DADA2 algorithm:', '\n')\" >> run.R" >> run.sh
  echo "echo \"for(sample in sample.names) {\" >> run.R" >> run.sh
  echo "echo \"  cat('\n', 'Processing:', sample, '\n')\" >> run.R" >> run.sh
  echo "echo \"  derepF <- derepFastq(filtFs[[sample]])\" >> run.R" >> run.sh
  echo "echo \"  ddF <- dada(derepF, err=errF, multithread=${CPU_5})\" >> run.R" >> run.sh
  echo "echo \"  cat('Forward reads remaining after denoising:',getN(ddF),'\n')\" >> run.R" >> run.sh
  echo "echo \"  derepR <- derepFastq(filtRs[[sample]])\" >> run.R" >> run.sh
  echo "echo \"  ddR <- dada(derepR, err=errR, multithread=${CPU_5})\" >> run.R" >> run.sh
  echo "echo \"  cat('Reverse reads remaining after denoising:',getN(ddR),'\n')\" >> run.R" >> run.sh
  echo "echo \"  merger <- mergePairs(ddF, derepF, ddR, derepR)\" >> run.R" >> run.sh
  echo "echo \"  cat('Reads remaining after merging:',getN(merger),'\n')\" >> run.R" >> run.sh
  echo "echo \"  mergers[[sample]] <- merger\" >> run.R" >> run.sh
  echo "echo \"}\" >> run.R" >> run.sh
  echo "echo \"rm(derepF); rm(derepR)\" >> run.R" >> run.sh
  echo "echo \"seqtab <- makeSequenceTable(mergers)\" >> run.R" >> run.sh
  echo "echo \"cat('\n','ASV sequence length distribution:', '\n')\" >> run.R" >> run.sh
  echo "echo \"table(nchar(getSequences(seqtab)))\" >> run.R" >> run.sh
  echo "echo \"lower_bound=round(mean(nchar(getSequences(seqtab)))-sd(nchar(getSequences(seqtab))), 0)\" >> run.R" >> run.sh
  echo "echo \"upper_bound=round(mean(nchar(getSequences(seqtab)))+sd(nchar(getSequences(seqtab))), 0)\" >> run.R" >> run.sh
  echo "echo \"cat('\n','Extracting ASVs with sequence length',lower_bound,'-',upper_bound, '\n')\" >> run.R" >> run.sh
  echo "echo \"seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(lower_bound,upper_bound)]\" >> run.R" >> run.sh
  echo "echo \"seqtab.nochim <- removeBimeraDenovo(seqtab2, method='consensus', multithread=${CPU_5})\" >> run.R" >> run.sh
  echo "echo \"cat('\n','Fraction of reads remaining after chimera removal:', sum(seqtab.nochim)/sum(seqtab2), '\n')\" >> run.R" >> run.sh
  echo "echo \"seqtab.df <- data.frame(SampleID=rownames(data.frame(seqtab.nochim)), data.frame(seqtab.nochim))\" >> run.R" >> run.sh
  echo "echo \"write.table(seqtab.df, '${RESULTS_DIR}/5.Amplicon_Sequence_Variants/ASV_table.txt', sep='\t', row.names=F, quote=F)\" >> run.R" >> run.sh
  echo "echo \"cat('\n','ASV detection complete. Assigning taxonomy...', '\n')\" >> run.R" >> run.sh
  echo "echo \"taxa <- assignTaxonomy(seqtab.nochim, '${TAXA_REF}', minBoot=80, tryRC=T, multithread=${CPU_5})\" >> run.R" >> run.sh
  echo "echo \"cat('\n','Genus level taxonomy assigned.', '\n')\" >> run.R" >> run.sh
  if [[ ! -z "$SPEC_REF" ]]; then
    echo "echo \"cat('\n','Adding species designations...', '\n')\" >> run.R" >> run.sh
    echo "echo \"taxa <- addSpecies(taxa, '${SPEC_REF}', allowMultiple = TRUE)\" >> run.R" >> run.sh
    echo "echo \"cat('\n','Species level taxonomy assigned.', '\n')\" >> run.R" >> run.sh
  else
    echo "echo \"cat('\n','No species reference file provided. Species designations will not be added to taxonomic assignments.', '\n')\" >> run.R" >> run.sh
  fi
  echo "echo \"taxa.df <- data.frame(ASV=rownames(data.frame(taxa)), data.frame(taxa))\" >> run.R" >> run.sh
  echo "echo \"write.table(taxa.df, '${RESULTS_DIR}/5.Amplicon_Sequence_Variants/ASV_taxa_IDs.txt', sep='\t', row.names=F, quote=F)\" >> run.R" >> run.sh
  if echo $SKIP | grep -q "phylogenetic_tree"; then
    echo "echo \"cat('\n','Skipping phylogenetic tree construction.', '\n')\" >> run.R" >> run.sh
    echo "echo 'phyloseq.object <- phyloseq(tax_table(taxa), otu_table(seqtab.nochim, taxa_are_rows=F))\" >> run.R" >> run.sh
    echo "echo \"cat('\n','ASV abundances and taxonomic classifications have been combined into a phyloseq object.', '\n')\" >> run.R" >> run.sh
  else
    echo "echo \"seqs <- colnames(seqtab.nochim)\" >> run.R" >> run.sh
    echo "echo \"names(seqs) <- seqs\" >> run.R" >> run.sh
    echo "echo \"cat('\n','Performing multiple sequence alignment with DECIPHER...', '\n')\" >> run.R" >> run.sh
    echo "echo \"alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA, processors=${CPU_5}, verbose=T)\" >> run.R" >> run.sh
    echo "echo \"writeXStringSet(alignment, file='${RESULTS_DIR}/5.Amplicon_Sequence_Variants/ASV_mult_seq_align.fa')\" >> run.R" >> run.sh
    echo "echo \"cat('\n','Multiple sequence alignment completed.', '\n')\" >> run.R" >> run.sh
    echo "echo \"phang.align <- as.phyDat(as(alignment, 'matrix'), type='DNA')\" >> run.R" >> run.sh
    echo "echo \"dm <- dist.ml(phang.align)\" >> run.R" >> run.sh
    echo "echo \"cat('\n','Distance matrix of alignment computed.','\n')\" >> run.R" >> run.sh
    echo "echo \"treeNJ <- NJ(dm)\" >> run.R" >> run.sh
    echo "echo 'treeNJ\$edge.length[which(treeNJ\$edge.length < 0)] <- 0' >> run.R" >> run.sh
    echo "echo \"cat('\n','Initial phylogenetic tree constructed using Neighbor-Joining. This will be used as starting tree for phylogenetic tree construction.', '\n')\" >> run.R" >> run.sh
    if [[ ! -z RAXML ]]; then
      echo "echo \"cat('\n','Constructing final phylogenetic tree with RAxML-NG using GTR+G+I model...','\n')\" >> run.R" >> run.sh
      echo "echo 'write.tree(treeNJ,\"${RESULTS_DIR}/5.Amplicon_Sequence_Variants/temp.start.tre\")' >> run.R" >> run.sh
      echo "echo \"system('raxml-ng --msa ${RESULTS_DIR}/5.Amplicon_Sequence_Variants/ASV_mult_seq_align.fa --tree ${RESULTS_DIR}/5.Amplicon_Sequence_Variants/temp.start.tre --model GTR+G+I --seed 2 --threads auto{${CPU_5}} --workers auto{${CPU_5}} --prefix ${RESULTS_DIR}/5.Amplicon_Sequence_Variants/temp')\" >> run.R" >> run.sh
      echo "echo \"system('mv ${RESULTS_DIR}/5.Amplicon_Sequence_Variants/temp.raxml.bestTree ${RESULTS_DIR}/5.Amplicon_Sequence_Variants/ASV_phylo_tree.tre')\" >> run.R" >> run.sh
      echo "echo \"system('rm ${RESULTS_DIR}/5.Amplicon_Sequence_Variants/temp*')\" >> run.R" >> run.sh
      echo "echo \"cat('Phylogenetic tree constructed.','\n')\" >> run.R" >> run.sh
      echo "echo \"phy.tree <- ape::read.tree('${RESULTS_DIR}/5.Amplicon_Sequence_Variants/ASV_phylo_tree.tre')\" >> run.R" >> run.sh
      echo "echo \"phyloseq.object <- phyloseq(tax_table(taxa), phy_tree(phy.tree),\" >> run.R" >> run.sh
      echo "echo \"                            otu_table(seqtab.nochim, taxa_are_rows=F))\" >> run.R" >> run.sh
      echo "echo \"cat('\n','ASV abundances, taxonomic classifications, and phylogenetic tree have been combined into a phyloseq object.', '\n')\" >> run.R" >> run.sh
    else
      echo "echo \"cat('\n','Constructing final phylogenetic tree with phangorn using GTR+G+I model...','\n')\" >> run.R" >> run.sh
      echo "echo \"fit <- pml(treeNJ, data=phang.align)\" >> run.R" >> run.sh
      echo "echo \"fitGTR <- update(fit, k=4, inv=0.2)\" >> run.R" >> run.sh
      echo "echo \"fitGTR <- optim.pml(fitGTR, model='GTR', optInv=TRUE, optGamma=TRUE,\" >> run.R" >> run.sh
      echo "echo \"                    rearrangement = 'stochastic', control = pml.control(trace=0))\" >> run.R" >> run.sh
      echo "echo \"cat('\n','Phylogenetic tree constructed.','\n')\" >> run.R" >> run.sh
      echo "echo 'write.tree(fitGTR\$tree,\"${RESULTS_DIR}/5.Amplicon_Sequence_Variants/ASV_phylo_tree.tre\")' >> run.R" >> run.sh
      echo "echo 'phyloseq.object <- phyloseq(tax_table(taxa), phy_tree(fitGTR\$tree),' >> run.R" >> run.sh
      echo "echo \"                           otu_table(seqtab.nochim, taxa_are_rows=F))\" >> run.R" >> run.sh
      echo "echo \"cat('\n','ASV abundances, taxonomic classifications, and phylogenetic tree have been combined into a phyloseq object.', '\n')\" >> run.R" >> run.sh
    fi
  fi
  echo "echo \"cat('\n', 'Phyloseq object summary:', '\n')\" >> run.R" >> run.sh
  echo "echo \"phyloseq.object\" >> run.R" >> run.sh
  echo "echo \"saveRDS(phyloseq.object, '${RESULTS_DIR}/5.Amplicon_Sequence_Variants/ASV_phyloseq.rds')\" >> run.R" >> run.sh
  echo "Rscript run.R > ${RESULTS_DIR}/5.Amplicon_Sequence_Variants/0.ASV_detection_classification.log 2>&1" >> run.sh
  echo "rm run.R" >> run.sh
  chmod +x run.sh
  
  sbatch run.sh > /dev/null
  
  rm run.sh
  
  #Signal jobs have ended
  echo "ASV detection / taxonomic classification complete"
  echo "Elapsed time: $(($SECONDS / 3600)) hr : $(($(($SECONDS % 3600)) / 60)) min : $(($SECONDS % 60)) sec"
  echo " "
fi
#################################################

###### REMOVE INTERMEDIATE SEQUENCE FILES #######

if [[ ! -z $DELETE ]]; then
  echo "*** Removing intermediate sequence files and reorganizing ***"
  echo " "
  
  #Create new directories needed
  mkdir ${RESULTS_DIR}/2.Processed_Sequences
  mkdir ${RESULTS_DIR}/2.Processed_Sequences/Log_Files

  #Create combined log file for each sample
  for file in ${SEQ_DIR}/*R1_001.${SEQ_EXT}; do
    FILE_NAME=$(echo $file | awk -F '/' '{print $NF}' | awk -F '_R1_001' '{print $1}')

    #Initialize log file
    echo "*** Log file for sample ${FILE_NAME} ***" > ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log
    echo " " >> ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log

    #Grab log file contents
    echo "### Log for primer removal using cutadapt" >> ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log
    echo " " >> ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log
    cat ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log \
        ${RESULTS_DIR}/2.Primer_Trimmed_Sequences/${FILE_NAME}.log \
        > temp
    mv temp ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log
    echo " " >> ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log

    echo "### Log for quality trimming/filtering using DADA2's filterAndTrim" >> ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log
    echo " " >> ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log
    cat ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log \
        ${RESULTS_DIR}/3.Quality_Controlled_Sequences/0.filterAndTrim.log \
        > temp
    mv temp ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log
    echo " " >> ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log
  done

  #Grab processed sequences
  cp ${RESULTS_DIR}/3.Quality_Controlled_Sequences/*fastq.gz \
     ${RESULTS_DIR}/2.Processed_Sequences/
  ERROR=$(diff <(ls ${RESULTS_DIR}/2.Processed_Sequences/*fastq.gz | awk -F'/' '{print $NF}') \
               <(ls ${RESULTS_DIR}/3.Quality_Controlled_Sequences/*fastq.gz | awk -F'/' '{print $NF}'))
  if [[ ! -z "$ERROR" ]]; then
    echo "ERROR: Something went wrong when copying processed sequences to new directory, copied sequences do not match sequences original sequences"
    exit 1
  fi

  #Modify numbering of post sequence processing directories
  mv ${RESULTS_DIR}/4.FastQC_Final_Reports ${RESULTS_DIR}/3.FastQC_Final_Reports
  mv ${RESULTS_DIR}/5.Amplicon_Sequence_Variants ${RESULTS_DIR}/4.Amplicon_Sequence_Variants
  
  #Remove directories with intermediate sequences
  rm -rf ${RESULTS_DIR}/2.Primer_Trimmed_Sequences
  rm -rf ${RESULTS_DIR}/3.Quality_Controlled_Sequences
fi
#################################################

echo "*** 16S amplicon pipeline complete ***"
