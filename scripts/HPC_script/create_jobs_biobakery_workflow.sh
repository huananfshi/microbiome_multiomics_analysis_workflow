#!/bin/bash

IFS=$'\n'

# step 1 - get sample ID
# step 2 - make sample id output directory
# step 3 - run star using IRFINDER reference genome


### paired fastq files in one subfolder:
## RF13987_S29_L001_R1_001.fastq.gz
## RF13987_S29_L001_R2_001.fastq.gz
## RF13987_L001-ds.51ebe69d5ad44a0eb1b9a5872ecbc317


### sample directory structure:
### |-- RF15788_L001-ds.bd32df1ae0b84eac80e26660aea121fd
### |   |-- RF15788_S2_L001_R1_001.fastq.gz
### |   `-- RF15788_S2_L001_R2_001.fastq.gz
### |-- RF15793_L001-ds.5e1392310f5345bcb4370f8852a46b01
### |   |-- RF15793_S7_L001_R1_001.fastq.gz
### |   `-- RF15793_S7_L001_R2_001.fastq.gz

### for each directory:
### 	step 1 - loop through each sample directory
###	step 2 - get sample filename prefix (sample ID + variable string)
###	step 2 - get sample ID
###	step 3 - sed sample directory, sample file & sample ID to template


### don't need to move or unzip... STAR can process .fastq.gz files natively
### using the "--readFilesCommand zcat" parameter


JOB_TEMPLATE=biobakery_workflow_TEMPLATE.pbs

### paired-end samples (R1 & R2) are saved in a subdirectory of the WORK_DIRECTORY
### the WORK_DIRECTORY will thus become PBS_O_WORKDIR for the PBS job

WORK_DIRECTORY=/mount/durgan/EODF/
#WORK_DIRECTORY=/mount/schiff/FASTQ_Generation_2019-05-31_17_55_34Z-185446111/Bt474mutationA
cd ${WORK_DIRECTORY}

### all sample subdirectories end with sequence
for sample_dir in `ls Durgan* -d`
do
	### print the directory
	echo ${sample_dir}

	#cd ${sample_dir}
	#sample_file=`ls *R1*.fastq.gz | cut -d'_' -f1-2`
	#sample_id=`echo ${sample_dir} | cut -d'_' -f1`
	#cd ..

	sed -e "s/SAMPLE_DIRECTORY/${sample_dir}/g" /mount/durgan/biobakery_workflow_TEMPLATE.pbs >> /mount/durgan/EODF/biobakery_workflow_${sample_dir}.pbs

	### submit the job
	#qsub ./Bt474mutationA_${sample_id}.pbs

done



