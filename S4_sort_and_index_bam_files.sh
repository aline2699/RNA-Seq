#!/bin/bash
#SBATCH --array=1-12
#SBATCH --time=20:00:00
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=3
#SBATCH --job-name=sort_and_index_bam_files
#SBATCH --output=array_%J.out
#SBATCH --error=array_%J.err
#SBATCH --partition=pibu_el8

# Specify the paths to the working directory and the outputdirectory
WORKDIR="/data/users/asteiner/RNA-seq_paper/3_Map_reads_to_the_reference_genome__and__4_Count_the_number_of_reads_per_gene"
OUTDIR="$WORKDIR/mapped_reads"

# Store the path to the samplelistfile in the variable SAMPLELIST
SAMPLELIST="/data/users/asteiner/RNA-seq_paper/2_Quality_checks/samplelist.txt"
# Use awk on the samplelist file to extract the names of each of my samples 
SAMPLE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`

# Go into the output directory
cd $OUTDIR

# I use apptainer to access the container containing the hisat2_samtools
# To be able to only call on the container once and use bash syntax i use bash -c "..." to pass all the commands in one line
# I use samtools sort to sort my bam files by the genomic coordinates of the reads. 
# The option -@ 12 specifies the use of 12 CPU threads to speed up processing. Then the path to the previously created bam files
# is provided. Then the option -o specifies the name and path of the output files.
# The second command samtools index is used to create the bai files too, this allows for quick access to specific regions of the BAM
# files by creating an index for the sorted bam files. The argument provided to this command is the path and name of the input files
# / mapped bam files created in the step befor.
apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif bash -c "samtools sort -@ 12 $OUTDIR/${SAMPLE}_mapped.bam -o $OUTDIR/${SAMPLE}_mapped_sorted.bam ; samtools index $OUTDIR/${SAMPLE}_mapped_sorted.bam"