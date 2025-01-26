#!/bin/bash
#SBATCH --array=1-12
#SBATCH --time=20:00:00
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=3
#SBATCH --job-name=mapped_bam_files
#SBATCH --output=array_%J.out
#SBATCH --error=array_%J.err
#SBATCH --partition=pibu_el8


# Save the path to the reads in the INDIR variable
INDIR="/data/courses/rnaseq_course/breastcancer_de/reads"
# Define the path to the directory i want to work in and store it as WORKDIR
WORKDIR="/data/users/asteiner/RNA-seq_paper/3_Map_reads_to_the_reference_genome__and__4_Count_the_number_of_reads_per_gene"
# Define in which directory the output should be saved and name it OUTDIR
OUTDIR="$WORKDIR/mapped_reads"

# With the variable INDEXFILES i store the path to the index files produced in the previous step.
INDEXFILES="$WORKDIR/genome_index/genome_index"

# Using the variable SAMPLELIST i can reference the paths to my individual reads as well as their names stored in a text document.
SAMPLELIST="/data/users/asteiner/RNA-seq_paper/2_Quality_checks/samplelist.txt"
# Using awk i can extract the name, and paths to read1 and read2 for my samples from the SAMPLELIST
SAMPLE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
READ1=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST`
READ2=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $3; exit}' $SAMPLELIST`

# To make sure my outputdirectory exists i use mkdir -p which creates the directory in case it doesn't exist
mkdir -p $OUTDIR

# I use apptainer to access the container containing the hisat2_samtools
# To be able to only call on the container once and use bash syntax i use bash -c "..." to pass all the commands in one line
# Using the command hisat2 i create the sam files. I use 12 CPU threads to speed up the alignment process. The option -x specifies
# which files should be used for indexing. The options -1 and -2 refer to the forward reads and reverse reads paths. With -S i
# store the output as samfiles in the provided path and file ($OUTDIR/$SAMPLE.sam). And lastly i redirect the standarderror to a 
# log file to get alignment statistics, warnings, and errors in a specified file.
# The second command used is samtools view to open the sam files created in the first step (option -S specifies the input files to
# be of sam format) and store them in the bam format (option -b makes the output to be saved as bam file). Then the actual path to 
# the input file is added and the file name after the assignment > will be the file the output gets stored in.
apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif bash -c "hisat2 -p 12 -x $INDEXFILES -1 $READ1 -2 $READ2 -S $OUTDIR/$SAMPLE.sam 2> $OUTDIR/${SAMPLE}_hisat2_summary.log ; samtools view -S -b $OUTDIR/$SAMPLE.sam > $OUTDIR/${SAMPLE}_mapped.bam"