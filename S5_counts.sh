#!/bin/bash
#SBATCH --time=14:00:00
#SBATCH --mem=8gb
#SBATCH --cpus-per-task=4
#SBATCH --job-name=counts
#SBATCH --output=array_%J.out
#SBATCH --error=array_%J.err
#SBATCH --partition=pibu_el8

# Store important directories and filepaths as variables
WORKDIR="/data/users/asteiner/RNA-seq_paper/3_Map_reads_to_the_reference_genome__and__4_Count_the_number_of_reads_per_gene"
INDIR="$WORKDIR/mapped_reads"
INFILE="$WORKDIR/genome.gtf"
OUTFILE="$WORKDIR/counts.txt"

# Using apptainer i can use the container that stores the command featureCounts.
# This command is used to generate read counts for each gene and output it in a table which can then be used for further analysis.
# Using the option -p i specify that the data provided is from paired end reads. -T specifies to use 12 threads to make the processing
# faster. With the option -t exon i specify that i would like to count the reads that are overlapping with an exon, using -g gene_id
# i can specify how the counts should be grouped in the table (here by using the gene ids form my gtf file). With the option -a i define
# the annotation file to be used, -o specifies the name and path of the output file containing the readcounts. And the last argument 
# provided refers to the bam files that should be used to calculate the counts.
apptainer exec --bind /data/ /containers/apptainer/subread_2.0.1--hed695b0_0.sif featureCounts -p -T 12 -t exon -g gene_id -a $INFILE -o $OUTFILE $INDIR/*_mapped_sorted.bam
