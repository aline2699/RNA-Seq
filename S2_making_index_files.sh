#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --mem=8gb
#SBATCH --cpus-per-task=16
#SBATCH --job-name=index_files
#SBATCH --output=array_%J.out
#SBATCH --error=array_%J.err
#SBATCH --partition=pibu_el8

WORKDIR="/data/users/asteiner/RNA-seq_paper/3_Map_reads_to_the_reference_genome__and__4_Count_the_number_of_reads_per_gene"
OUTDIR="$WORKDIR/genome_index"


# -----------------------------------------
# The following commands to get the reference genome and the index file have been commented out, since these steps were performed
# manually in the terminal instead of using a script.

# get ref genome
# wget ftp://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# check sum
# sum Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# get annotation genome
# wget ftp://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz

# check sum
# sum Homo_sapiens.GRCh38.113.gtf.gz

# unzip and rename the annotation file and the reference genome file
# gzip -d Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
# mv Homo_sapiens.GRCh38.dna.primary_assembly.fa genome.fa
# gzip -d Homo_sapiens.GRCh38.113.gtf.gz
# mv Homo_sapiens.GRCh38.113.gtf genome.gtf

# -------------------------------------------


# create a directory, where the index files will get stored
mkdir -p $WORKDIR/genome_index

# use the hisat2-build command from hisat2_samtools to create the indexfiles needed for indexing the reads
apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif hisat2-build -p 16 genome.fa genome_index/genome_index
