#!/bin/bash
#SBATCH --array=1-12
#SBATCH --time=02:00:00
#SBATCH --mem=2g
#SBATCH --cpus-per-task=2
#SBATCH --job-name=slurm_array
#SBATCH --output=array_%J.out
#SBATCH --error=array_%J.err
#SBATCH --partition=pibu_el8

# Here I define important paths that I will be using and store them in variables.
WORKDIR="/data/users/asteiner/RNA-seq_paper/2_Quality_checks"
OUTDIR="$WORKDIR/QC_results"
SAMPLELIST="$WORKDIR/samplelist.txt"

# Then i also get the paths to the Reads i will use and to the samplenames and store them in variables.
SAMPLE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
READ1=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST`
READ2=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $3; exit}' $SAMPLELIST`

# Here i specify where the quality control output should go.
OUTFILE="$OUTDIR/${SAMPLE}"

############################


# First i make sure the directories where i want to safe my results in exist. The -p command makes sure that all the
# needed parent directories will also be made in this step.
mkdir -p $OUTFILE

# Then it is important to load all the modules that will be used in the code later on.
module load FastQC/0.11.9-Java-11

# The fastqc command does the Quality Control. The first two arguments are the inputfiles it will use, the -o gives us the
# possibility to define the output file, and the argument after -o is the path and name to and of that outputfile.
fastqc $READ1 $READ2 -o $OUTFILE

#