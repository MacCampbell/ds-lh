#!/bin/bash -l
# NOTE the -l flag!

# This script is run some test associations on ds lh data.

# Name of the job 
#SBATCH -J doGeno

#Email myself
#SBATCH --mail-user=maccampbell@ucdavis.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

#Job level

#SBATCH --partition=high
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=8
#SBATCH --ntasks=24
#SBATCH --time=1-01:00:00 #run for a day and an hour

# Standard out and Standard Error output files with the job number in the name.
#SBATCH -o bench-%j.output
#SBATCH -e bench-%j.output

# hostname is just for debugging
hostname
export OMP_NUM_THREADS=$SLURM_NTASKS
module load benchmarks

# The main job executable to run: note the use of srun before it
# Using $HOME for better practice

srun $HOME/angsd/angsd -P 24 -b $HOME/ds-lh/bamlists/216.bamlist -minInd 194 -out $HOME/ds-lh/outputs/200/216-doGeno \
-minMaf 0.05 -minMapQ 30 -minQ 20 -GL 1 -doMajorMinor 1 -doMaf 1 \
-ref $HOME/ds-lh/genome/Hypomesus-transpacificus_10X_F_A.pseudohap2.1.fasta \
-SNP_pval 1e-6 -doGeno 2 -doPost 1 -postCutoff 0.95 -setMinDepthInd 5

srun $HOME/angsd/angsd -P 24 -b $HOME/ds-lh/bamlists/216.bamlist -minInd 194 -out $HOME/ds-lh/outputs/200/216-doPlink \
-minMaf 0.05 -minMapQ 30 -minQ 20 -GL 1 -doMajorMinor 1 -doMaf 1 \
-ref $HOME/ds-lh/genome/Hypomesus-transpacificus_10X_F_A.pseudohap2.1.fasta 
-SNP_pval 1e-6 -doGeno 2 -doPost 1 -postCutoff 0.95 -setMinDepthInd 5 -doPlink 2

