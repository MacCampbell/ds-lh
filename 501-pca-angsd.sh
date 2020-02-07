#!/bin/bash -l
# NOTE the -l flag!

# This script is run some test associations on ds lh data.

# Name of the job 
#SBATCH -J 202

#Email myself
#SBATCH --mail-user=maccampbell@ucdavis.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

#Job level

#SBATCH --partition=high
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=5
#SBATCH --ntasks=10
#SBATCH --time=1-08:00:00 #run for a day and an eight hours

# Standard out and Standard Error output files with the job number in the name.
#SBATCH -o bench-%j.output
#SBATCH -e bench-%j.output

# hostname is just for debugging
hostname
export OMP_NUM_THREADS=$SLURM_NTASKS
module load benchmarks

#Doing the calculations making sure sites are present in 90% of individuals.
srun $HOME/angsd/angsd -minInd 194 -GL 1 -out $HOME/ds-lh/outputs/500/216-pca -nThreads 10 -doGlf 2 -doMajorMinor 1 -doMaf 2 -SNP_pval 1e-6 -minMapQ 30 -minQ 20 \
-bam $HOME/ds-lh/bamlists/216.bamlist

# Generate a covariance matrix
# For Linux
srun python $HOME/pcangsd/pcangsd.py -beagle $HOME/ds-lh/outputs/500/216-pca.beagle.gz -admix -o $HOME/ds-lh/outputs/500/216-pca

