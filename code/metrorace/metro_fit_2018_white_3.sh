#! /bin/bash
##############################################
#  SBATCH CONFIG
##############################################
## resources
#SBATCH -p Lewis
#SBATCH --qos long
#SBATCH -N 1   # nodes
#SBATCH -n 28  # cores
#SBATCH --mem-per-cpu=8G
#SBATCH -t 7-00:00
#
## labels and outputs
#SBATCH -A stsn
#SBATCH -J metro_fit_2018_white_3
#SBATCH -o metro_fit_2018_white_3-%j.out
#SBATCH -e metro_fit_2018_white_3.e
#
## notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=themattsimpson@gmail.com


echo "### Starting at: $(date) ###"

# load modules then display what we have
module load miniconda3
module list

source activate r_env

# Serial operations - only runs on the first core
echo "### Starting at: $(date)"
echo "First core reporting from node:"
hostname

echo "Currently working in directory:"
pwd

echo "Files in this folder:"
ls -l

# start R
R --no-restore --no-save CMD BATCH metro_fit_2018_white_3.R

echo "### Ending at: $(date) ###"

