#!/bin/bash
#SBATCH --job-name=multi_matrix_krivonogova
#SBATCH --output=benchmark_cluster.csv
#SBATCH --ntasks=1
#
#SBATCH --cpus-per-task=20
#
#SBATCH --time=55:00
#SBATCH --mem-per-cpu=100


./main
