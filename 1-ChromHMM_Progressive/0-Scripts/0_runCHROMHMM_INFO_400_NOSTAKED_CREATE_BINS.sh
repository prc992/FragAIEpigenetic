#!/bin/bash
#SBATCH --job-name=0_runCHROMHMM_INFO_400_NOSTAKED_CREATE_BINS
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=ALL
#SBATCH --output=logs/0_runCHROMHMM_INFO_400_NOSTAKED_CREATE_BINS-%j.out
#SBATCH --error=logs/0_runCHROMHMM_INFO_400_NOSTAKED_CREATE_BINS-%j.err
#SBATCH --time=2-00:00:00  # Tempo limite: 2 dias

# Carregar módulo Java
module load java/14.0.1

MODELNAME="0-MODEL_2_STATES"

# Caminhos input
CHROMHMM_JAR="/data/baca/projects/ultima_deep/4-ChromHMM/0-program/ChromHMM.jar"
CHROM_SIZES="/data/baca/projects/ultima_deep/4-ChromHMM/0-program/CHROMSIZES/hg19.txt"
RAW_BAMS="/data/baca/projects/ultima_deep/4-ChromHMM/2-raw_files"
CELLMARK="/data/baca/projects/ultima_deep/4-ChromHMM/0_1-support_files/cellmarkfiletable_no_stacked.txt"

#Caminho output
BINARIZED="/data/baca/projects/ultima_deep/4-ChromHMM_Progressive/3-output/$MODELNAME/bins"

# Mostrar o comando como uma linha única
echo "🔹 Comando a ser executado:"
echo "java -jar \"$CHROMHMM_JAR\" BinarizeBam -b 400 \"$CHROM_SIZES\" \"$RAW_BAMS\" \"$CELLMARK\" \"$BINARIZED\" "

# Executar o comando como linha única
java -jar "$CHROMHMM_JAR" BinarizeBam -b 400 "$CHROM_SIZES" "$RAW_BAMS" "$CELLMARK" "$BINARIZED" 