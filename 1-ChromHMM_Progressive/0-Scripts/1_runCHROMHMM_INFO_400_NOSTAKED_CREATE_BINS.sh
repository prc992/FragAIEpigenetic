#!/bin/bash
#SBATCH --job-name=1_runCHROMHMM_INFO_400_NOSTAKED_CREATE_BINS
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=2-00:00:00  # Tempo limite: 2 dias

# Carregar módulo Java
module load java/14.0.1

MODELNAME="0-Bins_400bp"

# Caminhos input
CHROMHMM_JAR="/PHShome/pd004/FragAIEpigenetic/1-ChromHMM_Progressive/1-program/ChromHMM.jar"
CHROM_SIZES="/PHShome/pd004/FragAIEpigenetic/1-ChromHMM_Progressive/1-program/CHROMSIZES/hg19.txt"
RAW_BAMS="/PHShome/pd004/FragAIEpigenetic/1-ChromHMM_Progressive/2-raw-files-sl"
CELLMARK="/PHShome/pd004/FragAIEpigenetic/1-ChromHMM_Progressive/3-ChromHMM-Suport-Files/cellmarkfiletable_no_stacked.txt"

#Caminho output
BINARIZED="/PHShome/pd004/FragAIEpigenetic/1-ChromHMM_Progressive/3-ChromHMM-Suport-Files/0-BINS/$MODELNAME"
mkdir -p "$BINARIZED"


# Mostrar o comando como uma linha única
echo "🔹 Comando a ser executado:"
echo "java -jar \"$CHROMHMM_JAR\" BinarizeBam -b 400 \"$CHROM_SIZES\" \"$RAW_BAMS\" \"$CELLMARK\" \"$BINARIZED\" "

# Executar o comando como linha única
java -jar "$CHROMHMM_JAR" BinarizeBam -b 400 "$CHROM_SIZES" "$RAW_BAMS" "$CELLMARK" "$BINARIZED" 