#!/bin/bash
#SBATCH --job-name=2_runCHROMHMM_INFO_400_NOSTAKED_LEARNMODEL
#SBATCH --array=2-8
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=ALL
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --time=03:00:00

set -euo pipefail

# Carregar Java
module load java/14.0.1

# ===== Parâmetros =====
WINDOW_SIZE=400
MODELNAMEBIN="0-Bins_400bp"   # nome do diretório de bins já gerados

# Caminhos input
CHROMHMM_JAR="/PHShome/pd004/FragAIEpigenetic/1-ChromHMM_Progressive/1-program/ChromHMM.jar"
BINARIZED="/PHShome/pd004/FragAIEpigenetic/1-ChromHMM_Progressive/3-ChromHMM-Suport-Files/0-BINS/${MODELNAMEBIN}"

# Base do output
MODELS_BASE="/PHShome/pd004/FragAIEpigenetic/1-ChromHMM_Progressive/3-ChromHMM-Suport-Files/1-MODELS"

# Threads
N_THREADS="${SLURM_CPUS_PER_TASK:-1}"

# Valor de K vem do índice do array
K="${SLURM_ARRAY_TASK_ID}"

# Pastas de saída por K
MODELNAMELEARN="0-MODEL_${K}_STATES"
MODEL_OUT="${MODELS_BASE}/${MODELNAMELEARN}"
mkdir -p "${MODEL_OUT}"

# Sanity checks
[[ -f "${CHROMHMM_JAR}" ]] || { echo "ERRO: ChromHMM.jar não encontrado em ${CHROMHMM_JAR}"; exit 1; }
[[ -d "${BINARIZED}" ]] || { echo "ERRO: diretório de bins não encontrado em ${BINARIZED}"; exit 1; }

echo "============================================================="
echo "ChromHMM LearnModel | K=${K} estados"
echo "Bins:     ${BINARIZED}"
echo "Saída:    ${MODEL_OUT}"
echo "Jar:      ${CHROMHMM_JAR}"
echo "Threads:  ${N_THREADS}"
echo "Janela:   ${WINDOW_SIZE} bp"
echo "============================================================="
echo "Comando:"
echo "java -jar \"${CHROMHMM_JAR}\" LearnModel -b ${WINDOW_SIZE} -p ${N_THREADS} \"${BINARIZED}\" \"${MODEL_OUT}\" ${K} hg19"
echo "============================================================="

# Execução (opcional: -noautoopen para evitar abrir browser)
time java -jar "${CHROMHMM_JAR}" LearnModel \
    -b "${WINDOW_SIZE}" \
    -p "${N_THREADS}" \
    -noautoopen \
    "${BINARIZED}" \
    "${MODEL_OUT}" \
    "${K}" \
    hg19

echo "✅ Finalizado K=${K}: resultados em ${MODEL_OUT}"