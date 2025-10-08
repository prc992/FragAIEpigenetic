#!/bin/bash
#SBATCH --job-name=4_CreateBedWindow_400
#SBATCH --array=2-8
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --time=10:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL

set -euo pipefail

module load bedtools

# ======== PARÂMETROS GERAIS ========
MODEL_NUM="${SLURM_ARRAY_TASK_ID}"   # número de estados (2–8)
BASE_DIR="/PHShome/pd004/FragAIEpigenetic/1-ChromHMM_Progressive"
INPUT_DIR="${BASE_DIR}/3-ChromHMM-Suport-Files/1-MODELS/0-MODEL_${MODEL_NUM}_STATES"
OUTPUT_DIR="${INPUT_DIR}/1-WINDOWS-400BP"
BIN_SIZE=400

mkdir -p "${OUTPUT_DIR}"

[[ -d "${INPUT_DIR}" ]] || { echo "❌ Diretório de entrada não existe: ${INPUT_DIR}"; exit 1; }

echo "============================================================="
echo "🔹 Criando janelas de ${BIN_SIZE} bp para modelo ${MODEL_NUM} estados"
echo "Entrada: ${INPUT_DIR}"
echo "Saída:   ${OUTPUT_DIR}"
echo "============================================================="

# ======== LOCALIZAÇÃO DOS ARQUIVOS ========
mapfile -t BED_FILES < <(find "${INPUT_DIR}" -maxdepth 1 -type f -name "*_sorted.bed" | sort)
TOTAL_FILES=${#BED_FILES[@]}

if (( TOTAL_FILES == 0 )); then
    echo "⚠️ Nenhum arquivo *_sorted.bed encontrado em ${INPUT_DIR}"
    exit 0
fi

# ======== LOOP SOBRE ARQUIVOS ========
for INPUT_BED in "${BED_FILES[@]}"; do
    BASENAME=$(basename "${INPUT_BED}")
    OUTPUT_BED="${OUTPUT_DIR}/${BASENAME%.bed}_400bp.bed"

    echo "→ Criando janelas: ${BASENAME} -> $(basename "${OUTPUT_BED}")"
    # Mantido exatamente como você já usa:
    bedtools makewindows -b "${INPUT_BED}" -w "${BIN_SIZE}" -i src > "${OUTPUT_BED}"

    echo "✅ Arquivo criado: ${OUTPUT_BED}"
done

echo "🎯 Concluído para o modelo ${MODEL_NUM} estados."