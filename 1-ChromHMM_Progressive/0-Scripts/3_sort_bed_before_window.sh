#!/bin/bash
#SBATCH --job-name=3_sort_dense_bed
#SBATCH --array=2-8
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --time=10:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=ALL

set -euo pipefail

# ======== PARÂMETROS GERAIS ========
MODEL_NUM="${SLURM_ARRAY_TASK_ID}"   # número de estados (2–8)
BASE_DIR="/PHShome/pd004/FragAIEpigenetic/1-ChromHMM_Progressive"
INPUT_DIR="${BASE_DIR}/3-ChromHMM-Suport-Files/1-MODELS/0-MODEL_${MODEL_NUM}_STATES"
OUTPUT_DIR="${INPUT_DIR}/0-SORTED_HG19"
CHROM_ORDER="${BASE_DIR}/3-ChromHMM-Suport-Files/chrom_order_hg19.txt"

mkdir -p "${OUTPUT_DIR}" logs

echo "============================================================="
echo "🔹 Ordenando arquivos *dense.bed para o modelo ${MODEL_NUM} estados"
echo "Input:  ${INPUT_DIR}"
echo "Output: ${OUTPUT_DIR}"
echo "Chrom order: ${CHROM_ORDER}"
echo "============================================================="

# ======== LOCALIZAÇÃO DOS ARQUIVOS ========
mapfile -t BED_FILES < <(find "${INPUT_DIR}" -maxdepth 1 -type f -name "*dense.bed" | sort)
TOTAL_FILES=${#BED_FILES[@]}

if (( TOTAL_FILES == 0 )); then
    echo "⚠️ Nenhum arquivo *dense.bed encontrado em ${INPUT_DIR}"
    exit 0
fi

# ======== LOOP SOBRE ARQUIVOS ========
for (( i=0; i<${TOTAL_FILES}; i++ )); do
    INPUT_BED="${BED_FILES[$i]}"
    BASENAME=$(basename "${INPUT_BED}")
    OUTPUT_BED="${OUTPUT_DIR}/${BASENAME%.bed}_sorted.bed"

    echo "→ Ordenando: ${BASENAME}"

    # Ordenação respeitando a ordem cromossômica e posição
    TMP_FILE=$(mktemp)
    awk 'NR==FNR{a[$1]=++i; next} a[$1]{print a[$1], $0}' "${CHROM_ORDER}" "${INPUT_BED}" \
        | sort -k1,1n -k3,3n --parallel="${SLURM_CPUS_PER_TASK}" \
        | cut -d' ' -f2- > "${TMP_FILE}"

    mv "${TMP_FILE}" "${OUTPUT_BED}"
    echo "✅ Arquivo ordenado: ${OUTPUT_BED}"
done

echo "============================================================="
echo "🎯 Finalizado: todos os arquivos *dense.bed do modelo ${MODEL_NUM} ordenados"
echo "============================================================="