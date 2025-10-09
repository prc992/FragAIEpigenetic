#!/bin/bash
#SBATCH --job-name=5_extract-Features-windows
#SBATCH --array=2-3
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --time=5-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=BEGIN,END,FAIL

# Script unificado para extrair Coverage, GC Content e Size Features
# para cada conjunto de estados (K=2..8), em paralelo (job array).

set -euo pipefail
export LC_ALL=C

module load bedtools
module load samtools

# ======================= Configuração =======================
MODEL_NUM="${SLURM_ARRAY_TASK_ID}"   # número de estados (2–8)
BIN_SIZE=400
THREADS="${SLURM_CPUS_PER_TASK:-16}"


# Diretórios do modelo K
#DIR_CRAMS="/data/baca/projects/ultima_deep/2-data-samples_wgs_h19"
DIR_CRAMS="/data/baca/projects/ultima_deep/2-data-samples_wgs_h19"
DIR_BED_WINDOWS="/PHShome/pd004/FragAIEpigenetic/1-ChromHMM_Progressive/3-ChromHMM-Suport-Files/1-MODELS/0-MODEL_${MODEL_NUM}_STATES/1-WINDOWS-400BP"
#DIR_OUTPUT="/PHShome/pd004/FragAIEpigenetic/1-ChromHMM_Progressive/4-Features/MODEL_${MODEL_NUM}_STATES"
DIR_OUTPUT="/PHShome/pd004/FragAIEpigenetic/1-ChromHMM_Progressive/4-Features_small/MODEL_${MODEL_NUM}_STATES"
#CHECKPOINT_DIR="/PHShome/pd004/FragAIEpigenetic/1-ChromHMM_Progressive/4-Features/MODEL_${MODEL_NUM}_STATES/.checkpoints"
CHECKPOINT_DIR="/PHShome/pd004/FragAIEpigenetic/1-ChromHMM_Progressive/4-Features_small/MODEL_${MODEL_NUM}_STATES/.checkpoints"

mkdir -p logs "${DIR_OUTPUT}" "${CHECKPOINT_DIR}"

# Arquivos de referência
GENOME_SIZE="/data/baca/projects/ultima_deep/6-refs/hg19.genome"
GENOME_FASTA="/data/baca/projects/ultima_deep/6-refs/hg19.fa"
REF="/data/baca/projects/ultima_deep/6-refs/hg19.fa"

# Parâmetros para size features
MONO_MIN=120
MONO_MAX=180
DI_MIN=300
DI_MAX=380
PROG_DIVISOR=200
FLUSH_EVERY=1000
RESUME="${RESUME:-1}"

# ======================= PARES CRAM x BED =======================
# Liste aqui os CRAMs (em mesma ordem dos BEDs). Os BEDs são buscados por basename dentro de WINDOWS_DIR.
CRAM_FILES=(
  "${DIR_CRAMS}/M2673/M2673_merge_TrimAlignSort_hg19_small.cram"
  "${DIR_CRAMS}/M2707/M2707_merge_TrimAlignSort_hg19_small.cram"
  # Adicione mais CRAMs, um por linha
)

# Apenas os BASENAMES dos BEDs (devem existir em cada WINDOWS_DIR de cada K)
BED_BASENAMES=(
  "${DIR_BED_WINDOWS}/patient1_${MODEL_NUM}_dense_sorted_400bp.bed"
  "${DIR_BED_WINDOWS}/patient2_${MODEL_NUM}_dense_sorted_400bp.bed"
  # Adicione mais nomes, um por linha, na mesma ordem dos CRAMs
)

# ======================= Funções Auxiliares ======================
fmt_time() {
  s="${1:-0}"
  if [ "$s" -lt 0 ] 2>/dev/null; then echo "NA"; else
    h=$((s/3600)); m=$(((s%3600)/60)); sec=$((s%60))
    printf "%02dh%02dm%02ds" "$h" "$m" "$sec"
  fi
}
log_info(){ echo "[$(date '+%F %T')] [INFO] $*" >&2; }
log_error(){ echo "[$(date '+%F %T')] [ERROR] $*" >&2; exit 1; }
check_file(){ [[ -r "$1" ]] || log_error "$2 não encontrado/sem leitura: $1"; }

# ======================= Checagens =======================
[[ ${#CRAM_FILES[@]} -eq ${#BED_BASENAMES[@]} ]] || \
  log_error "CRAM_FILES e BED_BASENAMES têm tamanhos diferentes."

check_file "${GENOME_SIZE}" "genome size"
check_file "${GENOME_FASTA}" "genoma FASTA"
check_file "${REF}" "referência"

# Índices se faltarem
[[ -s "${REF}.fai" ]] || { log_info "Criando índice .fai..."; samtools faidx "${REF}"; }
# Índices CRAI por CRAM serão conferidos no loop

echo "============================================================="
echo "🔹 Extraindo features | K=${MODEL_NUM} estados"
echo "Windows dir: ${DIR_BED_WINDOWS}"
echo "Saída:       ${DIR_OUTPUT}"
echo "Threads:     ${THREADS}"
echo "============================================================="

# ======================= Loop Principal =======================
N=${#CRAM_FILES[@]}
for idx in $(seq 0 $((N-1))); do
  CRAM="${CRAM_FILES[$idx]}"
  BED="${BED_BASENAMES[$idx]}"

  SAMPLE_NAME=$(basename "$CRAM" .cram)

  # Saídas temporárias e final
  OUTPUT_COVERAGE_TMP="${CHECKPOINT_DIR}/${SAMPLE_NAME}_coverage_tmp.txt"
  OUTPUT_GC_TMP="${CHECKPOINT_DIR}/${SAMPLE_NAME}_gc_tmp.txt"
  OUTPUT_SIZE_TMP="${CHECKPOINT_DIR}/${SAMPLE_NAME}_size_tmp.txt"
  OUTPUT_FINAL="${DIR_OUTPUT}/${SAMPLE_NAME}_${MODEL_NUM}_STATES_${BIN_SIZE}bp_ALL_FEATURES.txt"

  log_info ">>> Par $((idx+1))/$N | K=${MODEL_NUM}"
  log_info "    CRAM: ${CRAM}"
  log_info "    BED : ${BED}"
  log_info "    OUT : ${OUTPUT_FINAL}"

  check_file "${CRAM}" "CRAM"
  check_file "${BED}" "BED"

  # Índice CRAI
  if [[ ! -s "${CRAM}.crai" ]]; then
    log_info "Criando índice CRAI..."
    samtools index -@ "${THREADS}" "${CRAM}"
  fi

  START_TS=$(date +%s)

  # ======================= 1. COVERAGE =======================
  log_info "Calculando coverage..."
  TS1=$(date +%s)
  bedtools coverage -a "${BED}" -b "${CRAM}" -sorted -g "${GENOME_SIZE}" -counts > "${OUTPUT_COVERAGE_TMP}"
  TS2=$(date +%s)
  log_info "Coverage: $(fmt_time $((TS2-TS1)))"

  # ======================= 2. GC CONTENT =======================
  log_info "Calculando GC content..."
  TS1=$(date +%s)
  bedtools nuc -fi "${GENOME_FASTA}" -bed "${BED}" > "${OUTPUT_GC_TMP}"
  TS2=$(date +%s)
  log_info "GC: $(fmt_time $((TS2-TS1)))"

  # ======================= 3. SIZE FEATURES =======================
  log_info "Calculando size features..."
  TS1=$(date +%s)

  echo -e "chrom\tstart\tend\tmean_fragment_len\tsd_fragment_len\tn_reads\tmono_frac\tdi_frac\tmono_to_di_ratio\tentropy_norm_1bp\tentropy_norm_5bp\tentropy_norm_15bp\tentropy_norm_20bp" > "${OUTPUT_SIZE_TMP}"

  # Cromossomos presentes nas janelas
  mapfile -t CHRS < <(awk '!/^#/ {print $1}' "${BED}" | awk 'prev!=$0{print; prev=$0}')
  TOTAL_WINDOWS=$(grep -vc '^#' "${BED}")
  PROG_STEP=$(( TOTAL_WINDOWS / PROG_DIVISOR )); (( PROG_STEP < 1 )) && PROG_STEP=1

  WINDOWS_DONE=0
  for CHR in "${CHRS[@]}"; do
    log_info "  → chr ${CHR}"
    CHR_BED=$(mktemp)
    awk -v C="$CHR" '!/^#/ && $1==C {print $0}' "${BED}" > "${CHR_BED}"
    CHR_WINDOWS=$(wc -l < "${CHR_BED}")
    if (( CHR_WINDOWS == 0 )); then rm -f "${CHR_BED}"; continue; fi

    samtools view -@ "${THREADS}" -T "${REF}" -L "${CHR_BED}" "${CRAM}" "${CHR}" \
    | awk -v BEDFILE="${CHR_BED}" -v OUT="${OUTPUT_SIZE_TMP}" \
          -v MONO_MIN="${MONO_MIN}" -v MONO_MAX="${MONO_MAX}" \
          -v DI_MIN="${DI_MIN}"     -v DI_MAX="${DI_MAX}" \
          -v TOTAL="${TOTAL_WINDOWS}" -v STEP="${PROG_STEP}" \
          -v OFFSET="${WINDOWS_DONE}" -v CHR="${CHR}" \
          -v FLUSH_EVERY="${FLUSH_EVERY}" '
        function log2(x){ return log(x)/log(2) }
        function bin_floor(len, w){ return int(len / w) * w }
        function do_flush(){ if(FLUSH_EVERY>0) fflush(OUT) }
        function flush_window(j,   mean,sd,var,mfrac,dfrac,ratio,H1,H1n,H5,H5n,H15,H15n,H20,H20n,p,k) {
          if(n>0){
            mean = sum / n
            if(n>1){ var=(sumsq/n)-(mean*mean); sd=(var>0?sqrt(var):0) } else sd="NA"
            mfrac = mono/n; dfrac=di/n; ratio=(mono>0 && di>0)?mono/di:"NA"
            H1=0;k=0; for(l in cnt1){p=cnt1[l]/n; H1-=p*log2(p); k++}   H1n=(k>1)?H1/log2(k):"NA"
            H5=0;k=0; for(l in cnt5){p=cnt5[l]/n; H5-=p*log2(p); k++}   H5n=(k>1)?H5/log2(k):"NA"
            H15=0;k=0;for(l in cnt15){p=cnt15[l]/n;H15-=p*log2(p);k++} H15n=(k>1)?H15/log2(k):"NA"
            H20=0;k=0;for(l in cnt20){p=cnt20[l]/n;H20-=p*log2(p);k++} H20n=(k>1)?H20/log2(k):"NA"
            printf("%s\t%d\t%d\t%.6f\t%s\t%d\t%.6f\t%.6f\t%s\t%s\t%s\t%s\t%s\n",
                   CHR, start[j], end[j], mean, sd, n, mfrac, dfrac, ratio, H1n, H5n, H15n, H20n) >> OUT
          } else {
            printf("%s\t%d\t%d\tNA\tNA\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", CHR, start[j], end[j]) >> OUT
          }
          delete cnt1; delete cnt5; delete cnt15; delete cnt20
          sum=0; sumsq=0; n=0; mono=0; di=0; printed++
          if(FLUSH_EVERY>0 && (printed % FLUSH_EVERY)==0) do_flush()
          done = OFFSET + printed
          if( (done % STEP)==0 || done==TOTAL ){
            pct = (100.0*done)/TOTAL
            printf("[PROGRESS] %d/%d (%.1f%%) | chr=%s\n", done, TOTAL, pct, CHR) > "/dev/stderr"
          }
        }
        BEGIN{
          j=1
          while( (getline line < BEDFILE) > 0 ){
            if(line ~ /^#/) continue
            split(line,a,"\t"); start[j]=a[2]+0; end[j]=a[3]+0; j++
          }
          maxj=j-1; close(BEDFILE)
          sum=0; sumsq=0; n=0; mono=0; di=0; printed=0; j=1
        }
        {
          pos = $4 + 0
          while(j<=maxj && pos >= end[j]) { flush_window(j); j++ }
          if(j>maxj) next
          if(pos >= start[j] && pos < end[j]){
            tm=0; a3=""
            for(i=12;i<=NF;i++){
              if($i ~ /^tm:Z:A/) tm=1
              else if($i ~ /^a3:i:/){ split($i,t,":"); a3=t[3] }
            }
            if(tm){
              len = (a3!="") ? a3 : length($10)
              sum += len; sumsq += len*len; n++
              if(len>=MONO_MIN && len<=MONO_MAX) mono++
              if(len>=DI_MIN   && len<=DI_MAX)   di++
              cnt1[len]++; cnt5[bin_floor(len,5)]++
              cnt15[bin_floor(len,15)]++; cnt20[bin_floor(len,20)]++
            }
          }
        }
        END{ while(j<=maxj){ flush_window(j); j++ } do_flush() }
    '
    WINDOWS_DONE=$(( WINDOWS_DONE + CHR_WINDOWS ))
    rm -f "${CHR_BED}"
  done

  TS2=$(date +%s)
  log_info "Size: $(fmt_time $((TS2-TS1)))"

  # ======================= 4. COMBINAR TODAS AS FEATURES =======================
  log_info "Combinando features..."
  TS1=$(date +%s)

  # Cabeçalho combinado
  echo -e "chrom\tstart\tend\tstate\tcoverage_count\tpct_gc\tmean_fragment_len\tsd_fragment_len\tn_reads\tmono_frac\tdi_frac\tmono_to_di_ratio\tentropy_norm_1bp\tentropy_norm_5bp\tentropy_norm_15bp\tentropy_norm_20bp" > "${OUTPUT_FINAL}"

  # pct_gc = coluna 6 de bedtools nuc (skip header)
  # size features = colunas 4- fim do OUTPUT_SIZE_TMP (já com header)
  paste "${OUTPUT_COVERAGE_TMP}" \
        <(tail -n +2 "${OUTPUT_GC_TMP}" | cut -f6) \
        <(tail -n +2 "${OUTPUT_SIZE_TMP}" | cut -f4-) \
  | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' >> "${OUTPUT_FINAL}"

  TS2=$(date +%s)
  log_info "Combinação: $(fmt_time $((TS2-TS1)))"

  ELAP=$(( $(date +%s) - START_TS ))
  log_info "✅ Completo em $(fmt_time "${ELAP}")"
  log_info "📄 Arquivo final: ${OUTPUT_FINAL}"
  log_info "👀 Preview:"
  head -27 "${OUTPUT_FINAL}" >&2

  # Limpeza (se quiser manter para debug, comente)
  # rm -f "${OUTPUT_COVERAGE_TMP}" "${OUTPUT_GC_TMP}" "${OUTPUT_SIZE_TMP}"

done

log_info "🎉 Finalizado para K=${MODEL_NUM}."