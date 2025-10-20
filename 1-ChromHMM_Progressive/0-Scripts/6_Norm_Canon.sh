#!/bin/bash
#SBATCH --job-name=6_Norm_Canon
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=2-00:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=BEGIN,END,FAIL

set +e
set +u
set +o pipefail
# (Opcional) logar cada comando para debug
set -x

export LC_ALL=C

# ====== DIRETÓRIOS ======
DIR_BASE="/PHShome/pd004/FragAIEpigenetic/1-ChromHMM_Progressive/4-Features/CANONICAL"
DIR_NORM="${DIR_BASE}/NORMALIZED"

mkdir -p logs "$DIR_NORM"

# Padrão dos arquivos a processar
PATTERN="*_FEATURES_BASE.txt"

shopt -s nullglob
files=( "${DIR_BASE}"/${PATTERN} )
shopt -u nullglob

if (( ${#files[@]} == 0 )); then
  echo "Nenhum arquivo encontrado em: ${DIR_BASE}/${PATTERN}" >&2
  exit 1
fi

echo "Processando ${#files[@]} arquivo(s) em: $DIR_BASE"
echo "Saída normalizada: $DIR_NORM"

for IN in "${files[@]}"; do
  [[ -r "$IN" ]] || { echo "Arquivo não legível: $IN" >&2; continue; }

  echo "→ $IN"

  # -------- 1) Ordenações para percentis (colunas 4 e 8) --------
  cov_sorted=$(mktemp)
  nrd_sorted=$(mktemp)
  tail -n +2 "$IN" | awk -F'\t' '{print $4}' | sort -n > "$cov_sorted"   # coverage_count
  tail -n +2 "$IN" | awk -F'\t' '{print $8}' | sort -n > "$nrd_sorted"   # n_reads

  N=$(wc -l < "$cov_sorted" | tr -d '[:space:]')

  # função percentil: índice = ceil(p*N)
  pctl() { # $1=sorted_file  $2=p (ex: 0.95)  $3=N
    local f="$1" p="$2" N="$3" idx
    idx=$(awk -v n="$N" -v p="$p" 'BEGIN{v=p*n; printf("%d",(v==int(v)?v:int(v)+1))}')
    (( idx < 1 )) && idx=1
    (( idx > N )) && idx=$N
    awk -v i="$idx" 'NR==i{print; exit}' "$f"
  }

  cov_p90=0; cov_p95=0; cov_p99=0
  nrd_p90=0; nrd_p95=0; nrd_p99=0
  if (( N > 0 )); then
    cov_p90=$(pctl "$cov_sorted" 0.90 "$N")
    cov_p95=$(pctl "$cov_sorted" 0.95 "$N")
    cov_p99=$(pctl "$cov_sorted" 0.99 "$N")
    nrd_p90=$(pctl "$nrd_sorted" 0.90 "$N")
    nrd_p95=$(pctl "$nrd_sorted" 0.95 "$N")
    nrd_p99=$(pctl "$nrd_sorted" 0.99 "$N")
  fi

  # -------- 2) Média e sd de mean_fragment_len (coluna 6) --------
  read mean sd nlines < <(awk -F'\t' 'NR>1{n++; s+=$6; ss+=$6*$6}
    END{
      if(n>0){ m=s/n; v=(n>1)?(ss/n - m*m):0; sd=(v>0?sqrt(v):0); printf "%.10f %.10f %d", m, sd, n }
      else    { printf "0 0 0" }
    }' "$IN")

  # -------- 3) Saídas --------
  base_name=$(basename "${IN%.txt}")
  OUT="${DIR_NORM}/${base_name}_normalized.txt"
  STATS="${DIR_NORM}/${base_name}_normalized.stats"

  # escreve header com novas colunas
  head -n1 "$IN" | awk 'BEGIN{OFS="\t"}{
    print $0,"coverage_norm_p95","coverage_norm_p99","coverage_log","nreads_norm","mean_fraglen_z"
  }' > "$OUT"

  # acrescenta linhas com colunas normalizadas no final
  awk -F'\t' -v OFS='\t' \
      -v cov95="$cov_p95" -v cov99="$cov_p99" \
      -v nrd95="$nrd_p95" -v nrd99="$nrd_p99" \
      -v mu="$mean" -v sig="$sd" \
      'NR==1{next}
       {
         # coverage norms e log
         covn95 = (cov95>0) ? ($4/cov95) : 0
         covn99 = (cov99>0) ? ($4/cov99) : 0
         covlog = log($4 + 1)

         # n_reads normalizado por p95
         nrm95  = (nrd95>0) ? ($8/nrd95) : 0

         # z-score de mean_fragment_len
         z      = (sig>0)  ? (($6-mu)/sig) : 0

         print $0, covn95, covn99, covlog, nrm95, z
       }' "$IN" >> "$OUT"

  # salva estatísticas usadas
  {
    echo "# arquivo: $IN"
    echo "# linhas (sem header): $N"
    echo "# coverage_count: p90=$cov_p90 p95=$cov_p95 p99=$cov_p99"
    echo "# n_reads       : p90=$nrd_p90 p95=$nrd_p95 p99=$nrd_p99"
    echo "# mean_fragment_len: mean=$mean sd=$sd (n=$nlines)"
  } > "$STATS"

  rm -f "$cov_sorted" "$nrd_sorted"
  echo "   ✅ Normalizado: $OUT"
  echo "   ℹ️  Stats:      $STATS"
done