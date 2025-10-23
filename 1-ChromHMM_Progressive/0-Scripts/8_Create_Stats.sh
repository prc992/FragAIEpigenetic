#!/bin/bash
#SBATCH --job-name=8_Create_Stats
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=2-00:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=BEGIN,END,FAIL


# ==========================
# CONFIG (hardcoded)
# ==========================
# CSV com as colunas: sample,mlfile
MANIFEST_CSV="/ABS/TO/your/mlready_list.csv"
# Saída (CSV long/tidy)
OUTPUT_CSV="$(dirname "$MANIFEST_CSV")/mlready_stats_long.csv"
OUTPUT_DIR="/data/baca/projects/ultima_deep/4-ChromHMM_Progressive/5-Features/1-ML_READY/"
# ==========================

if ! command -v gawk >/dev/null 2>&1; then
  echo "ERRO: gawk não encontrado. No macOS: brew install gawk" >&2
  exit 1
fi

echo "=== Gerando estatísticas em formato long ==="
echo "Manifest: $MANIFEST_CSV"
echo "Saída   : $OUTPUT_CSV"
echo "--------------------------------------------"

# Header do CSV de saída (long/tidy)
echo "section,sample,file,K,total_rows,label_set,metric,label_value,count,percent" > "$OUTPUT_CSV"

# Lê o manifest e processa cada arquivo MLREADY
{
  IFS=, read -r hdr1 hdr2 || true  # consome cabeçalho (sample,mlfile)
  lineno=1
  while IFS=, read -r sample mlfile; do
    lineno=$((lineno+1))
    # Pula linhas vazias ou comentários
    [ -z "${sample}${mlfile}" ] && continue
    case "$sample" in \#*|sample|"") continue ;; esac

    if [ ! -f "$mlfile" ]; then
      echo "ERRO (linha $lineno): arquivo não existe: $mlfile" >&2
      exit 1
    fi

    echo ">> Sample: $sample"
    echo "   File  : $mlfile"

    # gawk: calcula métricas e imprime linhas long
    gawk -v SAMPLE="$sample" -v FILEPATH="$mlfile" '
      BEGIN{
        FS=OFS="\t"
        # separador CSV para stdout
        sep=","
      }
      FNR==1 {
        # Mapeia cabeçalho
        for (i=1;i<=NF;i++) {
          H[$i]=i
          name=$i
          if (name ~ /^mask_/) nmask++
          if (name ~ /^label_[0-9]+_states$/) { labels[name]=i; nlabels++ }
        }
        # Infere K a partir da quantidade de máscaras: nmask = 2*K+1
        if (nmask>0) { K=(nmask-1)/2 } else { K="NA" }

        has_ctx = ("context_has_signal_count" in H) ? 1 : 0
        next
      }
      FNR>1 {
        total++

        # Contagem por conjunto de labels
        for (ls in labels) {
          v = $(labels[ls])
          if (v != "." && v != "") {
            rows_with_label[ls]++
            if (++cls_count[ls SUBSEP v] == 1) uniq_cls[ls]++
          }
        }

        # Histograma de contexto (opcional)
        if (has_ctx) {
          cv = $(H["context_has_signal_count"])
          ctx_count[cv]++
          ctx_total++
        }
      }
      END{
        # 1) Summary geral (total_rows)
        printf "summary%s%s%s%s%s%s%s%s%s%s%s%s%s\n", sep, SAMPLE, sep, FILEPATH, sep, K, sep, total, sep, "all", sep, "total_rows", sep, "", sep, total, sep, ""

        # 2) Summary por label_set: rows_with_label e unique_labels
        for (ls in labels) {
          rw = rows_with_label[ls] + 0
          uq = uniq_cls[ls] + 0
          pct = (total? 100.0*rw/total : 0)
          # rows_with_label
          printf "summary%s%s%s%s%s%s%s%s%s%s%s%s%s\n", sep, SAMPLE, sep, FILEPATH, sep, K, sep, total, sep, ls, sep, "rows_with_label", sep, "", sep, rw, sep, pct
          # unique_labels
          printf "summary%s%s%s%s%s%s%s%s%s%s%s%s%s\n", sep, SAMPLE, sep, FILEPATH, sep, K, sep, total, sep, ls, sep, "unique_labels", sep, "", sep, uq, sep, ""
        }

        # 3) Distribuição por classe (por label_set)
        for (key in cls_count) {
          split(key, a, SUBSEP); ls=a[1]; v=a[2]
          rw = rows_with_label[ls] + 0
          cnt = cls_count[key] + 0
          pct = (rw ? 100.0*cnt/rw : 0)
          printf "labels%s%s%s%s%s%s%s%s%s%s%s%s%s\n", sep, SAMPLE, sep, FILEPATH, sep, K, sep, total, sep, ls, sep, "", sep, v, sep, cnt, sep, pct
        }

        # 4) Histograma de context_has_signal_count (se existir)
        if (has_ctx) {
          for (c in ctx_count) {
            cnt = ctx_count[c] + 0
            pct = (ctx_total ? 100.0*cnt/ctx_total : 0)
            # section=context, label_set = context_has_signal_count, metric = hist_bin
            printf "context%s%s%s%s%s%s%s%s%s%s%s%s%s\n", sep, SAMPLE, sep, FILEPATH, sep, K, sep, total, sep, "context_has_signal_count", sep, "hist_bin", sep, c, sep, cnt, sep, pct
          }
        }
      }
    ' "$mlfile" >> "$OUTPUT_CSV"

  done
} < "$MANIFEST_CSV"

echo "✅ Pronto: $OUTPUT_CSV"