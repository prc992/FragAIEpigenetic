#!/bin/bash
#SBATCH --job-name=7_Create_Features
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=2-00:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=BEGIN,END,FAIL

# ==========================
# CONFIGURAÇÃO FIXA
# ==========================
# Caminho do CSV (hardcoded)
MANIFEST_CSV="/PHShome/pd004/FragAIEpigenetic/1-ChromHMM_Progressive/0-Scripts/7_Path_Samples.csv"

# Diretório de saída
OUTPUT_DIR="/data/baca/projects/ultima_deep/4-ChromHMM_Progressive/5-Features/1-ML_READY"

mkdir -p "$OUTPUT_DIR" logs

# Parâmetros do contexto e filtros
K=5
DESERT_MIN_RUN=10
MIN_NEIGHBORS=5

# ==========================
echo "=== Prep ML (via manifest CSV) ==="
echo "Manifest : $MANIFEST_CSV"
echo "Output   : $OUTPUT_DIR"
echo "K=±$K | desert≥$DESERT_MIN_RUN | vizinhos≥$MIN_NEIGHBORS"
echo "----------------------------------"

# Lê CSV linha a linha
# Espera cabeçalho: sample,features,label2,label3,label4,label5,label6,label7,label8
{
  read -r header || true
  line_no=1
  while IFS=, read -r sample features l2 l3 l4 l5 l6 l7 l8; do
    line_no=$((line_no+1))
    [ -z "${sample}${features}${l2}${l3}${l4}${l5}${l6}${l7}${l8}" ] && continue
    case "$sample" in \#*|sample) continue ;; esac

    # Validação mínima
    for f in "$features"; do
      if [ ! -f "$f" ]; then
        echo "ERRO (linha $line_no): arquivo de features não existe: $f" >&2
        exit 1
      fi
    done

    # Nome de saída
    if [ -n "${sample:-}" ]; then
      base="$sample"
    else
      fname=$(basename -- "$features")
      base="${fname%.*}"
    fi
    out="${OUTPUT_DIR}/${base}.tsv"

    echo ">> Sample: ${sample:-<basename>}"
    echo "   Features: $features"
    echo "   Labels  : $l2 | $l3 | $l4 | $l5 | $l6 | $l7 | $l8"
    echo "   Saída   : $out"

    gawk -v K="$K" -v DESERT_MIN_RUN="$DESERT_MIN_RUN" -v MIN_NEIGHBORS="$MIN_NEIGHBORS" '
    BEGIN{
      FS=OFS="\t"
      header_done=0
      split("coverage_count,pct_gc,mean_fragment_len,sd_fragment_len,n_reads,mono_frac,di_frac,mono_to_di_ratio,entropy_norm_1bp,entropy_norm_5bp,entropy_norm_15bp,entropy_norm_20bp,coverage_norm_p95,coverage_norm_p99,coverage_log,nreads_norm,mean_fraglen_z", FEATNAMES, ",")
      nfeat=length(FEATNAMES)
    }

    # ----- labels 2..8 -----
    ARGIND==1 && FILENAME!="" { key=$1 OFS $2 OFS $3; lab2[key]=$4; next }
    ARGIND==2 && FILENAME!="" { key=$1 OFS $2 OFS $3; lab3[key]=$4; next }
    ARGIND==3 && FILENAME!="" { key=$1 OFS $2 OFS $3; lab4[key]=$4; next }
    ARGIND==4 && FILENAME!="" { key=$1 OFS $2 OFS $3; lab5[key]=$4; next }
    ARGIND==5 && FILENAME!="" { key=$1 OFS $2 OFS $3; lab6[key]=$4; next }
    ARGIND==6 && FILENAME!="" { key=$1 OFS $2 OFS $3; lab7[key]=$4; next }
    ARGIND==7 && FILENAME!="" { key=$1 OFS $2 OFS $3; lab8[key]=$4; next }

    # ----- FEATURES -----
    ARGIND==8 && FNR==1 {
      for(i=1;i<=NF;i++) H[$i]=i
      req="chrom start end coverage_count n_reads"
      n=split(req, R, " ")
      for(i=1;i<=n;i++) if(!(R[i] in H)) { print "ERRO: falta coluna " R[i] > "/dev/stderr"; exit 2 }
      for (f=1; f<=nfeat; f++) if (!(FEATNAMES[f] in H)) {
        print "ERRO: feature ausente: " FEATNAMES[f] > "/dev/stderr"; exit 3
      }
      next
    }

    ARGIND==8 {
      c=$(H["chrom"]); s=$(H["start"]); e=$(H["end"])
      if (chrom_buf!="" && c!=chrom_buf) flush_chrom_()
      if (chrom_buf=="") chrom_buf=c

      idx_top++
      start[idx_top]=s; end[idx_top]=e
      cov=$(H["coverage_count"]); nrd=$(H["n_reads"])
      if (cov==""||cov=="NA") cov=0
      if (nrd==""||nrd=="NA") nrd=0
      has_signal[idx_top]=((cov+0)>0||(nrd+0)>0)?1:0

      for (f=1; f<=nfeat; f++) {
        col=H[FEATNAMES[f]]
        V[idx_top,col]=$(col)
      }
      next
    }

    function flush_chrom_(   i,j,run_start,run_len,run,neigh_on,posname,fname,col,val,ii,key,
                             y2,y3,y4,y5,y6,y7,y8,h2,h3,h4,h5,h6,h7,h8,hasAny) {
      if (idx_top==0) return
      for (i=1;i<=idx_top;i++) desert[i]=0
      i=1
      while (i<=idx_top) {
        if (has_signal[i]==0) {
          run_start=i; run_len=0; j=i
          while (j<=idx_top && has_signal[j]==0) { run_len++; j++ }
          if (run_len>=DESERT_MIN_RUN) for (run=run_start; run<j; run++) desert[run]=1
          i=j
        } else i++
      }

      if (!header_done) {
        printf "chrom\tstart\tend\tlabel_2_states\tlabel_3_states\tlabel_4_states\tlabel_5_states\tlabel_6_states\tlabel_7_states\tlabel_8_states"
        for (o=-K;o<=K;o++) {
          posname=(o<0)?sprintf("m%d",-o):(o==0?"p0":sprintf("p%d",o))
          printf "\tmask_%s", posname
        }
        for (f=1; f<=nfeat; f++) {
          fname=FEATNAMES[f]
          for (o=-K;o<=K;o++) {
            posname=(o<0)?sprintf("m%d",-o):(o==0?"p0":sprintf("p%d",o))
            printf "\t%s_%s", fname, posname
          }
        }
        printf "\tcontext_has_signal_count\n"
        header_done=1
      }

      for (i=1;i<=idx_top;i++) {
        if (desert[i]==1) continue
        key = chrom_buf OFS start[i] OFS end[i]
        h2=(key in lab2); y2=h2?lab2[key]:"."
        h3=(key in lab3); y3=h3?lab3[key]:"."
        h4=(key in lab4); y4=h4?lab4[key]:"."
        h5=(key in lab5); y5=h5?lab5[key]:"."
        h6=(key in lab6); y6=h6?lab6[key]:"."
        h7=(key in lab7); y7=h7?lab7[key]:"."
        h8=(key in lab8); y8=h8?lab8[key]:"."
        hasAny=(h2||h3||h4||h5||h6||h7||h8)
        if (!hasAny) continue

        neigh_on=0
        for (o=1;o<=K;o++) {
          if ((i-o)>=1 && has_signal[i-o]==1) neigh_on++
          if ((i+o)<=idx_top && has_signal[i+o]==1) neigh_on++
        }
        if (neigh_on < MIN_NEIGHBORS) continue

        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", chrom_buf, start[i], end[i], y2, y3, y4, y5, y6, y7, y8
        for (o=-K;o<=K;o++) {
          ii=i+o
          if (ii>=1 && ii<=idx_top) printf "\t%d", has_signal[ii]
          else printf "\t0"
        }
        for (f=1; f<=nfeat; f++) {
          col=H[FEATNAMES[f]]
          for (o=-K;o<=K;o++) {
            ii=i+o
            if (ii>=1 && ii<=idx_top) {
              val=V[ii,col]
              if (val==""||val=="NA") printf "\tNA"; else printf "\t%s", val
            } else printf "\tNA"
          }
        }
        printf "\t%d\n", neigh_on
      }
      idx_top=0; chrom_buf=""
      delete start; delete end; delete has_signal; delete desert
    }

    END{ flush_chrom_() }
    ' "$l2" "$l3" "$l4" "$l5" "$l6" "$l7" "$l8" "$features" > "$out"

    echo "✅ Gerado: $out"
  done
} < "$MANIFEST_CSV"

echo "== Concluído =="