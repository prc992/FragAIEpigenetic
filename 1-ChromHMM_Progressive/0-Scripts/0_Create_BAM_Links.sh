# Diretórios de origem e destino
src_dir="/data/baca/projects/ultima_deep/2-data-samples-chip-seq"
dest_dir="/PHShome/pd004/FragAIEpigenetic/1-ChromHMM_Progressive/2-raw-files-sl"

# Patient 1 (25M2673)
ln -s "$src_dir/DV_25M2673K27_1/DV_25M2673K27_1_unique.sorted.dedup.bam" "$dest_dir/patient1_H3K27ac.bam"
ln -s "$src_dir/DV_25M2673K4_1/DV_25M2673K4_1_unique.sorted.dedup.bam" "$dest_dir/patient1_H3K4me3.bam"
ln -s "$src_dir/mHS_GR_M2673/mHS_GR_M2673_unique.sorted.dedup.bam" "$dest_dir/patient1_Medip.bam"

# Patient 2 (25M2707)
ln -s "$src_dir/DV_25M2707K27_1/DV_25M2707K27_1_unique.sorted.dedup.bam" "$dest_dir/patient2_H3K27ac.bam"
ln -s "$src_dir/DV_25M2707K4_1/DV_25M2707K4_1_unique.sorted.dedup.bam" "$dest_dir/patient2_H3K4me3.bam"
ln -s "$src_dir/mHS_GR_M2707/mHS_GR_M2707_unique.sorted.dedup.bam" "$dest_dir/patient2_Medip.bam"

# Patient 3 (30M2594)
ln -s "$src_dir/DV_30M2594K27_1/DV_30M2594K27_1_unique.sorted.dedup.bam" "$dest_dir/patient3_H3K27ac.bam"
ln -s "$src_dir/DV_30M2594K4_1/DV_30M2594K4_1_unique.sorted.dedup.bam" "$dest_dir/patient3_H3K4me3.bam"
ln -s "$src_dir/mHS_GR_M2594/mHS_GR_M2594_unique.sorted.dedup.bam" "$dest_dir/patient3_Medip.bam"

# Patient 4_60 (41M2560_60)
ln -s "$src_dir/DV_41M2560_60K27_1/DV_41M2560_60K27_1_unique.sorted.dedup.bam" "$dest_dir/patient4_60_H3K27ac.bam"
ln -s "$src_dir/DV_41M2560_60K4_1/DV_41M2560_60K4_1_unique.sorted.dedup.bam" "$dest_dir/patient4_60_H3K4me3.bam"
ln -s "$src_dir/mRN_GR_M2560_60/mRN_GR_M2560_60_unique.sorted.dedup.bam" "$dest_dir/patient4_60_Medip.bam"

# Patient 4_120 (41M2560_120)
ln -s "$src_dir/DV_41M2560_120K27_1/DV_41M2560_120K27_1_unique.sorted.dedup.bam" "$dest_dir/patient4_120_H3K27ac.bam"
ln -s "$src_dir/DV_41M2560_120K4_1/DV_41M2560_120K4_1_unique.sorted.dedup.bam" "$dest_dir/patient4_120_H3K4me3.bam"
ln -s "$src_dir/mRN_GR_M2560_120/mRN_GR_M2560_120_unique.sorted.dedup.bam" "$dest_dir/patient4_120_Medip.bam"

# Patient 5 (41M2671_120)
ln -s "$src_dir/DV_41M2671_120K27_1/DV_41M2671_120K27_1_unique.sorted.dedup.bam" "$dest_dir/patient5_H3K27ac.bam"
ln -s "$src_dir/DV_41M2671_120K4_1/DV_41M2671_120K4_1_unique.sorted.dedup.bam" "$dest_dir/patient5_H3K4me3.bam"
ln -s "$src_dir/mRN_GR_M2671_120/mRN_GR_M2671_120_unique.sorted.dedup.bam" "$dest_dir/patient5_Medip.bam"

