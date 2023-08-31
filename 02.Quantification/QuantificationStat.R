library(tidyverse)

sampleInfo <- read_tsv(file = "../00.data/samples.txt", col_names = c("Group", "Sample", "fq_1", "fq_2"))

df <- data.frame(row.names = c("Assigned","Unassigned_Unmapped","Unassigned_Read_Type","Unassigned_Singleton","Unassigned_MappingQuality","Unassigned_Chimera","Unassigned_FragmentLength","Unassigned_Duplicate","Unassigned_MultiMapping","Unassigned_Secondary","Unassigned_NonSplit","Unassigned_NoFeatures","Unassigned_Overlapping_Length","Unassigned_Ambiguity"))
for (i in sampleInfo$Sample) {
  df_tmp <- read.table(file = paste(i, "log", sep = "."), row.names = 1, header = FALSE)
  colnames(df_tmp) <- i
  df <- cbind(df, df_tmp)
}
df <- t(df) %>% as.data.frame() %>% rownames_to_column(var = "Sample")
write_tsv(x = df, file = "Quantification_Stat.txt")
write_csv(x = df, file = "Quantification_Stat.csv")
