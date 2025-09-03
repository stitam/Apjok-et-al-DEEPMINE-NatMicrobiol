# GitHub does not store files larger than 100MB. Run this script in its own
# directory to merge chunks into a single table that can be used in the analysis.

df <- data.frame()
files <- paste0("blastx_chunk", 1:4, ".tsv")

for (i in files) {
  chunk <- read.csv(i, sep = "\t")
  df <- df |> dplyr::bind_rows(chunk)
}

write.table(
  df,
  file = "Nanopore sequences Run2 blastx.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
