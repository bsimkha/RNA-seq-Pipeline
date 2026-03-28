#!/usr/bin/env Rscript

library(yaml)

# ================================
# Part 1: Read config + find files
# ================================

cfg <- yaml::read_yaml("/Volumes/Bibhu's SSD/RNA_seq/Pipeline/config.yaml")

input_dir <- cfg$input_dir
pattern   <- cfg$files$pattern

files <- Sys.glob(file.path(input_dir, pattern))

if (length(files) == 0) {
  stop("No FASTQ files found")
}

# ================================
# Part 2: Classify files (R1 / R2)
# ================================

base <- basename(files)

mate <- ifelse(
  grepl("(_R?1|_1)\\.(f(ast)?q)(\\.gz)?$", base, ignore.case = TRUE), "1",
  ifelse(
    grepl("(_R?2|_2)\\.(f(ast)?q)(\\.gz)?$", base, ignore.case = TRUE), "2",
    NA
  )
)

sample_id <- sub(
  "(_R?1|_R?2|_1|_2)\\.(f(ast)?q)(\\.gz)?$",
  "",
  base,
  ignore.case = TRUE
)

df <- data.frame(
  file = files,
  sample_id = sample_id,
  mate = mate,
  stringsAsFactors = FALSE
)

df <- df[!is.na(df$mate), ]

# ================================
# Part 3: Group into paired samples
# ================================

groups <- split(df, df$sample_id)

paired <- list()
unpaired <- character()

for (sid in names(groups)) {
  g <- groups[[sid]]
  
  if (all(c("1", "2") %in% g$mate)) {
    paired[[sid]] <- g
  } else {
    unpaired <- c(unpaired, g$file)
  }
}

# ================================
# Print summary + preview
# ================================

cat("\nPaired samples detected:", length(paired), "\n")
cat("Unmatched files:", length(unpaired), "\n\n")

for (sid in names(paired)) {
  g <- paired[[sid]]
  
  r1 <- basename(g$file[g$mate == "1"])
  r2 <- basename(g$file[g$mate == "2"])
  
  cat(sid, "\n")
  cat("  R1:", paste(r1, collapse = ", "), "\n")
  cat("  R2:", paste(r2, collapse = ", "), "\n\n")
}