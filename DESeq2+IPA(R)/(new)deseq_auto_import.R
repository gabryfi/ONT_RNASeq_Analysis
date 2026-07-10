#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(stringr)
  library(tools)
  library(pheatmap)
  library(ggplot2)
  library(matrixStats)
  library(RColorBrewer)
  library(dplyr)
})

# ==============================================================================
# PARAMETRI DA TERMINALE
#
# Uso:
#   Rscript deseq_auto_import.R <WORK_DIR> <TREATMENT_LABEL> <CONTROL_LABEL>
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    "Uso corretto: Rscript deseq_auto_import.R ",
    "<WORK_DIR> <TREATMENT_LABEL> <CONTROL_LABEL>"
  )
}

work_dir  <- normalizePath(args[1], mustWork = TRUE)
treatment <- args[2]
control   <- args[3]

cat("[INFO] Cartella di lavoro: ", work_dir, "\n", sep = "")
cat("[INFO] Contrasto: ", treatment, " vs ", control, "\n", sep = "")

# ==============================================================================
# 1. LETTURA DEI FILE DI CONTEGGIO
# ==============================================================================

count_files <- sort(
  list.files(
    work_dir,
    pattern = "_counts\\.txt$",
    full.names = TRUE
  )
)

if (length(count_files) == 0) {
  stop(
    "Nessun file *_counts.txt trovato nella cartella: ",
    work_dir
  )
}

# Estrae i metadati da nomi nel formato:
#   <BATCH>_barcode_<BARCODE>_<CONDITION>_counts.txt
#
# Il nome del batch può contenere underscore.
extract_metadata <- function(filename) {
  base_name <- basename(filename)
  stem <- file_path_sans_ext(base_name)

  match <- str_match(
    stem,
    "^(.+)_barcode_([^_]+)_([^_]+)_counts$"
  )

  if (any(is.na(match))) {
    stop(
      "Il nome del file non segue il formato atteso ",
      "<BATCH>_barcode_<BARCODE>_<CONDITION>_counts.txt: ",
      base_name
    )
  }

  batch <- match[1, 2]
  barcode <- match[1, 3]
  condition <- match[1, 4]
  raw_sample <- stem
  sample <- paste(batch, "bc", barcode, condition, sep = "_")

  data.frame(
    raw_sample = raw_sample,
    sample = sample,
    batch = batch,
    barcode = barcode,
    condition = condition,
    file = filename,
    stringsAsFactors = FALSE
  )
}

sample_table <- do.call(
  rbind,
  lapply(count_files, extract_metadata)
)

if (anyDuplicated(sample_table$raw_sample)) {
  stop("Sono presenti nomi di campione duplicati nei file di conteggio.")
}

if (anyDuplicated(sample_table$sample)) {
  stop("Sono presenti identificatori di campione duplicati.")
}

# Legge un file prodotto da featureCounts.
read_featurecounts <- function(file_path) {
  df <- read.table(
    file_path,
    header = TRUE,
    skip = 1,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    comment.char = ""
  )

  if (nrow(df) == 0 || ncol(df) == 0) {
    stop("File featureCounts vuoto o non valido: ", file_path)
  }

  count_values <- df[[ncol(df)]]

  if (!is.numeric(count_values)) {
    stop(
      "La colonna dei conteggi non è numerica nel file: ",
      file_path
    )
  }

  if (anyNA(count_values)) {
    stop("Sono presenti conteggi NA nel file: ", file_path)
  }

  if (any(count_values < 0)) {
    stop("Sono presenti conteggi negativi nel file: ", file_path)
  }

  result <- data.frame(
    count = count_values,
    row.names = rownames(df),
    check.names = FALSE
  )

  colnames(result) <- file_path_sans_ext(basename(file_path))
  result
}

count_list <- lapply(sample_table$file, read_featurecounts)

# Verifica che tutti i file contengano gli stessi geni nello stesso ordine.
reference_gene_ids <- rownames(count_list[[1]])

for (i in seq_along(count_list)) {
  if (!identical(rownames(count_list[[i]]), reference_gene_ids)) {
    stop(
      "L'ordine o l'insieme dei Gene_ID non coincide nel file: ",
      sample_table$file[i]
    )
  }
}

count_matrix <- do.call(cbind, count_list)

storage.mode(count_matrix) <- "numeric"
count_matrix <- round(count_matrix)

if (anyNA(count_matrix)) {
  stop("La matrice dei conteggi contiene valori NA.")
}

if (any(count_matrix < 0)) {
  stop("La matrice dei conteggi contiene valori negativi.")
}

# Allinea sample_table e count_matrix.
sample_table <- sample_table[
  match(colnames(count_matrix), sample_table$raw_sample),
  ,
  drop = FALSE
]

if (anyNA(sample_table$raw_sample)) {
  stop(
    "Impossibile allineare i nomi delle colonne della matrice ",
    "con i metadati dei campioni."
  )
}

rownames(sample_table) <- sample_table$raw_sample

if (!identical(rownames(sample_table), colnames(count_matrix))) {
  stop(
    "I nomi delle righe dei metadati non coincidono con ",
    "le colonne della matrice dei conteggi."
  )
}

# Controlla che trattamento e controllo siano presenti.
available_conditions <- unique(sample_table$condition)

if (!control %in% available_conditions) {
  stop(
    "La condizione di controllo '",
    control,
    "' non è presente nei file di conteggio."
  )
}

if (!treatment %in% available_conditions) {
  stop(
    "La condizione di trattamento '",
    treatment,
    "' non è presente nei file di conteggio."
  )
}

sample_table$batch <- factor(sample_table$batch)

sample_table$condition <- relevel(
  factor(sample_table$condition),
  ref = control
)

# ==============================================================================
# 2. COSTRUZIONE DELL'OGGETTO DESEQ2 E FILTRAGGIO
# ==============================================================================

if (nlevels(sample_table$batch) > 1) {
  design_formula <- ~ batch + condition

  cat(
    "[INFO] Rilevati più batch. ",
    "Modello applicato: ~ batch + condition\n",
    sep = ""
  )
} else {
  design_formula <- ~ condition

  cat(
    "[INFO] Rilevato un singolo batch. ",
    "Modello applicato: ~ condition\n",
    sep = ""
  )
}

col_data <- sample_table[
  ,
  c("sample", "batch", "condition"),
  drop = FALSE
]

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = col_data,
  design = design_formula
)

cat(
  "[INFO] Applicazione filtro: ",
  "rimozione dei geni con rowSums < 10...\n",
  sep = ""
)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

if (nrow(dds) == 0) {
  stop(
    "Nessun gene supera il filtro rowSums >= 10. ",
    "Controllare i file di conteggio."
  )
}

cat(
  "[INFO] Geni mantenuti dopo il filtro: ",
  nrow(dds),
  "\n",
  sep = ""
)

# ==============================================================================
# 3. ANALISI STATISTICA
# ==============================================================================

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("condition", treatment, control)
)

res_df <- as.data.frame(res)

write.csv(
  res_df,
  file = file.path(work_dir, "deseq_results.csv")
)

write.csv(
  sample_table,
  file = file.path(work_dir, "samples_metadata.csv"),
  row.names = FALSE
)

# ==============================================================================
# 4. MATRICE PER IPA QIAGEN
# ==============================================================================

cat(
  "[INFO] Generazione del file di input ottimizzato per IPA Qiagen...\n"
)

res_df$Fold_Change_Linear <- ifelse(
  is.na(res_df$log2FoldChange),
  NA_real_,
  ifelse(
    res_df$log2FoldChange >= 0,
    2^(res_df$log2FoldChange),
    -1 / (2^(res_df$log2FoldChange))
  )
)

contrast_suffix <- paste0(treatment, "_vs_", control)

ipa_input <- data.frame(
  Gene_ID = rownames(res_df),
  Expr_Log2_Fold_Change = res_df$log2FoldChange,
  Expr_Fold_Change = res_df$Fold_Change_Linear,
  Expr_p_value = res_df$pvalue,
  Expr_FDR = res_df$padj,
  stringsAsFactors = FALSE
)

colnames(ipa_input) <- c(
  "Gene_ID",
  paste0("Expr Log2 Fold Change ", contrast_suffix),
  paste0("Expr Fold Change ", contrast_suffix),
  paste0("Expr p-value ", contrast_suffix),
  paste0("Expr FDR ", contrast_suffix)
)

write.table(
  ipa_input,
  file = file.path(
    work_dir,
    paste0("IPA_ready_input_", contrast_suffix, ".txt")
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ==============================================================================
# 5. PCA E HEATMAP
# ==============================================================================

rld <- rlog(dds, blind = FALSE)

# PCA.
pca_data <- plotPCA(
  rld,
  intgroup = "condition",
  returnData = TRUE
)

percent_var <- round(
  100 * attr(pca_data, "percentVar")
)

p_pca <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    color = condition
  )
) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percent_var[1], "%")) +
  ylab(paste0("PC2: ", percent_var[2], "%")) +
  theme_minimal() +
  ggtitle("PCA Plot")

ggsave(
  filename = file.path(work_dir, "PCA_plot.png"),
  plot = p_pca,
  width = 6,
  height = 5
)

# Heatmap dei geni più variabili.
variance_values <- matrixStats::rowVars(
  assay(rld)
)

number_top_genes <- min(
  30L,
  length(variance_values)
)

if (number_top_genes >= 2L) {
  top_genes <- head(
    order(
      variance_values,
      decreasing = TRUE
    ),
    number_top_genes
  )

  mat <- assay(rld)[top_genes, , drop = FALSE]
  mat <- mat - rowMeans(mat)

  annotation_columns <- as.data.frame(
    colData(rld)[
      ,
      c("condition", "batch"),
      drop = FALSE
    ]
  )

  heatmap_colors <- colorRampPalette(
    c("red", "white", "blue")
  )(50)

  pheatmap(
    mat,
    annotation_col = annotation_columns,
    color = heatmap_colors,
    fontsize_row = 8,
    cluster_cols = TRUE,
    filename = file.path(
      work_dir,
      "Heatmap_top30.png"
    )
  )
} else {
  warning(
    "Heatmap non generata: meno di due geni disponibili dopo il filtro."
  )
}

# ==============================================================================
# 6. VOLCANO PLOT
# ==============================================================================

cat("[INFO] Generazione Volcano Plot...\n")

vol <- res_df
vol$ID <- rownames(vol)
vol$lfc <- vol$log2FoldChange

vol$gene_type <- "Not Significant"

vol$gene_type[
  !is.na(vol$lfc) &
  vol$lfc > 0.58 &
  !is.na(vol$padj) &
  vol$padj < 0.05
] <- "Upregulated"

vol$gene_type[
  !is.na(vol$lfc) &
  vol$lfc < -0.58 &
  !is.na(vol$padj) &
  vol$padj < 0.05
] <- "Downregulated"

vol$gene_type <- factor(
  vol$gene_type,
  levels = c(
    "Upregulated",
    "Downregulated",
    "Not Significant"
  )
)

finite_lfc <- abs(
  vol$lfc[is.finite(vol$lfc)]
)

if (length(finite_lfc) == 0) {
  x_limit <- 5
} else {
  max_lfc <- max(finite_lfc)

  x_limit <- ceiling(
    min(
      max(max_lfc, 5),
      30
    )
  )
}

volcano_data <- vol %>%
  filter(
    !is.na(padj),
    !is.na(lfc),
    is.finite(lfc)
  ) %>%
  mutate(
    minus_log10_padj = -log10(
      pmax(padj, .Machine$double.xmin)
    )
  )

volcano_plot <- ggplot(
  volcano_data,
  aes(
    x = lfc,
    y = minus_log10_padj,
    colour = gene_type,
    size = gene_type,
    alpha = gene_type
  )
) +
  geom_point(shape = 16) +
  coord_cartesian(
    xlim = c(-x_limit, x_limit)
  ) +
  scale_colour_manual(
    values = c(
      "Upregulated" = "blue",
      "Downregulated" = "red",
      "Not Significant" = "grey"
    ),
    drop = FALSE
  ) +
  scale_size_manual(
    values = c(
      "Upregulated" = 2,
      "Downregulated" = 2,
      "Not Significant" = 1
    ),
    drop = FALSE
  ) +
  scale_alpha_manual(
    values = c(
      "Upregulated" = 1,
      "Downregulated" = 1,
      "Not Significant" = 0.4
    ),
    drop = FALSE
  ) +
  theme_minimal() +
  xlab("log2 Fold Change") +
  ylab("-log10 adjusted p-value") +
  ggtitle(
    paste0(
      "Volcano Plot: ",
      contrast_suffix
    )
  ) +
  theme(
    legend.title = element_blank()
  )

ggsave(
  filename = file.path(
    work_dir,
    "Volcano_plot.png"
  ),
  plot = volcano_plot,
  width = 7,
  height = 6
)

svg(
  filename = file.path(
    work_dir,
    "volcano_labelled.svg"
  ),
  width = 8,
  height = 8
)

print(volcano_plot)
dev.off()

cat(
  "[SUCCESS] Script R completato con successo. ",
  "Tutti i file sono stati salvati in: ",
  work_dir,
  "\n",
  sep = ""
)
