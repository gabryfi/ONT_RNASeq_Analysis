#!/usr/bin/env Rscript

# ==============================================================================
# DESEQ2 AUTOMATICO — VERSIONE SEMPLICE E ROBUSTA
#
# Uso:
#   Rscript deseq_auto_import_robusto.r \
#       <WORK_DIR> <TREATMENT_LABEL> <CONTROL_LABEL>
#
# Output essenziali (un errore qui interrompe lo script):
#   - deseq_results.csv
#   - samples_metadata.csv
#   - IPA_ready_input_<TRATTAMENTO>_vs_<CONTROLLO>.txt
#
# Report facoltativi (un errore qui produce un avviso ma NON interrompe DESeq2):
#   - PCA_plot.png
#   - Heatmap_top30.png
#   - Volcano_plot.png
#   - Volcano_plot.svg
#   - sessionInfo.txt
#
# Parametri opzionali mediante variabili d'ambiente:
#   DESEQ_MIN_COUNT=10
#   DESEQ_MIN_SAMPLES=LEGACY   # conserva il filtro originale rowSums >= 10
#                              # AUTO = almeno 10 conteggi nel gruppo piu' piccolo
#   DESEQ_FDR=0.05
#   DESEQ_LFC=0.5849625        # log2(1.5)
# ==============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
})

info <- function(...) {
  cat("[INFO] ", ..., "\n", sep = "")
}

warn <- function(...) {
  cat("[ATTENZIONE] ", ..., "\n", sep = "", file = stderr())
}

stop_clean <- function(...) {
  stop(paste0(...), call. = FALSE)
}

read_env_number <- function(name, default, lower = -Inf, upper = Inf) {
  raw <- Sys.getenv(name, unset = "")
  value <- if (nzchar(raw)) suppressWarnings(as.numeric(raw)) else default

  if (length(value) != 1L || !is.finite(value) || value < lower || value > upper) {
    stop_clean(
      "Valore non valido per ", name, ": '", raw,
      "'. Valore atteso tra ", lower, " e ", upper, "."
    )
  }

  value
}

atomic_write <- function(target, writer) {
  tmp <- tempfile(
    pattern = paste0(".", basename(target), "."),
    tmpdir = dirname(target)
  )

  on.exit(unlink(tmp, force = TRUE), add = TRUE)
  writer(tmp)

  if (!file.rename(tmp, target)) {
    if (!file.copy(tmp, target, overwrite = TRUE)) {
      stop_clean("Impossibile salvare il file: ", target)
    }
    unlink(tmp, force = TRUE)
  }

  invisible(target)
}

safe_report <- function(label, expression) {
  tryCatch(
    {
      force(expression)
      info(label, ": completato")
      TRUE
    },
    error = function(e) {
      warn(label, " non generato: ", conditionMessage(e))
      FALSE
    }
  )
}

# ==============================================================================
# 1. ARGOMENTI E PARAMETRI
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3L) {
  stop_clean(
    "Uso corretto: Rscript deseq_auto_import_robusto.r ",
    "<WORK_DIR> <TREATMENT_LABEL> <CONTROL_LABEL>"
  )
}

work_dir <- normalizePath(args[1], mustWork = TRUE)
treatment <- args[2]
control <- args[3]

if (!nzchar(treatment) || !nzchar(control)) {
  stop_clean("Le etichette di trattamento e controllo non possono essere vuote.")
}

if (identical(treatment, control)) {
  stop_clean("Trattamento e controllo devono essere differenti.")
}

fdr_cutoff <- read_env_number(
  "DESEQ_FDR",
  0.05,
  lower = .Machine$double.eps,
  upper = 1 - .Machine$double.eps
)
lfc_cutoff <- read_env_number("DESEQ_LFC", log2(1.5), lower = 0)
min_count <- read_env_number("DESEQ_MIN_COUNT", 10, lower = 0)

if (min_count != round(min_count)) {
  stop_clean("DESEQ_MIN_COUNT deve essere un numero intero.")
}
min_count <- as.integer(round(min_count))

info("Cartella di lavoro: ", work_dir)
info("Contrasto: ", treatment, " vs ", control)
info("FDR: ", fdr_cutoff, "; soglia |log2FC|: ", signif(lfc_cutoff, 4))

# ==============================================================================
# 2. LETTURA DEI FILE FEATURECOUNTS
# ==============================================================================

count_files <- sort(
  list.files(
    work_dir,
    pattern = "_counts\\.txt$",
    full.names = TRUE
  )
)

if (length(count_files) == 0L) {
  stop_clean("Nessun file *_counts.txt trovato in: ", work_dir)
}

extract_metadata <- function(filename) {
  base_name <- basename(filename)
  stem <- tools::file_path_sans_ext(base_name)

  # Formato: <BATCH>_barcode_<BARCODE>_<CONDITION>_counts
  # Batch e condizione possono contenere underscore.
  match_object <- regexec(
    "^(.+)_barcode_([^_]+)_(.+)_counts$",
    stem,
    perl = TRUE
  )
  fields <- regmatches(stem, match_object)[[1L]]

  if (length(fields) != 4L) {
    stop_clean(
      "Nome file non valido: ", base_name,
      ". Formato atteso: ",
      "<BATCH>_barcode_<BARCODE>_<CONDITION>_counts.txt"
    )
  }

  batch <- fields[2L]
  barcode <- fields[3L]
  condition <- fields[4L]

  data.frame(
    raw_sample = stem,
    sample = paste(batch, "bc", barcode, condition, sep = "_"),
    batch = batch,
    barcode = barcode,
    condition = condition,
    file = filename,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

sample_table <- do.call(rbind, lapply(count_files, extract_metadata))

if (anyDuplicated(sample_table$raw_sample)) {
  stop_clean("Sono presenti nomi di campione duplicati.")
}

if (anyDuplicated(sample_table$sample)) {
  stop_clean("Sono presenti identificatori di campione duplicati.")
}

read_featurecounts <- function(file_path) {
  df <- tryCatch(
    read.delim(
      file_path,
      header = TRUE,
      skip = 1,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = ""
    ),
    error = function(e) {
      stop_clean(
        "Impossibile leggere il file featureCounts ",
        basename(file_path), ": ", conditionMessage(e)
      )
    }
  )

  if (nrow(df) == 0L || ncol(df) < 2L) {
    stop_clean("File featureCounts vuoto o non valido: ", file_path)
  }

  gene_ids <- as.character(df[[1L]])
  raw_counts <- df[[ncol(df)]]
  count_values <- suppressWarnings(as.numeric(raw_counts))

  if (anyNA(gene_ids) || any(!nzchar(gene_ids))) {
    stop_clean("Gene_ID mancanti nel file: ", file_path)
  }

  if (anyDuplicated(gene_ids)) {
    stop_clean("Gene_ID duplicati nel file: ", file_path)
  }

  if (anyNA(count_values) || any(!is.finite(count_values))) {
    stop_clean("Conteggi mancanti o non numerici nel file: ", file_path)
  }

  if (any(count_values < 0)) {
    stop_clean("Conteggi negativi nel file: ", file_path)
  }

  if (any(abs(count_values - round(count_values)) > 1e-8)) {
    stop_clean("Conteggi non interi nel file: ", file_path)
  }

  stats::setNames(round(count_values), gene_ids)
}

info("File di conteggio trovati: ", length(count_files))
count_vectors <- lapply(sample_table$file, read_featurecounts)

reference_gene_ids <- names(count_vectors[[1L]])

for (i in seq_along(count_vectors)) {
  current_ids <- names(count_vectors[[i]])

  if (!setequal(current_ids, reference_gene_ids)) {
    missing_n <- sum(!reference_gene_ids %in% current_ids)
    extra_n <- sum(!current_ids %in% reference_gene_ids)

    stop_clean(
      "L'insieme dei Gene_ID non coincide nel file ",
      basename(sample_table$file[i]),
      " (mancanti: ", missing_n, "; aggiuntivi: ", extra_n, ")."
    )
  }

  # Un ordine diverso non e' un errore: viene riallineato automaticamente.
  count_vectors[[i]] <- count_vectors[[i]][reference_gene_ids]
}

count_matrix <- do.call(
  cbind,
  lapply(count_vectors, unname)
)

rownames(count_matrix) <- reference_gene_ids
colnames(count_matrix) <- sample_table$raw_sample

if (any(count_matrix > .Machine$integer.max)) {
  stop_clean("Sono presenti conteggi superiori al limite intero supportato da R.")
}

storage.mode(count_matrix) <- "integer"

# L'ordine deve essere esplicito: DESeq2 non deve indovinarlo.
sample_table <- sample_table[
  match(colnames(count_matrix), sample_table$raw_sample),
  ,
  drop = FALSE
]

if (anyNA(sample_table$raw_sample)) {
  stop_clean("Impossibile allineare conteggi e metadati dei campioni.")
}

rownames(sample_table) <- sample_table$raw_sample

if (!identical(rownames(sample_table), colnames(count_matrix))) {
  stop_clean("Conteggi e metadati non sono nello stesso ordine.")
}

available_conditions <- unique(sample_table$condition)

if (!control %in% available_conditions) {
  stop_clean("Condizione di controllo assente: ", control)
}

if (!treatment %in% available_conditions) {
  stop_clean("Condizione di trattamento assente: ", treatment)
}

contrast_samples <- table(
  factor(
    sample_table$condition[
      sample_table$condition %in% c(control, treatment)
    ],
    levels = c(control, treatment)
  )
)

if (any(contrast_samples < 2L)) {
  stop_clean(
    "DESeq2 richiede replicati biologici. Campioni disponibili: ",
    control, "=", contrast_samples[[control]], ", ",
    treatment, "=", contrast_samples[[treatment]], "."
  )
}

sample_table$batch <- droplevels(factor(sample_table$batch))
sample_table$condition <- relevel(
  droplevels(factor(sample_table$condition)),
  ref = control
)

# ==============================================================================
# 3. MODELLO, FILTRO E ANALISI DESEQ2
# ==============================================================================

if (nlevels(sample_table$batch) > 1L) {
  design_formula <- ~ batch + condition
  info("Modello: ~ batch + condition")
} else {
  design_formula <- ~ condition
  info("Modello: ~ condition")
}

col_data <- sample_table[
  ,
  c("sample", "batch", "condition"),
  drop = FALSE
]

model_matrix <- model.matrix(design_formula, data = col_data)

if (qr(model_matrix)$rank < ncol(model_matrix)) {
  batch_condition_table <- capture.output(
    print(table(sample_table$batch, sample_table$condition))
  )

  stop_clean(
    "Il modello ", deparse(design_formula),
    " non e' stimabile: batch e condizione sono confusi o mancano combinazioni.\n",
    paste(batch_condition_table, collapse = "\n"),
    "\nNon e' sicuro eliminare automaticamente il batch dal modello."
  )
}

raw_min_samples <- Sys.getenv("DESEQ_MIN_SAMPLES", unset = "LEGACY")
filter_mode <- toupper(raw_min_samples)

if (identical(filter_mode, "LEGACY")) {
  # Mantiene il comportamento dello script originale.
  min_samples <- NA_integer_
  keep <- rowSums(count_matrix) >= min_count
  info("Filtro compatibile: somma dei conteggi >= ", min_count)
} else {
  if (identical(filter_mode, "AUTO")) {
    min_samples <- as.integer(min(contrast_samples))
  } else {
    min_samples <- suppressWarnings(as.integer(raw_min_samples))

    if (
      length(min_samples) != 1L ||
      is.na(min_samples) ||
      min_samples < 1L ||
      min_samples > ncol(count_matrix)
    ) {
      stop_clean(
        "DESEQ_MIN_SAMPLES non valido: '", raw_min_samples,
        "'. Usare LEGACY, AUTO oppure un intero tra 1 e ",
        ncol(count_matrix), "."
      )
    }
  }

  keep <- rowSums(count_matrix >= min_count) >= min_samples
  info(
    "Filtro per replicati: almeno ", min_count,
    " conteggi in almeno ", min_samples, " campioni"
  )
}
info("Geni prima del filtro: ", nrow(count_matrix))
info("Geni dopo il filtro: ", sum(keep))

if (!any(keep)) {
  stop_clean(
    "Nessun gene supera il filtro. Verificare i conteggi oppure usare ",
    "DESEQ_MIN_SAMPLES=LEGACY."
  )
}

count_matrix <- count_matrix[keep, , drop = FALSE]

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = col_data,
  design = design_formula
)

info("Avvio DESeq2...")

dds <- tryCatch(
  DESeq(dds),
  error = function(e) {
    message_text <- conditionMessage(e)

    if (grepl("every gene contains at least one zero", message_text, fixed = TRUE)) {
      warn(
        "Stima standard dei size factor non possibile; ",
        "nuovo tentativo con sfType='poscounts'."
      )
      return(DESeq(dds, sfType = "poscounts"))
    }

    stop(e)
  }
)

res <- results(
  dds,
  contrast = c("condition", treatment, control),
  alpha = fdr_cutoff
)

res_df <- as.data.frame(res)
res_table <- data.frame(
  Gene_ID = rownames(res_df),
  res_df,
  row.names = NULL,
  check.names = FALSE
)

# Ordinamento pratico: risultati significativi e p-value piu' piccoli per primi.
result_order <- order(
  is.na(res_table$padj),
  res_table$padj,
  is.na(res_table$pvalue),
  res_table$pvalue
)
res_table <- res_table[result_order, , drop = FALSE]

atomic_write(
  file.path(work_dir, "deseq_results.csv"),
  function(path) {
    write.csv(res_table, file = path, row.names = FALSE, na = "")
  }
)

metadata_output <- sample_table
rownames(metadata_output) <- NULL

atomic_write(
  file.path(work_dir, "samples_metadata.csv"),
  function(path) {
    write.csv(metadata_output, file = path, row.names = FALSE, na = "")
  }
)

# ==============================================================================
# 4. FILE IPA QIAGEN
# ==============================================================================

contrast_suffix <- paste0(treatment, "_vs_", control)
lfc <- res_df$log2FoldChange

# Fold change con segno, formato comunemente usato per IPA:
# log2FC = +1 -> +2; log2FC = -1 -> -2.
signed_fold_change <- rep(NA_real_, length(lfc))
finite_lfc <- is.finite(lfc)
signed_fold_change[finite_lfc & lfc >= 0] <- 2^lfc[finite_lfc & lfc >= 0]
signed_fold_change[finite_lfc & lfc < 0] <- -2^(-lfc[finite_lfc & lfc < 0])
signed_fold_change[!is.finite(signed_fold_change)] <- NA_real_

ipa_input <- data.frame(
  Gene_ID = rownames(res_df),
  log2_fold_change = res_df$log2FoldChange,
  signed_fold_change = signed_fold_change,
  p_value = res_df$pvalue,
  FDR = res_df$padj,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

colnames(ipa_input) <- c(
  "Gene_ID",
  paste0("Expr Log2 Fold Change ", contrast_suffix),
  paste0("Expr Fold Change ", contrast_suffix),
  paste0("Expr p-value ", contrast_suffix),
  paste0("Expr FDR ", contrast_suffix)
)

atomic_write(
  file.path(work_dir, paste0("IPA_ready_input_", contrast_suffix, ".txt")),
  function(path) {
    write.table(
      ipa_input,
      file = path,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      na = ""
    )
  }
)

significant <- !is.na(res_df$padj) &
  res_df$padj < fdr_cutoff &
  !is.na(res_df$log2FoldChange) &
  abs(res_df$log2FoldChange) >= lfc_cutoff

info("Geni testati: ", nrow(res_df))
info(
  "Geni significativi (FDR < ", fdr_cutoff,
  ", |log2FC| >= ", signif(lfc_cutoff, 4), "): ",
  sum(significant)
)
info("Risultati essenziali salvati correttamente.")

# ==============================================================================
# 5. REPORT FACOLTATIVI E NON BLOCCANTI
# ==============================================================================

vsd <- tryCatch(
  {
    transformed <- tryCatch(
      vst(dds, blind = FALSE),
      error = function(e) {
        warn(
          "vst() non riuscita; uso varianceStabilizingTransformation(): ",
          conditionMessage(e)
        )
        varianceStabilizingTransformation(dds, blind = FALSE)
      }
    )
    info("Trasformazione VST: completata")
    transformed
  },
  error = function(e) {
    warn("Trasformazione VST non generata: ", conditionMessage(e))
    NULL
  }
)

if (!is.null(vsd)) {
  transformed_matrix <- assay(vsd)

  safe_report("PCA", {
    pca <- stats::prcomp(t(transformed_matrix), center = TRUE, scale. = FALSE)
    variance_percent <- 100 * pca$sdev^2 / sum(pca$sdev^2)
    condition_codes <- as.integer(sample_table$condition)
    condition_levels <- levels(sample_table$condition)

    grDevices::png(
      filename = file.path(work_dir, "PCA_plot.png"),
      width = 1400,
      height = 1100,
      res = 180
    )

    tryCatch(
      {
        graphics::plot(
          pca$x[, 1L],
          pca$x[, 2L],
          pch = 19,
          col = condition_codes,
          xlab = paste0("PC1: ", round(variance_percent[1L]), "%"),
          ylab = paste0("PC2: ", round(variance_percent[2L]), "%"),
          main = "PCA dei campioni"
        )
        graphics::text(
          pca$x[, 1L],
          pca$x[, 2L],
          labels = sample_table$sample,
          pos = 3,
          cex = 0.65
        )
        graphics::legend(
          "topright",
          legend = condition_levels,
          col = seq_along(condition_levels),
          pch = 19,
          bty = "n"
        )
      },
      finally = grDevices::dev.off()
    )
  })

  safe_report("Heatmap top 30", {
    if (!requireNamespace("pheatmap", quietly = TRUE)) {
      stop("pacchetto opzionale 'pheatmap' non installato")
    }

    variance_values <- apply(transformed_matrix, 1L, stats::var)
    variance_values[!is.finite(variance_values)] <- 0
    number_top_genes <- min(30L, length(variance_values))

    if (number_top_genes < 2L) {
      stop("meno di due geni disponibili")
    }

    top_genes <- head(order(variance_values, decreasing = TRUE), number_top_genes)
    heatmap_matrix <- transformed_matrix[top_genes, , drop = FALSE]
    heatmap_matrix <- heatmap_matrix - rowMeans(heatmap_matrix)

    annotation_columns <- data.frame(
      condition = sample_table$condition,
      batch = sample_table$batch,
      row.names = colnames(heatmap_matrix),
      check.names = FALSE
    )

    pheatmap::pheatmap(
      heatmap_matrix,
      annotation_col = annotation_columns,
      color = grDevices::colorRampPalette(c("blue", "white", "red"))(50),
      fontsize_row = 8,
      cluster_cols = TRUE,
      filename = file.path(work_dir, "Heatmap_top30.png")
    )
  })
}

draw_volcano <- function() {
  valid <- !is.na(res_df$padj) &
    is.finite(res_df$padj) &
    !is.na(res_df$log2FoldChange) &
    is.finite(res_df$log2FoldChange)

  if (!any(valid)) {
    stop("nessun risultato finito disponibile")
  }

  x <- res_df$log2FoldChange[valid]
  y <- -log10(pmax(res_df$padj[valid], .Machine$double.xmin))

  classes <- rep("Non significativo", length(x))
  classes[x >= lfc_cutoff & res_df$padj[valid] < fdr_cutoff] <- "Upregulated"
  classes[x <= -lfc_cutoff & res_df$padj[valid] < fdr_cutoff] <- "Downregulated"

  point_colors <- c(
    "Upregulated" = "blue",
    "Downregulated" = "red",
    "Non significativo" = "grey"
  )

  finite_abs_lfc <- abs(x[is.finite(x)])
  x_limit <- if (length(finite_abs_lfc) == 0L) {
    5
  } else {
    ceiling(min(max(max(finite_abs_lfc), 5), 30))
  }

  graphics::plot(
    x,
    y,
    pch = 16,
    cex = ifelse(classes == "Non significativo", 0.55, 0.8),
    col = unname(point_colors[classes]),
    xlim = c(-x_limit, x_limit),
    xlab = "log2 Fold Change",
    ylab = "-log10 adjusted p-value",
    main = paste0("Volcano Plot: ", contrast_suffix)
  )

  graphics::abline(
    v = c(-lfc_cutoff, lfc_cutoff),
    h = -log10(fdr_cutoff),
    lty = 2,
    col = "black"
  )

  graphics::legend(
    "topright",
    legend = names(point_colors),
    col = unname(point_colors),
    pch = 16,
    bty = "n"
  )
}

safe_report("Volcano PNG", {
  grDevices::png(
    filename = file.path(work_dir, "Volcano_plot.png"),
    width = 1400,
    height = 1200,
    res = 180
  )
  tryCatch(draw_volcano(), finally = grDevices::dev.off())
})

safe_report("Volcano SVG", {
  grDevices::svg(
    filename = file.path(work_dir, "Volcano_plot.svg"),
    width = 8,
    height = 7
  )
  tryCatch(draw_volcano(), finally = grDevices::dev.off())
})

safe_report("Session info", {
  atomic_write(
    file.path(work_dir, "sessionInfo.txt"),
    function(path) {
      writeLines(capture.output(sessionInfo()), con = path)
    }
  )
})

cat(
  "\n[SUCCESS] Analisi DESeq2 completata.\n",
  "Output: ", work_dir, "\n",
  sep = ""
)
