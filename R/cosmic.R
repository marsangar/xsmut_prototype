# COSMIC has no public API: data must be downloaded (free for academic use) from
# https://cancer.sanger.ac.uk/cosmic/download. This module reads those files.
#
# Supported inputs
#  * Mutations:  Cosmic_MutantCensus_v*_GRCh38.tsv(.gz)   (current format, v98+)
#                Cosmic_GenomeScreensMutant_v*_GRCh38.tsv(.gz)
#                CosmicMutantExport.tsv(.gz)                (legacy format, <= v97)
#  * Gene list:  Cosmic_CancerGeneCensus_v*_GRCh38.tsv(.gz) or cancer_gene_census.csv
#
# Everything is normalised to a common schema so downstream code is format-agnostic.

.cosmic_schema <- c(
  # new format                 # legacy format
  GENE_SYMBOL              = "Gene name",
  TRANSCRIPT_ACCESSION     = "Accession Number",
  MUTATION_CDS             = "Mutation CDS",
  MUTATION_AA              = "Mutation AA",
  MUTATION_DESCRIPTION     = "Mutation Description",
  GENOMIC_MUTATION_ID      = "GENOMIC_MUTATION_ID",
  COSMIC_SAMPLE_ID         = "ID_sample",
  PRIMARY_SITE             = "Primary site",
  PRIMARY_HISTOLOGY        = "Primary histology",
  MUTATION_SOMATIC_STATUS  = "Mutation somatic status",
  CHROMOSOME               = NA,   # legacy: parsed from "Mutation genome position"
  GENOME_START             = NA,
  GENOME_STOP              = NA,
  STRAND                   = "Mutation strand",
  GENOMIC_WT_ALLELE        = NA,
  GENOMIC_MUT_ALLELE       = NA
)

.cosmic_std_names <- c(
  GENE_SYMBOL = "gene", TRANSCRIPT_ACCESSION = "transcript", MUTATION_CDS = "cds_change",
  MUTATION_AA = "aa_change", MUTATION_DESCRIPTION = "consequence",
  GENOMIC_MUTATION_ID = "mutation_id", COSMIC_SAMPLE_ID = "sample_id",
  PRIMARY_SITE = "primary_site", PRIMARY_HISTOLOGY = "histology",
  MUTATION_SOMATIC_STATUS = "somatic_status", CHROMOSOME = "chrom",
  GENOME_START = "start", GENOME_STOP = "end", STRAND = "strand",
  GENOMIC_WT_ALLELE = "ref", GENOMIC_MUT_ALLELE = "alt"
)

#' Load and normalise a COSMIC mutation export
#'
#' @param path path to a COSMIC mutation TSV (optionally gzipped)
#' @param genes optional character vector; only rows for these genes are kept
#'   (much faster / lower memory for the full genome-screens file)
#' @param somatic_only keep only confirmed/reported somatic variants
#' @return data.table with standard columns: gene, transcript, cds_change,
#'   aa_change, aa_pos, consequence, mutation_id, sample_id, primary_site,
#'   histology, chrom, start, end, ref, alt
#' @export
cosmic_load <- function(path, genes = NULL, somatic_only = TRUE) {
  stopifnot(file.exists(path))
  hdr <- names(fread(path, nrows = 0))
  legacy <- "Gene name" %in% hdr

  if (!is.null(genes)) {
    gene_col <- if (legacy) "Gene name" else "GENE_SYMBOL"
    # Stream with grep pre-filter when possible (fast on the multi-GB file)
    pat <- paste(genes, collapse = "|")
    cmd <- sprintf("%s %s | grep -E -w '%s'",
                   if (grepl("\\.gz$", path)) "gzip -cd" else "cat", shQuote(path), pat)
    dt <- tryCatch(fread(cmd = cmd, header = FALSE, col.names = hdr, quote = ""),
                   error = function(e) NULL)
    if (is.null(dt)) dt <- fread(path, quote = "")
    dt <- dt[get(gene_col) %in% genes]
  } else {
    dt <- fread(path, quote = "")
  }

  if (legacy) {
    # rename legacy -> new names
    for (new in names(.cosmic_schema)) {
      old <- .cosmic_schema[[new]]
      if (!is.na(old) && old %in% names(dt)) setnames(dt, old, new)
    }
    # "Mutation genome position" = "17:7577120-7577120"
    if ("Mutation genome position" %in% names(dt)) {
      pos <- dt[["Mutation genome position"]]
      dt[, CHROMOSOME := sub(":.*", "", pos)]
      dt[, GENOME_START := as.integer(sub("^[^:]+:(\\d+)-.*", "\\1", pos))]
      dt[, GENOME_STOP  := as.integer(sub(".*-(\\d+)$", "\\1", pos))]
    }
    dt[, GENOMIC_WT_ALLELE := NA_character_]
    dt[, GENOMIC_MUT_ALLELE := NA_character_]
  }

  keep <- intersect(names(.cosmic_std_names), names(dt))
  dt <- dt[, ..keep]
  setnames(dt, keep, .cosmic_std_names[keep])
  for (col in setdiff(.cosmic_std_names, names(dt))) dt[, (col) := NA]

  if (somatic_only && "somatic_status" %in% names(dt)) {
    dt <- dt[grepl("somatic", tolower(somatic_status)) & !grepl("not", tolower(somatic_status))]
  }
  dt[, aa_pos := parse_aa_pos(aa_change)]
  dt[, gene := sub("_ENST.*", "", gene)]  # legacy files use GENE_ENSTxxx for alt transcripts
  dt[]
}

#' Parse the (first) amino-acid position from an HGVS protein string
#' e.g. "p.R175H" -> 175, "p.Q331*" -> 331, "p.E285_K286del" -> 285, "p.?" -> NA
#' @param x character
#' @return integer
#' @export
parse_aa_pos <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- sub("^p\\.\\(?", "", x)
  r <- regexpr("[A-Z][a-z]{0,2}(\\d+)", x)
  out <- rep(NA_integer_, length(x))
  hit <- r > 0
  if (any(hit)) out[hit] <- as.integer(gsub("\\D", "", regmatches(x, r)))
  out
}

#' Check whether a gene is listed in COSMIC (Cancer Gene Census)
#'
#' @param genes character vector of symbols
#' @param census_path path to the Cancer Gene Census file (tsv/csv). If NULL and
#'   `mutations` is given, the check falls back to "has >= 1 COSMIC mutation".
#' @param mutations optional normalised COSMIC mutation table from [cosmic_load()]
#' @return data.table: gene, in_census, tier, role, n_cosmic_mutations
#' @export
cosmic_gene_listed <- function(genes, census_path = NULL, mutations = NULL) {
  out <- data.table(gene = genes, in_census = NA, tier = NA_character_,
                    role = NA_character_, n_cosmic_mutations = NA_integer_)
  if (!is.null(census_path)) {
    cg <- fread(census_path)
    sym_col  <- intersect(c("GENE_SYMBOL", "Gene Symbol"), names(cg))[1]
    tier_col <- intersect(c("TIER", "Tier"), names(cg))[1]
    role_col <- intersect(c("ROLE_IN_CANCER", "Role in Cancer"), names(cg))[1]
    if (is.na(sym_col)) stop("Census file has no GENE_SYMBOL / 'Gene Symbol' column")
    idx <- match(toupper(genes), toupper(cg[[sym_col]]))
    # compute outside the data.table call to avoid column/variable name clashes
    tier_val <- if (!is.na(tier_col)) as.character(cg[[tier_col]][idx]) else rep(NA_character_, length(genes))
    role_val <- if (!is.na(role_col)) as.character(cg[[role_col]][idx]) else rep(NA_character_, length(genes))
    data.table::set(out, j = "in_census", value = !is.na(idx))
    data.table::set(out, j = "tier", value = tier_val)
    data.table::set(out, j = "role", value = role_val)
  }
  if (!is.null(mutations)) {
    n <- mutations[, .N, by = gene]
    n_val <- n$N[match(toupper(genes), toupper(n$gene))]
    n_val[is.na(n_val)] <- 0L
    data.table::set(out, j = "n_cosmic_mutations", value = as.integer(n_val))
    if (is.null(census_path)) data.table::set(out, j = "in_census", value = n_val > 0)
  }
  out[]
}

#' Retrieve all human COSMIC mutations for a gene
#'
#' @param gene symbol
#' @param cosmic normalised table from [cosmic_load()] or a file path
#' @param coding_only drop synonymous / unknown / non-coding entries
#' @return data.table (same schema as [cosmic_load()]) plus `n_samples` per
#'   unique mutation available via `attr(x, "recurrence")`
#' @export
cosmic_gene_mutations <- function(gene, cosmic, coding_only = FALSE) {
  if (is.character(cosmic)) cosmic <- cosmic_load(cosmic, genes = gene)
  g <- toupper(gene)
  m <- cosmic[toupper(cosmic$gene) == g]
  if (coding_only) {
    m <- m[!grepl("synonymous|coding silent|unknown|intronic", tolower(consequence)) &
             !is.na(aa_pos)]
  }
  rec <- m[!is.na(start), .(n_samples = .N), by = .(chrom, start, end, aa_change, consequence)]
  attr(m, "recurrence") <- rec[order(-n_samples)]
  m[]
}
