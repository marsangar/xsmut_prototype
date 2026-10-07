#' Retrieve a gene model (canonical transcript exons/CDS) from Ensembl
#'
#' @param gene gene symbol (e.g. "TP53")
#' @param species Ensembl production name; default human
#' @param transcript_id optional Ensembl transcript ID to use instead of the
#'   Ensembl canonical transcript
#' @param human_assembly "GRCh38" (default) or "GRCh37"; only affects the REST host
#' @return list of class `gene_model` with `gene` (list), `transcript` (list) and
#'   `exons` (data.table: chrom, start, end, strand, exon_rank, feature in
#'   c("UTR","CDS"))
#' @export
get_gene_model <- function(gene, species = "homo_sapiens", transcript_id = NULL,
                           human_assembly = c("GRCh38", "GRCh37")) {
  host <- ens_host(match.arg(human_assembly))
  g <- ens_lookup_symbol(gene, species, host = host)
  tx <- g$Transcript
  if (!length(tx)) stop("No transcripts returned for ", gene, " in ", species)

  if (!is.null(transcript_id)) {
    t <- Filter(function(x) x$id == transcript_id, tx)
    if (!length(t)) stop("Transcript ", transcript_id, " not found for ", gene)
    t <- t[[1]]
  } else {
    canon <- Filter(function(x) isTRUE(as.logical(x$is_canonical)), tx)
    if (length(canon)) {
      t <- canon[[1]]
    } else {
      # fall back: protein-coding transcript with the longest CDS
      pc <- Filter(function(x) !is.null(x$Translation), tx)
      if (!length(pc)) pc <- tx
      len <- vapply(pc, function(x) if (!is.null(x$Translation)) x$Translation$length else 0, 0)
      t <- pc[[which.max(len)]]
    }
  }

  ex <- rbindlist(lapply(t$Exon, function(e)
    data.table(chrom = e$seq_region_name, start = e$start, end = e$end, strand = e$strand)))
  ex <- ex[order(start)]
  ex[, exon_rank := if (t$strand == 1) seq_len(.N) else rev(seq_len(.N))]

  # Split each exon into UTR / CDS segments using the translation's genomic bounds
  if (!is.null(t$Translation)) {
    cds_start <- min(t$Translation$start, t$Translation$end)
    cds_end   <- max(t$Translation$start, t$Translation$end)
    segs <- rbindlist(lapply(seq_len(nrow(ex)), function(i) {
      s <- ex$start[i]; e <- ex$end[i]
      out <- list()
      if (e < cds_start || s > cds_end) {
        out[[1]] <- data.table(start = s, end = e, feature = "UTR")
      } else {
        if (s < cds_start) out[[length(out) + 1]] <- data.table(start = s, end = cds_start - 1, feature = "UTR")
        out[[length(out) + 1]] <- data.table(start = max(s, cds_start), end = min(e, cds_end), feature = "CDS")
        if (e > cds_end) out[[length(out) + 1]] <- data.table(start = cds_end + 1, end = e, feature = "UTR")
      }
      r <- rbindlist(out)
      r[, `:=`(chrom = ex$chrom[i], strand = ex$strand[i], exon_rank = ex$exon_rank[i])]
      r
    }))
    ex <- segs[, .(chrom, start, end, strand, exon_rank, feature)]
  } else {
    ex[, feature := "UTR"]
  }

  structure(list(
    gene = list(symbol = g$display_name, id = g$id, species = species,
                chrom = g$seq_region_name, start = g$start, end = g$end,
                strand = g$strand, biotype = g$biotype, assembly = g$assembly_name),
    transcript = list(id = t$id, display_name = t$display_name, biotype = t$biotype,
                      is_canonical = isTRUE(as.logical(t$is_canonical)),
                      protein_id = if (!is.null(t$Translation)) t$Translation$id else NA_character_,
                      protein_length = if (!is.null(t$Translation)) t$Translation$length else NA_integer_),
    exons = ex,
    host = host
  ), class = "gene_model")
}

#' @export
print.gene_model <- function(x, ...) {
  cat(sprintf("<gene_model> %s (%s) %s chr%s:%d-%d strand %d\n",
              x$gene$symbol, x$gene$id, x$gene$species, x$gene$chrom,
              x$gene$start, x$gene$end, x$gene$strand))
  cat(sprintf("  transcript %s (%s%s), protein %s (%s aa), %d exon segments\n",
              x$transcript$id, x$transcript$biotype,
              if (x$transcript$is_canonical) ", canonical" else "",
              x$transcript$protein_id, x$transcript$protein_length, nrow(x$exons)))
  invisible(x)
}
