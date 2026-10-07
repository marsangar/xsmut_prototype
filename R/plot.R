#' @import ggplot2
NULL

# Piecewise-linear transform that keeps exons at native scale and squeezes each
# intron to a fixed width (fraction of total exon length). Returns a function.
.make_transform <- function(exons, compress_introns = TRUE, intron_frac = 0.04) {
  ex <- as.data.table(exons)[, .(start = min(start), end = max(end)), by = exon_rank][order(start)]
  if (!compress_introns || nrow(ex) < 2) return(list(f = identity, exons = ex))
  intron_w <- intron_frac * sum(ex$end - ex$start + 1)
  # cumulative plot offsets
  ex[, plot_start := 0]
  for (i in seq_len(nrow(ex))) {
    ex$plot_start[i] <- if (i == 1) ex$start[1] else ex$plot_start[i - 1] + (ex$end[i - 1] - ex$start[i - 1] + 1) + intron_w
  }
  ex[, plot_end := plot_start + (end - start)]
  f <- function(x) {
    out <- rep(NA_real_, length(x))
    for (i in seq_len(nrow(ex))) {
      in_ex <- !is.na(x) & x >= ex$start[i] & x <= ex$end[i]
      out[in_ex] <- ex$plot_start[i] + (x[in_ex] - ex$start[i])
      if (i < nrow(ex)) {
        in_int <- !is.na(x) & x > ex$end[i] & x < ex$start[i + 1]
        frac <- (x[in_int] - ex$end[i]) / (ex$start[i + 1] - ex$end[i])
        out[in_int] <- ex$plot_end[i] + frac * intron_w
      }
    }
    # outside gene: clamp linearly beyond the ends
    lo <- !is.na(x) & x < ex$start[1]; out[lo] <- ex$plot_start[1] - (ex$start[1] - x[lo])
    hi <- !is.na(x) & x > ex$end[nrow(ex)]; out[hi] <- ex$plot_end[nrow(ex)] + (x[hi] - ex$end[nrow(ex)])
    out
  }
  list(f = f, exons = ex)
}

.consequence_class <- function(x) {
  x <- tolower(as.character(x))
  data.table::fcase(
    grepl("missense", x), "Missense",
    grepl("nonsense|stop_gained|stop gained", x), "Nonsense",
    grepl("frameshift", x), "Frameshift",
    grepl("inframe|in frame|deletion|insertion", x), "In-frame indel",
    grepl("splice", x), "Splice",
    grepl("synonymous|silent", x), "Synonymous",
    is.na(x) | x == "", "Unknown",
    default = "Other"
  )
}

.pal <- c(Missense = "#2E86AB", Nonsense = "#C73E1D", Frameshift = "#F18F01",
          "In-frame indel" = "#7B2D8E", Splice = "#3B8B5A", Synonymous = "#9E9E9E",
          Other = "#5E5E5E", Unknown = "#BDBDBD")

.lollipop <- function(dt, tf, title, flip = FALSE, min_label_n = 3) {
  # dt: pos, aa_change, consequence, n
  dt <- as.data.table(dt)
  if (!nrow(dt)) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "no mappable mutations", colour = "grey50") +
             theme_void() + labs(title = title))
  }
  dt[, class := .consequence_class(consequence)]
  dt[, x := tf(pos)]
  dt[, ymax := if (flip) -n else n]
  lab <- dt[n >= min_label_n]
  p <- ggplot(dt) +
    geom_segment(aes(x = x, xend = x, y = 0, yend = ymax), colour = "grey60", linewidth = 0.3) +
    geom_point(aes(x = x, y = ymax, colour = class, size = n), alpha = 0.85) +
    scale_colour_manual(values = .pal, drop = TRUE, name = NULL) +
    scale_size_area(max_size = 6, guide = "none") +
    labs(title = title, y = "samples", x = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "right", plot.title = element_text(face = "bold"))
  if (nrow(lab)) {
    p <- p + geom_text(data = lab, aes(x = x, y = ymax, label = aa_change),
                       vjust = if (flip) 1.6 else -0.8, size = 2.7, check_overlap = TRUE)
  }
  if (flip) p <- p + scale_y_continuous(labels = function(v) abs(v))
  p
}

.gene_track <- function(model, tf, show_exon_numbers = TRUE) {
  ex <- as.data.table(model$exons)
  ex[, `:=`(px = tf(start), pxend = tf(end))]
  ex[, h := ifelse(feature == "CDS", 0.5, 0.25)]
  merged <- ex[, .(px = min(px), pxend = max(pxend), start = min(start), end = max(end)), by = exon_rank][order(px)]
  strand_lab <- if (model$gene$strand == 1) "5' \u2192 3' (+)" else "3' \u2190 5' (\u2212)"
  p <- ggplot() +
    geom_segment(aes(x = min(merged$px), xend = max(merged$pxend), y = 0, yend = 0), colour = "grey40") +
    geom_rect(data = ex, aes(xmin = px, xmax = pxend, ymin = -h, ymax = h, fill = feature), colour = NA) +
    scale_fill_manual(values = c(CDS = "#1F3A5F", UTR = "#9DB4D0"), name = NULL) +
    coord_cartesian(ylim = c(-1, 1)) +
    labs(x = sprintf("%s  chr%s:%s-%s (%s, %s)  %s", model$gene$symbol, model$gene$chrom,
                     format(model$gene$start, big.mark = ","), format(model$gene$end, big.mark = ","),
                     model$transcript$id, model$gene$assembly, strand_lab), y = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
          legend.position = "right")
  if (show_exon_numbers) {
    p <- p + geom_text(data = merged, aes(x = (px + pxend) / 2, y = -0.75, label = exon_rank), size = 2.6, colour = "grey30")
  }
  p
}

#' Plot human COSMIC mutations (top) vs species mutations (bottom) on a gene model
#'
#' @param model `gene_model` from [get_gene_model()] (human)
#' @param cosmic data.table from [cosmic_gene_mutations()] (needs `start`,
#'   `aa_change`, `consequence`)
#' @param species_mut data.table from [map_mutations_to_human()] (needs
#'   `human_pos`; optional `aa_change`, `consequence`)
#' @param species_label label for the lower panel, e.g. "Macaque (Mmul_10)"
#' @param compress_introns squeeze introns to a fixed width so exons dominate
#' @param coding_only drop COSMIC rows without a protein change
#' @param min_label_n label mutations recurring in at least this many samples
#' @return a patchwork object (ggplot); use `print()` or `ggsave()`
#' @export
plot_gene_mutations <- function(model, cosmic, species_mut, species_label = "Species",
                                compress_introns = TRUE, coding_only = TRUE, min_label_n = 3) {
  tr <- .make_transform(model$exons, compress_introns)
  tf <- tr$f

  cs <- as.data.table(cosmic)[!is.na(start)]
  if (coding_only) cs <- cs[!is.na(aa_pos)]
  cs_agg <- cs[, .(n = .N), by = .(pos = start, aa_change, consequence)]

  sp <- as.data.table(species_mut)[!is.na(human_pos)]
  if (!"aa_change" %in% names(sp)) sp[, aa_change := paste0("p.", aa_pos)]
  if (!"consequence" %in% names(sp)) sp[, consequence := NA_character_]
  sp_agg <- sp[, .(n = .N), by = .(pos = human_pos, aa_change, consequence)]

  n_unmapped <- sum(is.na(as.data.table(species_mut)$human_pos))
  bottom_title <- sprintf("%s mutations (n = %d mapped%s)", species_label, nrow(sp),
                          if (n_unmapped) sprintf(", %d unmapped", n_unmapped) else "")

  top <- .lollipop(cs_agg, tf, sprintf("Human COSMIC: %s (%d mutations, %d samples)",
                                       model$gene$symbol, nrow(cs_agg), nrow(cs)), min_label_n = min_label_n)
  mid <- .gene_track(model, tf)
  bot <- .lollipop(sp_agg, tf, bottom_title, flip = TRUE, min_label_n = 1)

  xr <- range(c(tf(model$exons$start), tf(model$exons$end), tf(cs_agg$pos), tf(sp_agg$pos)), na.rm = TRUE)
  pad <- diff(xr) * 0.02
  lock <- function(p) p + coord_cartesian(xlim = xr + c(-pad, pad), clip = "off")
  top <- lock(top); bot <- lock(bot)
  mid <- mid + coord_cartesian(xlim = xr + c(-pad, pad), ylim = c(-1, 1))

  patchwork::wrap_plots(top, mid, bot, ncol = 1, heights = c(4, 1, 3)) +
    patchwork::plot_layout(guides = "collect")
}
