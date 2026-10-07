#' Run the full cross-species comparison for one gene
#'
#' @param gene human gene symbol
#' @param mutations data.frame of species mutations (see [map_mutations_to_human()])
#' @param species species name in any form accepted by [ens_resolve_species()]
#' @param assembly expected species assembly (e.g. "Mmul_10"); checked against Ensembl
#' @param cosmic COSMIC mutation table from [cosmic_load()] or a file path
#' @param census_path optional Cancer Gene Census file
#' @param human_assembly "GRCh38" or "GRCh37" (must match your COSMIC download)
#' @param ... passed to [plot_gene_mutations()]
#' @return list: gene_check, model, cosmic_mutations, species_mutations, plot
#' @export
xsmut_pipeline <- function(gene, mutations, species, assembly = NULL, cosmic,
                           census_path = NULL, human_assembly = c("GRCh38", "GRCh37"), ...) {
  human_assembly <- match.arg(human_assembly)
  host <- ens_host(human_assembly)
  sp <- ens_resolve_species(species, host = host)
  asm <- ens_assembly(sp, assembly, host = host)

  if (is.character(cosmic)) cosmic <- cosmic_load(cosmic, genes = gene)
  cos_mut <- cosmic_gene_mutations(gene, cosmic)
  check <- cosmic_gene_listed(gene, census_path = census_path, mutations = cosmic)
  if (!isTRUE(check$in_census) && check$n_cosmic_mutations == 0)
    warning(gene, " has no COSMIC entries in the supplied files")

  model <- get_gene_model(gene, species = "homo_sapiens", human_assembly = human_assembly)
  sp_mut <- map_mutations_to_human(mutations, gene, sp, human_assembly = human_assembly)

  label <- sprintf("%s (%s)", tools::toTitleCase(gsub("_", " ", sp)),
                   if (!is.null(assembly)) assembly else asm$assembly_name)
  p <- plot_gene_mutations(model, cos_mut, sp_mut, species_label = label, ...)

  list(gene_check = check, species = sp, assembly = asm, model = model,
       cosmic_mutations = cos_mut, species_mutations = sp_mut, plot = p)
}

#' Launch the Shiny app
#' @param ... passed to [shiny::runApp()]
#' @export
run_app <- function(...) {
  shiny::runApp(system.file("app", package = "xsmut"), ...)
}
