#' @importFrom data.table data.table as.data.table rbindlist setnames fread :=
NULL

# ---------------------------------------------------------------------------
# Low-level REST client
# ---------------------------------------------------------------------------

ens_host <- function(human_assembly = c("GRCh38", "GRCh37")) {
  human_assembly <- match.arg(human_assembly)
  if (human_assembly == "GRCh37") "https://grch37.rest.ensembl.org" else "https://rest.ensembl.org"
}

.ens_get_raw <- function(path, query = list(), host = ens_host()) {
  req <- httr2::request(host)
  req <- httr2::req_url_path_append(req, path)
  req <- httr2::req_headers(req, Accept = "application/json")
  if (length(query)) req <- httr2::req_url_query(req, !!!query)
  req <- httr2::req_retry(req, max_tries = 4, backoff = ~ 2)
  req <- httr2::req_user_agent(req, "xsmut (https://github.com/your-org/xsmut)")
  resp <- httr2::req_perform(req)
  jsonlite::fromJSON(httr2::resp_body_string(resp), simplifyVector = FALSE)
}

# Cache identical calls for the R session (Ensembl asks clients to be polite).
ens_get <- memoise::memoise(.ens_get_raw)

# ---------------------------------------------------------------------------
# Species & assembly
# ---------------------------------------------------------------------------

#' Resolve a user-supplied species name to an Ensembl production name
#'
#' Accepts scientific names ("Macaca mulatta"), Ensembl names ("macaca_mulatta"),
#' common names ("macaque", "rhesus") or Ensembl aliases.
#'
#' @param species character(1)
#' @param host Ensembl REST host
#' @return character(1) Ensembl production name, e.g. "macaca_mulatta"
#' @export
ens_resolve_species <- function(species, host = ens_host()) {
  s <- tolower(trimws(species))
  s_us <- gsub("[ .-]+", "_", s)
  info <- ens_get("info/species", host = host)$species
  for (sp in info) {
    cands <- tolower(c(sp$name, sp$common_name, sp$display_name, unlist(sp$aliases)))
    if (s %in% cands || s_us %in% cands) return(sp$name)
  }
  # partial match on common name (e.g. "macaque" -> "Macaque" alias for macaca_mulatta)
  for (sp in info) {
    cands <- tolower(c(sp$common_name, sp$display_name, unlist(sp$aliases)))
    if (any(grepl(s, cands, fixed = TRUE))) return(sp$name)
  }
  stop("Could not resolve species '", species, "' in Ensembl. Try the scientific name.")
}

#' Get the current Ensembl assembly for a species and optionally check it
#'
#' @param species Ensembl production name (see [ens_resolve_species()])
#' @param assembly optional expected assembly name (e.g. "Mmul_10"); a warning is
#'   issued if it does not match what Ensembl serves.
#' @param host Ensembl REST host
#' @return list with `assembly_name`, `assembly_accession`, `matches`
#' @export
ens_assembly <- function(species, assembly = NULL, host = ens_host()) {
  a <- ens_get(paste0("info/assembly/", species), host = host)
  out <- list(assembly_name = a$assembly_name, assembly_accession = a$assembly_accession,
              default_coord_system = a$default_coord_system_version, matches = NA)
  if (!is.null(assembly)) {
    norm <- function(x) tolower(gsub("[^a-z0-9]", "", tolower(x)))
    out$matches <- norm(assembly) %in% norm(c(a$assembly_name, a$default_coord_system_version, a$assembly_accession))
    if (!out$matches) {
      warning(sprintf("Requested assembly '%s' but Ensembl REST serves '%s' for %s. ",
                      assembly, a$assembly_name, species),
              "Coordinates in your mutation table must match the served assembly ",
              "(for human GRCh37 use human_assembly = 'GRCh37').")
    }
  }
  out
}

# ---------------------------------------------------------------------------
# Gene lookup
# ---------------------------------------------------------------------------

ens_lookup_symbol <- function(symbol, species, host = ens_host()) {
  ens_get(paste0("lookup/symbol/", species, "/", symbol),
          query = list(expand = 1, `content-type` = "application/json"), host = host)
}

ens_lookup_id <- function(id, host = ens_host(), expand = TRUE) {
  ens_get(paste0("lookup/id/", id), query = list(expand = as.integer(expand)), host = host)
}

# ---------------------------------------------------------------------------
# Orthology & sequences
# ---------------------------------------------------------------------------

#' Find the ortholog of a human gene in another species
#'
#' @param gene human gene symbol
#' @param target_species Ensembl production name
#' @param host Ensembl REST host
#' @return list with `human` (id, protein_id) and `target` (id, protein_id, perc_id)
#' @export
get_ortholog <- function(gene, target_species, host = ens_host()) {
  res <- ens_get(paste0("homology/symbol/human/", gene),
                 query = list(target_species = target_species, type = "orthologues",
                              sequence = "protein"), host = host)
  homs <- res$data[[1]]$homologies
  if (!length(homs)) stop("No Ensembl ortholog of ", gene, " found in ", target_species)
  # prefer one2one; else highest %id
  type <- vapply(homs, function(h) h$type, "")
  pid  <- vapply(homs, function(h) as.numeric(h$target$perc_id), 0)
  ord  <- order(type != "ortholog_one2one", -pid)
  h <- homs[[ord[1]]]
  list(
    human  = list(id = h$source$id, protein_id = h$source$protein_id,
                  seq = h$source$align_seq),
    target = list(id = h$target$id, protein_id = h$target$protein_id,
                  seq = h$target$align_seq, perc_id = h$target$perc_id,
                  species = h$target$species, type = h$type)
  )
}

ens_protein_seq <- function(protein_id, host = ens_host()) {
  r <- ens_get(paste0("sequence/id/", protein_id), query = list(type = "protein"), host = host)
  r$seq
}

# Map protein residues -> genomic coordinates (Ensembl /map/translation)
ens_map_translation <- function(protein_id, start, end, host = ens_host()) {
  r <- ens_get(paste0("map/translation/", protein_id, "/", start, "..", end), host = host)
  m <- r$mappings
  if (!length(m)) return(data.table(chrom = character(), start = integer(), end = integer(), strand = integer()))
  rbindlist(lapply(m, function(x) data.table(chrom = x$seq_region_name, start = x$start,
                                              end = x$end, strand = x$strand)))
}
