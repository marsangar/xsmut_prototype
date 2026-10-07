# xsmut — cross-species mutations vs. human COSMIC

`xsmut` takes somatic mutations annotated in **any Ensembl species** (macaque, dog, mouse, …),
maps them onto **human reference coordinates**, and plots them beneath the human COSMIC
mutation landscape for the same gene.

```
 ┌──────────────────────────────────────────────┐
 │  Human COSMIC lollipops (n samples)          │   upper panel
 ├──────────────────────────────────────────────┤
 │  ▬▬█████▬▬▬████▬▬▬▬███▬▬  gene model (CDS/UTR)│   shared axis
 ├──────────────────────────────────────────────┤
 │  Species lollipops, mapped to human          │   lower panel
 └──────────────────────────────────────────────┘
```

## How the species-agnostic mapping works

1. `ens_resolve_species("macaque")` → `macaca_mulatta`; `ens_assembly()` checks the requested
   assembly (e.g. `Mmul_10`) against what Ensembl serves.
2. `get_ortholog("TP53", "macaca_mulatta")` fetches the human↔species ortholog pair **with its
   protein alignment** from Ensembl Compara.
3. Each species amino-acid position is translated to the aligned human residue, then to a human
   genomic coordinate with Ensembl `/map/translation`.
4. Optionally, rows without a protein position (non-coding) are lifted with a UCSC chain file
   via `rtracklayer::liftOver`.

Nothing in steps 1–4 is macaque-specific: swap the species string and it works for any genome in
Ensembl. If your mutations were called against a specific species transcript, pass its protein ID
(`species_protein_id = "ENSMMUP..."`) and a fresh pairwise alignment is used instead.

## COSMIC data

COSMIC has **no public API**; downloads require a (free academic) account at
<https://cancer.sanger.ac.uk/cosmic/download>. Place the files in `data/cosmic/` (git-ignored):

| purpose            | file (current format)                        | legacy format             |
|--------------------|----------------------------------------------|---------------------------|
| mutations          | `Cosmic_MutantCensus_v*_GRCh38.tsv.gz`       | `CosmicMutantExport.tsv`  |
| "is gene listed?"  | `Cosmic_CancerGeneCensus_v*_GRCh38.tsv.gz`   | `cancer_gene_census.csv`  |

Both formats are normalised by `cosmic_load()`. Use the GRCh38 files with the default
`human_assembly = "GRCh38"`, or GRCh37 files with `human_assembly = "GRCh37"`.

## Install

```r
# Bioconductor deps
BiocManager::install(c("Biostrings", "pwalign"))
# rtracklayer/GenomicRanges only if you want chain-file liftOver
remotes::install_github("your-org/xsmut")
```

## Quick start (R)

```r
library(xsmut)

muts <- read.csv(system.file("extdata/example_macaque_TP53.csv", package = "xsmut"))
#   sample_id gene aa_change       consequence
#   MM001     TP53 p.R175H   missense_variant ...

res <- xsmut_pipeline(
  gene        = "TP53",
  mutations   = muts,
  species     = "macaque",
  assembly    = "Mmul_10",
  cosmic      = "data/cosmic/Cosmic_MutantCensus_v101_GRCh38.tsv.gz",
  census_path = "data/cosmic/Cosmic_CancerGeneCensus_v101_GRCh38.tsv.gz"
)

res$gene_check         # in_census, tier, role, n_cosmic_mutations
res$model              # human canonical transcript, exons, CDS
res$cosmic_mutations   # all human COSMIC mutations for TP53
res$species_mutations  # your mutations + human_aa_pos / human_chrom / human_pos
res$plot               # patchwork; ggsave("TP53.pdf", res$plot, width = 12, height = 8)
```

Several genes:

```r
cos <- cosmic_load("data/cosmic/Cosmic_MutantCensus_v101_GRCh38.tsv.gz", genes = c("TP53","KRAS","PTEN"))
res <- lapply(c("TP53","KRAS","PTEN"), xsmut_pipeline, mutations = muts,
              species = "macaque", assembly = "Mmul_10", cosmic = cos)
```

Step-by-step (same thing the pipeline does):

```r
sp    <- ens_resolve_species("Macaca mulatta")
model <- get_gene_model("TP53")                                 # human
cos   <- cosmic_gene_mutations("TP53", cos)
mm    <- map_mutations_to_human(muts, "TP53", sp)
plot_gene_mutations(model, cos, mm, species_label = "Macaque (Mmul_10)")
```

## Shiny app

```r
xsmut::run_app()
```

Upload the COSMIC export, the (optional) census file and your mutation table; enter genes,
species and assembly; press **Run**. One tab per gene, with plot + both tables, PDF/CSV export.

## Input mutation table

| column        | required | notes                                                        |
|---------------|----------|--------------------------------------------------------------|
| `gene`        | yes      | human symbol (the ortholog is looked up for you)             |
| `aa_change`   | one of   | HGVS protein, e.g. `p.R175H`, `p.Q136*`, `p.P72fs`           |
| `aa_pos`      | these    | integer protein position in the species protein              |
| `consequence` | no       | drives colouring (missense / nonsense / frameshift / splice) |
| `sample_id`   | no       | counted per unique mutation                                  |
| `chrom`,`pos` | no       | species genomic coords; used only with a `chain` file        |

## Repository layout

```
R/ensembl.R      Ensembl REST client, species/assembly resolution, orthologs, protein→genome
R/cosmic.R       COSMIC file normalisation, gene-listed check, per-gene retrieval
R/gene_model.R   canonical transcript, exon/CDS/UTR model
R/crossmap.R     species aa → human aa → human genome (+ optional liftOver)
R/plot.R         gene track + upper/lower lollipop panels (ggplot2 + patchwork)
R/pipeline.R     xsmut_pipeline(), run_app()
inst/app/app.R   Shiny app
inst/extdata/    example macaque TP53 table
tests/           offline unit tests (network tests should be skipped on CRAN/CI without access)
```

## Caveats

* Ensembl REST serves only the current assembly per species (plus GRCh37 for human). For an
  older species assembly, annotate against it, then supply a chain file for genomic rows.
* Positions are drawn at the **first base of the human codon**; indels are drawn at their start.
* Ensembl rate-limits at ~15 req/s; requests are memoised for the session.
