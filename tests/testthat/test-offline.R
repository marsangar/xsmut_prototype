# Tests that need no network access

test_that("parse_aa_pos handles common HGVS forms", {
  expect_equal(parse_aa_pos(c("p.R175H", "p.Q331*", "p.E285_K286del", "p.?", "", NA, "p.Met1?")),
               c(175L, 331L, 285L, NA, NA, NA, 1L))
})

test_that("alignment map handles gaps on both sides", {
  # human:  A-CDE F
  # target: ABCD-EF
  m <- xsmut:::.aln_map("A-CDEF", "ABCD-EF")
  # target residues: A B C D E F  -> human 1 NA 2 3 4 5
  expect_equal(m, c(1L, NA, 2L, 3L, 4L, 5L))
})

test_that("intron compression keeps exons contiguous and monotone", {
  ex <- data.table::data.table(exon_rank = 1:3, start = c(100, 1000, 5000), end = c(200, 1100, 5100), feature = "CDS")
  tr <- xsmut:::.make_transform(ex, compress_introns = TRUE)
  x  <- tr$f(c(100, 150, 200, 600, 1000, 5100))
  expect_true(all(diff(x) > 0))
  expect_equal(x[2] - x[1], 50)
})

test_that("legacy and new COSMIC formats normalise to the same schema", {
  new <- tempfile(fileext = ".tsv")
  data.table::fwrite(data.table::data.table(
    GENE_SYMBOL = "TP53", TRANSCRIPT_ACCESSION = "ENST00000269305", MUTATION_CDS = "c.524G>A",
    MUTATION_AA = "p.R175H", MUTATION_DESCRIPTION = "missense_variant", GENOMIC_MUTATION_ID = "COSV52661038",
    COSMIC_SAMPLE_ID = "S1", PRIMARY_SITE = "lung", PRIMARY_HISTOLOGY = "carcinoma",
    MUTATION_SOMATIC_STATUS = "Confirmed somatic variant", CHROMOSOME = "17",
    GENOME_START = 7675088L, GENOME_STOP = 7675088L, STRAND = "-", GENOMIC_WT_ALLELE = "C", GENOMIC_MUT_ALLELE = "T"), new, sep = "\t")
  old <- tempfile(fileext = ".tsv")
  data.table::fwrite(data.table::data.table(
    `Gene name` = "TP53", `Accession Number` = "ENST00000269305", `Mutation CDS` = "c.524G>A",
    `Mutation AA` = "p.R175H", `Mutation Description` = "Substitution - Missense",
    GENOMIC_MUTATION_ID = "COSV52661038", ID_sample = "S1", `Primary site` = "lung",
    `Primary histology` = "carcinoma", `Mutation somatic status` = "Confirmed somatic variant",
    `Mutation genome position` = "17:7675088-7675088", `Mutation strand` = "-"), old, sep = "\t")
  a <- cosmic_load(new); b <- cosmic_load(old)
  expect_setequal(names(a), names(b))
  expect_equal(a$start, b$start)
  expect_equal(a$aa_pos, 175L)
  expect_equal(cosmic_gene_listed("TP53", mutations = a)$n_cosmic_mutations, 1L)
})
