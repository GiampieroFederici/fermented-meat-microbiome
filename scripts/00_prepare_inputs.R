# Materialise disposable inputs from the journal's Tables S1-S10.
# Usage from the repository root:
# Rscript --vanilla scripts/00_prepare_inputs.R /path/to/supplementary_tables .
# No supplementary workbook is modified or copied into this repository.
args <- commandArgs(TRUE)
if (length(args) != 2L) stop('Usage: Rscript scripts/00_prepare_inputs.R TABLE_DIRECTORY REPOSITORY_ROOT')
if (!requireNamespace('readxl', quietly = TRUE)) stop('Install the R package readxl before preparing inputs.')
tables <- normalizePath(args[1], winslash = '/', mustWork = TRUE)
root <- normalizePath(args[2], winslash = '/', mustWork = TRUE)
options(stringsAsFactors = FALSE)
for (p in c('02_CONTROLLI', 'checks', 'outputs', 'work/data/04',
            'work/derived', '03_RETTIFICHE_PRONTE/03_MANOSCRITTO')) {
  dir.create(file.path(root, p), recursive = TRUE, showWarnings = FALSE)
}
stems <- c('Table_S1_Literature_inventory', 'Table_S2_Literature_matrix',
           'Table_S3_Samples_accessions_QC', 'Table_S4_Taxonomic_profiles',
           'Table_S5_Profiler_concordance', 'Table_S6_Alpha_diversity',
           'Table_S7_Beta_diversity', 'Table_S8_Taxon_groups_and_screening',
           'Table_S9_Community_sensitivities', 'Table_S10_StrainPhlAn_Lsakei')
files <- file.path(tables, paste0(stems, '.xlsx'))
if (any(!file.exists(files))) stop('Missing supplementary workbook(s): ', paste(basename(files[!file.exists(files)]), collapse = ', '))
write_tsv <- function(d, relative) {
  path <- file.path(root, relative)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  # The plotting scripts use the title-case boolean spelling in the journal exports.
  for (j in seq_along(d)) if (is.logical(d[[j]])) d[[j]] <- ifelse(is.na(d[[j]]), NA_character_, ifelse(d[[j]], 'True', 'False'))
  write.table(d, path, sep = '\t', row.names = FALSE, quote = TRUE, na = 'NA')
}
read_sheet <- function(n, sheet) {
  as.data.frame(readxl::read_xlsx(files[n], sheet = sheet, .name_repair = 'minimal', guess_max = 100000), check.names = FALSE)
}
# These are runtime text views of the accompanying workbooks, excluded by .gitignore.
for (n in seq_along(files)) {
  for (s in readxl::excel_sheets(files[n])) {
    if (grepl('[/\\\\]', s)) stop('Unexpected worksheet name: ', s)
    write_tsv(read_sheet(n, s), file.path('02_CONTROLLI', stems[n], paste0(s, '.tsv')))
  }
}
put <- function(d, n, s) write_tsv(d, file.path('02_CONTROLLI', stems[n], paste0(s, '.tsv')))
alpha <- read_sheet(6, 'Depth_samples')[, 1:10]
stopifnot(nrow(alpha) == 100L, identical(names(alpha)[c(1, 8:10)], c('sample_name', 'shannon', 'PC1', 'PC2')))
put(alpha, 6, 'Alpha_Kraken')
community <- read_sheet(7, 'PERMANOVA_PERMDISP')
figure5 <- community[community$tool == 'kraken', ]
figure5 <- figure5[match(c('product', 'country', 'raw_material'), figure5$factor), ]
figure5$PERMANOVA_q <- p.adjust(figure5$PERMANOVA_p, 'BH')
put(figure5, 7, 'Figure5_statistics')

# Translate only the historical display labels expected by the retained plotting code.
# Sample IDs, values, grouping memberships and statistical methods stay unchanged.
historical <- c('I_Ciauscolo salami_P' = 'Ciauscolo',
                'I_Felino salami_P' = 'Felino salami',
                'I_Piedmont salami_P' = 'Piedmont salami',
                'C_Miao Sour Meat_P' = 'Chinese sausage (miao)',
                'C_Mutton sausage_M' = 'Mutton sausage',
                'S_Cecina meat_B' = 'cecina', 'S_Cured tongue_B' = 'cured tongue',
                'S_Dry aged beef_B' = 'dry-aged beef',
                'S_Chorizo salami_P' = 'DFS chorizo',
                'S_Salchichon salami_P' = 'DFS salchichon')
rename_product <- function(x) {
  out <- unname(historical[as.character(x)])
  if (anyNA(out)) stop('Unrecognised product label in supplementary tables.')
  out
}
ann <- read_sheet(3, 'Samples_100')[, c('sample_name', 'srr', 'batch', 'country', 'product', 'raw_material', 'day')]
ann$product <- rename_product(ann$product)
write_tsv(ann, 'work/data/sample_annotation.tsv')
alpha$product <- rename_product(alpha$product)
write_tsv(alpha, 'work/derived/samples_kraken.tsv')
write_tsv(community, 'work/derived/community_statistics.tsv')
write_tsv(read_sheet(6, 'Omnibus_KW'), 'work/derived/paper_alpha_omnibus.tsv')
contrasts <- read_sheet(6, 'Conover_Shannon')
for (g in c('group1', 'group2')) {
  i <- contrasts$factor == 'product'
  contrasts[[g]][i] <- rename_product(contrasts[[g]][i])
}
write_tsv(contrasts, 'work/derived/paper_conover_BH.tsv')
# Check the current ten-product table rather than asserting an unavailable historical audit.
bh_error <- max(vapply(split(contrasts, contrasts$factor), function(d) max(abs(p.adjust(d$p, 'BH') - d$p.adj)), numeric(1)))
stopifnot(is.finite(bh_error), bh_error < 1e-10)
write_tsv(data.frame(check = 'Current ten-product Conover/BH matrix', pairs = nrow(contrasts),
                     max_absolute_error = bh_error, result = 'PASS'), 'work/derived/paper_conover_reproduction.tsv')

k <- read_sheet(4, 'Kraken_100')
K <- as.matrix(k[-1]); storage.mode(K) <- 'double'; rownames(K) <- k[[1]]
stopifnot(nrow(K) == 100L, ncol(K) == 1380L, all(is.finite(K)), all(K >= 0))
pseudocount <- min(K[K > 0]) / 2
logK <- log(K + pseudocount)
CLR <- logK - rowMeans(logK)
write_tsv(data.frame(sample_name = rownames(CLR), CLR, check.names = FALSE), 'work/data/04/clr_matrix_kraken_species.tsv')
# These two batches have no values removed by the published filtering rule.
# Their percent abundance matrices are therefore recoverable from Table S4.
# The five other batches retain pre-filter information and are supplied separately.
live_ann <- read_sheet(3, 'Samples_100')
for (batch in c('ciau', 'mutton')) {
  ids <- live_ann$sample_name[live_ann$batch == batch]
  B <- 100 * K[ids, colSums(K[ids, , drop = FALSE]) > 0, drop = FALSE]
  write_tsv(data.frame(species = colnames(B), t(B), check.names = FALSE),
            file.path('work/data/kraken_batches', batch, 'species_relabund_matrix.tsv'))
}
mapping <- read_sheet(8, 'Group_mapping')
mapping <- mapping[, c('historical_species_label', 'guild', 'mean_abund_allsamples_pct')]
names(mapping)[1] <- 'species'
write.csv(mapping, file.path(root, 'work/data/species_to_group_original.csv'), row.names = FALSE, na = 'NA')

refs <- read_sheet(10, 'References_182')
labels <- read_sheet(10, 'FigureS1_tip_labels')
idx <- match(refs$Accession, labels$original_tip)
stopifnot(nrow(refs) == 182L, !anyNA(idx), !anyDuplicated(labels$original_tip))
parts <- strsplit(labels$printed_label[idx], ' | ', fixed = TRUE)
stopifnot(all(lengths(parts) == 5L))
parts <- do.call(rbind, parts)
stopifnot(identical(unname(parts[, 5]), refs$Accession))
source_codes <- c(meat = 'M', vegetal = 'V', animal_host = 'A', dairy = 'D',
                  starter = 'T', cereal = 'C', environment = 'E', seafood = 'F')
refs$Source_code <- unname(source_codes[labels$source_category[idx]])
stopifnot(!anyNA(refs$Source_code))
refs$Geography_label <- parts[, 2]
refs$Year_label <- parts[, 3]
refs$Strain_label <- parts[, 4]
display_category <- parts[, 1]
display_category[display_category == 'animal_host'] <- 'animal host'
refs$Figure_label_header <- paste(display_category, parts[, 2], parts[, 3], sep = ' | ')
refs$Figure_label_name <- ifelse(parts[, 4] == 'NR', refs$Accession, parts[, 4])
refs$Full_human_label <- labels$printed_label[idx]
write_tsv(refs, 'sources/reference_labels_verified.tsv')
palette <- labels[!duplicated(labels$source_category), c('source_category', 'marker_colour', 'label_colour')]
stopifnot(setequal(palette$source_category, names(source_codes)))
palette <- palette[match(names(source_codes), palette$source_category), ]
palette$code <- unname(source_codes)
palette$source_category[palette$source_category == 'animal_host'] <- 'animal host'
names(palette)[2:3] <- c('marker_triangle_colour', 'label_terminal_branch_colour')
palette <- palette[, c('code', 'source_category', 'marker_triangle_colour', 'label_terminal_branch_colour')]
write_tsv(palette, 'sources/Figure6_source_palette.tsv')
cat('PREPARE_INPUTS_PASS: Tables S1-S10 read; runtime inputs materialised; no workbook modified.\n')
