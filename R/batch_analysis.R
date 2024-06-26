### Load libraries and functions ----
library(rnaseq)
library(tidyverse)
library(org.Hs.eg.db)
library(ensembldb)

source("anno/R/prepare_anno.R")
source('./R/kallisto_to_deseq2_util.R')
source('./R/get_files_from_bundle.R')
source('./R/batch_functions.R')
source('./R/intercepts_values.R')
source('./R/get_FCtable.R')
source('../rnaseq/R/batch_volcanos.R')
source('../rnaseq/R/volcano.R')
source('../rnaseq/R/batch_DE.R')
source('../rnaseq/R/gprofiler2_analysis.R')


### Set paths ----
# Kallisto
path_kallistos <- './deseq2_files/'

# Util files
path_util <- './files_util'
path_anno <- file.path(path_util, 'anno')
path_design <- file.path(path_util, 'simple_design.csv')
path_de_infos <- file.path(path_util, 'de_infos.csv')
path_pca_infos <- file.path(path_util, 'pca_infos.csv')
path_pca_meta <- file.path(path_util, 'pca_meta.csv')
path_volcano_infos <- file.path(path_util, 'volcano_infos.csv')

# DE output
path_out <- './dgea_output'
path_batch_de <- file.path(path_out,'de')
if (!dir.exists(path_batch_de)) dir.create(path = path_batch_de, recursive = TRUE)
path_txi_info <- file.path(path_out,'txi_infos')
if (!dir.exists(path_txi_info)) dir.create(path = path_txi_info, recursive = TRUE)
path_gene_list <- file.path(path_out,'gene_list')
if (!dir.exists(path_gene_list)) dir.create(path = path_gene_list, recursive = TRUE)
path_batch_deg_list <- file.path(path_out,'deg_list')
if (!dir.exists(path_batch_deg_list)) dir.create(path = path_batch_deg_list, recursive = TRUE)


# PCA output
path_pca <- file.path(path_out,'batch_pca')
if (!dir.exists(path_pca)) dir.create(path = path_pca, recursive = TRUE)
path_batch_pca_rds <- file.path(path_pca,'pca_rds')
if (!dir.exists(path_batch_pca_rds)) dir.create(path = path_batch_pca_rds, recursive = TRUE)
path_batch_pca_pdf <- file.path(path_pca,'pca_pdf')
if (!dir.exists(path_batch_pca_pdf)) dir.create(path = path_batch_pca_pdf, recursive = TRUE)

# Volcano output
path_batch_volcano <- file.path(path_out, 'volcano')
if (!dir.exists(path_batch_volcano)) dir.create(path = path_batch_volcano, recursive = TRUE)
path_batch_volcano_rds <- file.path(path_batch_volcano, 'volcano_rds')
if (!dir.exists(path_batch_volcano_rds)) dir.create(path = path_batch_volcano_rds, recursive = TRUE)
path_batch_volcano_pdf <- file.path(path_batch_volcano, 'volcano_pdf')
if (!dir.exists(path_batch_volcano_pdf)) dir.create(path = path_batch_volcano_pdf, recursive = TRUE)

### Make simplified annotations ----
prepare_anno(
    org = 'Homo sapiens',
    db = 'Ensembl',
    outdir = path_anno,
    force_download = TRUE
    )

anno <- list.files(
    path = path_anno,
    pattern='*.protein_coding.csv',
    full.names = TRUE
    )

### Load kallisto files with tximport ----
h5_files <- list.files(
    path = path_kallistos,
    pattern='*.h5',
    recursive = TRUE,
    full.names = TRUE
    )

txi <- rnaseq::import_kallisto(
    filenames = h5_files,
    anno = anno,
    ignoreTxVersion = TRUE
    )

samples <- list.dirs(path = path_kallistos, full.names = FALSE)[-1]
txi <- modify_colnames(df = txi, names = samples)
rownames(txi$length) <- gsub(
    pattern = '\\.[0-9]+',
    replacement = '',
    x = rownames(txi$length)
)

generate_info_txi(txi = txi, path = path_txi_info)

### Make design dataframe ----
group <- gsub(pattern = '_Rep[0-9]', replacement = '', x = samples)
simple_design <- data.frame(sample = samples, group = group)
readr::write_csv(x = simple_design, file = path_design)

### Make analysis ----
### Make batch analysis for differential expression
elem <- list(c('DU145', 'VCaP','LAPC4'), c('DMSO', 'R1881', '11KT', '11OHT'))
repo <- generate_comparison_file(elem = elem, path_write = path_de_infos)
batch_de(
    de_infos = path_de_infos,
    txi = txi,
    design = path_design,
    outdir = path_batch_de,
    force = TRUE,
    deg_output = path_gene_list
)

modify_batch_de_output(path = path_batch_de)
order_de_files(path = path_batch_de, directories = repo, file_type = '.csv', volcano = FALSE)

### Make batch application of principal component analysis (PCA)
generate_pca_files(
    samples = samples,
    out_meta = path_pca_meta,
    out_infos = path_pca_infos
)

rnaseq::batch_pca(
    pca_infos = path_pca_infos,
    txi = txi,
    metadata = path_pca_meta,
    outdir = path_batch_pca_pdf,
    r_objects = path_batch_pca_rds,
    force = TRUE
)

### Make batch application of volcanos
de_results <- generate_batch_volcano(
    path_de = path_batch_de,
    output = path_volcano_infos
)

batch_volcano(
    volcano_infos = path_volcano_infos,
    de_results = de_results,
    outdir = path_batch_volcano_pdf,
    r_objects = path_batch_volcano_rds,
    force = TRUE
)

order_de_files(
    path = path_batch_volcano_pdf,
    directories = repo,
    file_type = '.pdf',
    volcano = TRUE
)

fold.changes <- compute_foldchange(
    in_path = path_batch_de,
    fold.changes = c(NA, 1.1, 1.3, 1.5, 1.7, 2, 2.3, 2.5, 2.7, 3)
)

get_deg_list(
    outfile = path_batch_deg_list,
    infile = path_batch_de,
    fc.seuil = 1.3,
    p.adj = 0.05
)

sub_anno <- txi$anno %>% dplyr::select(c('id','symbol'))
make_intercepts_list(anno = sub_anno, files_dir = path_batch_deg_list)


