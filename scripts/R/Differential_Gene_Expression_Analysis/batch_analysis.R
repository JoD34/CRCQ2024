### Load libraries and functions ----
library(rnaseq)
library(tidyverse)

### Load scripts ----
scripts <- c(
    # My code
    './scripts/R/Differential_Gene_Expression_Analysis/kallisto_to_deseq2_util.R',
    './scripts/R/Differential_Gene_Expression_Analysis/get_files_from_bundle.R',
    './scripts/R/Differential_Gene_Expression_Analysis/batch_functions.R',
    './scripts/R/Venn_Diagram/intercepts_values.R',
    './scripts/R/get_FCtable.R',
    './scripts/R/Venn_Diagram/venndiagram.R',
    './scripts/R/Venn_Diagram/unicGenes_from_Intersect.R',

    # rnaseq package from Charles Joly-Beauparlant - Differential Expression Analysis
    '../rnaseq/R/batch_volcanos.R',
    '../rnaseq/R/volcano.R',
    '../rnaseq/R/batch_DE.R',
    '../rnaseq/R/gprofiler2_analysis.R',

    # anno package from Arnaud Droit labo
    "anno/R/prepare_anno.R"
)

# Loading
lapply(scripts, source)


### Set paths ----
# Base paths
path_kallistos <- './deseq2_files/'
path_util <- './files_util'
path_out <- './output_deseq2_pval0.01'

# Utility files
paths_util <- list(
    anno = file.path(path_util, 'anno'),
    design = file.path(path_util, 'simple_design.csv'),
    de_infos = file.path(path_util, 'de_infos.csv'),
    pca_infos = file.path(path_util, 'pca_infos.csv'),
    pca_meta = file.path(path_util, 'pca_meta.csv'),
    volcano_infos = file.path(path_util, 'volcano_infos.csv')
)

# Create directories if they don't exist
paths <- list(
    de = file.path(path_out, 'de'),
    rds = file.path(path_out, 'rds'),
    txi_infos = file.path(path_out, 'txi_infos'),
    gene_list = file.path(path_out, 'gene_list'),
    deg_list = file.path(path_out, 'deg_list'),
    pca = file.path(path_out, 'batch_pca'),
    volcano = file.path(path_out, 'volcano'),
    venn_diagrams = file.path(path_out, 'venn_diagrams')
)

# Create directories
lapply(paths, function(p) {
    if (!dir.exists(p)) dir.create(p, recursive = TRUE)
})

# Specific subdirectories
sub_dirs <- list(
    pca_rds = file.path(paths$pca, 'pca_rds'),
    pca_pdf = file.path(paths$pca, 'pca_pdf'),
    volcano_rds = file.path(paths$volcano, 'volcano_rds'),
    volcano_pdf = file.path(paths$volcano, 'volcano_pdf'),
    venn_ARpos = file.path(paths$venn_diagrams, 'Hormones_ARp'),
    gene_list_de = file.path(paths$gene_list, 'batch_de_output')
)

# Create specific subdirectories
lapply(sub_dirs, function(p) {
    if (!dir.exists(p)) dir.create(p, recursive = TRUE)
})


### Make simplified annotations ----
prepare_anno(
    org = 'Homo sapiens',
    db = 'Ensembl',
    outdir = paths_util$anno,
    force_download = TRUE
    )

anno <- list.files(
    path = paths_util$anno,
    pattern='*.protein_coding.csv',
    full.names = TRUE
    )

### Load kallisto files with tximport ----
## List all files ending in .h5 - exclusive to kallisto files in this directory
h5_files <- list.files(
    path = path_kallistos,
    pattern='*.h5',
    recursive = TRUE,
    full.names = TRUE
    )

# Import kallisto data files
txi <- rnaseq::import_kallisto(
    filenames = h5_files, # Kallisto data outpu
    anno = anno, # Annotation file
    ignoreTxVersion = TRUE # Ignore the ENSEMBL version to simplify annotation
    )

# Customize kallisto dataframe - add sample IDs to columns
samples <- list.dirs(path = path_kallistos, full.names = FALSE)[-1]
txi <- modify_colnames(df = txi, names = samples)
rownames(txi$length) <- gsub('\\.[0-9]+', '', x = rownames(txi$length))

# Extract and customize the intermediary data generated by the tximport function
#   Summarize the trireplicates by giving only the mean and std
generate_info_txi(txi = txi, path = paths$txi_infos, cells = c('DU145','LAPC4','VCaP'))

### Make design dataframe ----
group <- gsub(pattern = '_Rep[0-9]', replacement = '', x = samples)
simple_design <- data.frame(sample = samples, group = group)
readr::write_csv(x = simple_design, file = paths_util$design)

### Make analysis ----
### Make batch analysis for differential expression

## Generating file describing the various comparison to do
elem <- list(c('DU145', 'VCaP','LAPC4'), c('DMSO', 'R1881', '11KT', '11OHT'))
repo <- generate_comparison_file(elem = elem, path_write = paths_util$de_infos)

# Execute the Differential Genes Expression Analysis using DESeq2
batch_de(
    de_infos = paths_util$de_infos,   # Relative path describing the analysis to do
    txi = txi,                        # Kallisto data importations
    design = paths_util$design,       # Relative path to the design experiment
    outdir = paths$de,                # Relative path to output (.csv) directory
    force = FALSE,                    # Should the files be overwritten?
    r_objects = paths$rds,            # Relative path to output (.csv) directory
    deg_output = sub_dirs$gene_list_de # Where to store the gene lists of DEG
)

# Modify DESeq2 output files - add information on counts (raw & TPM)
modify_batch_de_output(path = paths$de)

# Order DESeq2 output fils into subdirectories
order_de_files(path = paths$de, directories = repo, file_type = '.csv', volcano = FALSE)

### Make batch application of principal component analysis (PCA)

## Generate a file stating various parameters needed for a PCA analysis
generate_pca_files(
    samples = samples,
    out_meta = paths_util$pca_meta,     # Relative path to store the PCA metadata
    out_infos = paths_util$pca_infos    # Relative path to store the PCA information
)

## Make PCA analysis as stated in the metadata generated in the previous step
rnaseq::batch_pca(
    pca_infos = paths_util$pca_infos,   # Relative path to the description of the PCA analysis
    txi = txi,                          # Kallisto information
    metadata = paths_util$pca_meta,     # Relative path to metadata of the PCA analysis
    outdir = sub_dirs$pca_pdf,          # Relative path to store the image of the results
    r_objects = sub_dirs$pca_rds,       # Relative path to store the RDS object
    force = FALSE                       # Should the files be overwritten?
)

### Make batch application of volcanos
de_results <- generate_batch_volcano(
    path_de = paths$de,
    output = paths_util$volcano_infos
)

# Generate Volcano plots base on the information generated in the earlier step
batch_volcano(
    volcano_infos = paths_util$volcano_infos,
    de_results = de_results,
    outdir = sub_dirs$volcano_pdf,
    r_objects = sub_dirs$pca_rds,
    force = FALSE
)

order_de_files(
    path = sub_dirs$volcano_pdf,
    directories = repo,
    file_type = '.pdf',
    volcano = TRUE
)

# Generate a table quantifying the number fo significant genes per fold change threshold
fold.changes <- compute_foldchange(
    in_path = paths_util$de_infos,
    fold.changes = c(NA, 1.1, 1.3, 1.5, 1.7, 2, 2.3, 2.5, 2.7, 3)
)

# Indicate pairwise subsetting of DEG data you want to do
fc.seuil.values <- c(1.3, NA, NA)
pAdj.values <- c(0.05, 0.05, NA)

# Iterate over both previously constructed vectors
purrr::walk2(fc.seuil.values, pAdj.values, ~ {
    # Generate files representing the ENSEMBL Gene IDs of the gene respecting
    #   the Fold Change threshold as well as the pValue adjusted
    get_deg_list(
        outfile = paths$deg_list,
        infile = paths$de,
        fc.seuil = .x,
        p.adj = .y
    )
})

# Generate Venn Diagrams for various genes lists from 'get_deg_list' output
constraintes <- list.dirs(path = paths$deg_list, recursive = FALSE, full.names = FALSE)
map(constraintes, ~ {

    # Get files name
    output <- file.path(paths$venn_diagrams, .x)
    input <-  file.path(paths$deg_list, .x)

    # Generate Venn diagrams and corresponding list
    triple_set_venn(
        output.dir = output,
        input.dir = input,
        cell = c('DU145', 'VCaP', 'LAPC4'),
        hormones = c('R1881', '11KT', '11OHT')
    )

    # Get gene list for genes exclusive to 11oxy hormones
    intersect2unic(
        input.dir = file.path(output, 'H_vs_DMSO', 'H_vs_DMSO_GeneLists.xlsx'),
        output.dir = file.path(paths$gene_list, .x),
        sheets = c('LAPC4', 'VCaP'),
        cols = c('11KT', '11OHT'),
        intersect = TRUE
        )

    output <- file.path(sub_dirs$venn_ARpos, .x)
    triple_set_venn(
        output.dir = output,
        input.dir = input,
        cell = c('LAPC4', 'VCaP'),
        hormones = c('R1881', '11KT', '11OHT'),
        byCell = FALSE
    )
    })

sub_anno <- txi$anno %>% dplyr::select(c('id','symbol'))
make_intercepts_list(anno = sub_anno, files_dir = paths$deg_list)


