library(tidyverse)
library(fs)
library(openxlsx)

source('./scripts/R/util.R')
source('./scripts/R/Venn_Diagram/heatmap_ARpos_11oxyVsEnza_functions.R')

# Set variables of interest ----
# Hard code variables
cell_lineage <- c('LAPC4','VCaP')
hormones <- c('11KT', '11OHT')

# Associated regex
regex.cell_lineage <-paste(cell_lineage, collapse = '|')
regex.hormones <-paste(hormones, collapse = '|')
pattern.files <- paste0("^(",regex.cell_lineage,")__(",regex.hormones,").*__DMSO_de\\.rds$")

# Update cell version
version <- 'v6'

# Set directories' path to variables ----
dir.output <- file.path(getwd(), 'dgea_output/heatmap_ARpos_11oxy')
file.output <- file.path(dir.output, paste(version, 'heatmap_11oxy_exclusif', sep = '_'))
dir.de <- file.path(getwd(), 'dgea_output/de')
dir.venns <- file.path(getwd(), 'dgea_output/venn_diagrams')
dir.commun_genes <- file.path(dir.venns, 'all_venn', c('H_vs_DMSO', 'HE_vs_DMSO'))
dir.de.rds <- file.path(path = 'dgea_output/rds')

# Create directories that don't already exists
if(!dir.exists(dir.de)) dir.create(dir.de)
if(!dir.exists(dir.output)) dir.create(dir.output)
if(!dir.exists(dir.de.rds)) dir.create(dir.de.rds)

# Get comparison of interest for both cell lineage, shorten file names ----
r_object <- list.files(path = dir.de.rds, pattern = pattern.files) %>%
    # Split by cell lineage
    split(x = ., f = sub( paste0('^(', regex.cell_lineage, ').*'), '\\1', .))

# Get Excel file absolute path to directories ----
file.commun_genes <- fs::dir_ls(path = dir.commun_genes, regexp = '*.xlsx')

# Generate HeatMap ----
purrr::walk(r_object, ~ {

    # Extract information from filename
    infos <- sub("^(.*)__(.*)__.*de.rds$", "\\1,\\2", .x) %>%
        str_split(',', simplify = TRUE) %>%
        apply(2, FUN = unique) %>%
        setNames(c('name.lineage','sets'))

    # Get list of genes exclusive to 11oxy
    genes.11oxy <- file.commun_genes %>%
        purrr::map(~ get_genes_excel_table(
            xlsx.file = .x,
            sheet = infos$name.lineage,
            cols = hormones)) %>%
        unlist %>%
        unique

    # Generate the dataset used by the heatmap
    data.heatmap <- fs::path(dir.de.rds, .x) %>%
        # Read DESeq2 output data and select the ENSEMBL IDs and the log2FC
        purrr::map(~ {
            base::readRDS(.x) %>%
                dplyr::filter(ensembl_gene %in% genes.11oxy) %>%
                dplyr::select(ensembl_gene, log2FoldChange)
        }) %>%
        # Join all the extracted information by their ENSEMBL IDs
        purrr::reduce(~ full_join(.x, .y, by = 'ensembl_gene')) %>%
        # Remove rows having NAs
        na.omit() %>%
        # Switch de type from data.frame to tibble
        as_tibble() %>%
        # Replace the ENSEMBL IDs for the corresponding SYMBOLS
        dplyr::mutate(ensembl_gene = switchIds(ids = ensembl_gene)[ensembl_gene]) %>%
        # Set the newly formed SYMBOLS' column as rownames
        tibble::column_to_rownames('ensembl_gene') %>%
        # Change the column names
        setNames(infos$sets) %>%
        # Shuffle columns' order to cluster those with 'ENZA' on one side
        dplyr::select(!contains('ENZA'), contains('ENZA'))

    # Build output directory
    filedir <- paste(file.output, infos$name.lineage, sep = '_')

    # Set variables values
    shorten <- c(FALSE, FALSE, TRUE)
    number <- c(TRUE, FALSE, FALSE)
    title <-  paste0('Fold Change of 11 oxy in ', infos$name.lineage, ' - Effect of enzalutamide')

    walk2(shorten, number, ~{
        # Generate the various filenames
        filename <- paste(
            filedir,
            ifelse(.x, 'shorten', 'long'),
            ifelse(.y, 'with_number', ''),
            'unscaled.pdf',
            sep = '_'
        )

        # Generate general heatmap
        Make_HeatMap(
            plot.data = data.heatmap,
            filename = filename,
            title = title,
            scale = FALSE,
            shorten = .x,
            display_numbers = .y,
            ccluster = FALSE
        )
    })

    if(FALSE){
    # Addition representation for top 100 genes up- or down- regulated
    if(nrow(plot.data) >= 200){
        lapply(c(TRUE, FALSE), function(.dec){ # Look at upregulated and downregulated solely
            lapply(c('11KT', '11OHT'), function(.pattern){ # Look at both 11oxy separatly
                subset_heatmap_100(
                    plot.data = data.heatmap,
                    pattern = .pattern,
                    file.dir = filedir,
                    decreasing = .dec
                )
            })
        })
        }
    }
    })

shorten <- c(FALSE, FALSE, TRUE)
number <- c(TRUE, FALSE, FALSE)
title <- paste('Fold Change of 11oxy genes - Combined AR+ cell lineage')

walk2(shorten, number, ~ {
    file.out <- paste0(
        version,
        '_ARpos',
        ifelse(.x,'_shorten', ''),
        ifelse(.y, '_withNumbers', ''),
        '_unscaled.pdf'
    )

    # Generate HeatMap for gene expression on jointe on VCaP & LAPC4 for 11oxy
    generate_whole_heatmap(
        pattern.files = pattern.files,
        rds_files = dir.de.rds,
        dir.commun = fs::path(dir.venns, 'Hormones_ARp'),
        xlsx.sheetname = hormones,
        xlsx.col = 'Both.ENSEMBL',
        shorten = .x,
        display_numbers = .y,
        scale = FALSE,
        filename.heatmap = fs::path(dir.output, file.out),
        title.heatmap = title
    )
})

    # Generate HeatMap for intersection VCaP-LAPC4 for AR+ question
generate_whole_heatmap(
    pattern.file = pattern.files,
    rds_files = dir.de.rds,
    dir.commun = dir.commun_genes,
    xlsx.sheetname = cell_lineage,
    xlsx.col = c('11KT','11OHT', '11KT.11OHT'),
    filename.heatmap = file.path(dir.output, paste(version, '11oxy_vs_11oxyENZA','scaled.pdf',sep = '_')),
    title.heatmap = 'Fold Change - 11oxy specific'
)

# Generate heatmap for intersection 11oxy only for VCaP and LAPC4
generate_whole_heatmap(
    pattern.file = pattern.files,
    rds_files = dir.de.rds,
    dir.commun = dir.commun_genes,
    xlsx.sheetname = cell_lineage,
    xlsx.col = c('11KT.11OHT'),
    filename.heatmap = file.path(dir.output, paste(version, 'jonction11KT11OHT','scaled.pdf',sep = '_')),
    title.heatmap = 'Fold Change - 11oxy specific'
)
