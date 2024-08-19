### Load libraries ----
library(tidyverse)
library(clusterProfiler)
library(RColorBrewer)
library(org.Hs.eg.db)
library(enrichplot)
library(ggupset)
library(openxlsx)

source('./scripts/R/util.R')

# Functions ----

###
### prep.genes.enricher(cells, paths)
###
##  <description de ce que fait la fonction en une phrase>
##
##  Arguments
##
##  <argument>: <type (vecteur, liste, etc.), signification
##               et valeurs admissibles, le cas échéant>
##
##  Valeur
##
##  <type (vecteur, liste, etc.) et valeur retournée par la
##   fonction>
##
##  Exemples
##
##  <expression(s) démontrant l'utilisation de la fonction>
##

prep.genes.enricher <- function(cells, paths){
    df <- process.data.11oxy(
        cells = cells,
        path.genes = paths$in.genes,
        path.rds = paths$in.rds,
        col.filter = c('ensembl_gene', 'pvalue', 'padj', 'log2FoldChange')
        ) %>%
        dplyr::bind_rows() %>%
        # Group rows based on ensembl_gene (eg. ids)
        dplyr::group_by(ensembl_gene) %>%
        # Simplify the data.frame by taking the maximum FC for each genes
        dplyr::slice_max(order_by = abs(log2FoldChange), n = 1) %>%
        # Ordered the dataframe by descending order of the log2FoldChange column
        dplyr::arrange(desc(log2FoldChange)) %>%
        # Edit column name
        dplyr::rename(ENSEMBL = ensembl_gene, log2FC = log2FoldChange) %>%
        # Add a column to represented the direction of change
        dplyr::mutate(diffExpr = dplyr::case_when(
            log2FC > 0 & pvalue < 0.05 ~  'UP',
            log2FC < 0 & pvalue < 0.05 ~  'DOWN',
            pvalue > 0.05 ~  'NO'
        )) %>%
        # Remove genes annotated as 'NO'
        #Shouldn't be useful; Genes are filtered in the making of the dataset
        dplyr::filter(!diffExpr == 'NO')

    # Change IDs from ENSEMBL to SYMBOLS
    symbols <- switchIds(ids = df$ENSEMBL, to = 'SYMBOL')
    df$ENSEMBL <- symbols

    # Reformat Maped data
    symbols <- data.frame(
        ENSEMBL = names(symbols),
        SYMBOLS = symbols,
        row.names = NULL
    )

    # Matching the appropriate type of IDs
    df %<>% dplyr::rename(SYMBOL = ENSEMBL)

    # Split by Fold change changes (upregulated or downregulated)
    df.split <- split(df, df$diffExpr)
    df.split$BOTH <- df

        # Return a list with all information
    list(symbols = symbols, gene.express = df.split)
}

###
### <signature de la fonction>
###
##  <description de ce que fait la fonction en une phrase>
##
##  Arguments
##
##  <argument>: <type (vecteur, liste, etc.), signification
##               et valeurs admissibles, le cas échéant>
##
##  Valeur
##
##  <type (vecteur, liste, etc.) et valeur retournée par la
##   fonction>
##
##  Exemples
##
##  <expression(s) démontrant l'utilisation de la fonction>
##

prep.gs.files <- function(paths){

    # Get file's path containing gene sets
    gs.files <- list.files(paths$gs, pattern = '.*\\.gmt$', full.names = TRUE)

    # Get background genes
    bg.genes <- readr::read_csv(paths$bg, col_names = "background") %>%
        dplyr::pull() %>%
        switchIds(from = 'ENSEMBL', to = 'SYMBOL')

    # Susbet gene sets to the genes present in my background
    gs <- lapply(X = gs.files, FUN = function(f.path){

        # Read .gmt file containing gene sets and subset to fit background genes
        f <- clusterProfiler::read.gmt(gmtfile = f.path) %>%
            dplyr::filter(gene %in% bg.genes)

        # Construct the output path
        f.filename <- f.path %>%
            stringr::str_replace(pattern = 'raw', replacement = 'modified') %>%
            stringr::str_replace(pattern = '\\.gmt$', replacement = '\\.rds')

        # Save RDS object
        if(file.exists(f.filename)) file.remove(f.filename)
        saveRDS(object = f, file = f.filename)

        # Prompt for verbose
        file <- f.filename %>% stringr::str_extract('(c[0-9]+|h).*')
        print(paste0('Modifcation completed : ', file))

        # Return the modified gene set
        f

    })

    # Return the modified gene sets
    setNames(gs, stringr::str_extract(gs.files, '(c[0-9]+|h).*'))
}

###
### <signature de la fonction>
###
##  <description de ce que fait la fonction en une phrase>
##
##  Arguments
##
##  <argument>: <type (vecteur, liste, etc.), signification
##               et valeurs admissibles, le cas échéant>
##
##  Valeur
##
##  <type (vecteur, liste, etc.) et valeur retournée par la
##   fonction>
##
##  Exemples
##
##  <expression(s) démontrant l'utilisation de la fonction>
##

visualize_data.pathway_enrichment_analysis <- function(df, outdir, setsId, genes, bg){
    df <- openxlsx::read.xlsx(file.path(paths$out,'BOTHregulated_MSigDb_enrichR.xlsx'), sheet = 'c3') %>%
        tibble::column_to_rownames(var = 'ID')
    df$Description <- gsub('_', ' ', x = df$Description) %>% tolower() %>% toTitleCase()


    # Create an EnrichResult
    enrichres <- new(
        Class = "enrichResult",
        redeable = FALSE,
        result = df,
        pvalueCutoff = 0.05,
        pAdjustedMethod = 'BH',
        qvalueCutoff = 0.2,
        organism = 'human',
        ontology = 'UNKNOWN',
        gene = genes,
        keytype = 'UNKNOWN',
        universe = bg,
        gene2symbol = character(0),
        geneSets = NULL)
    )
}

## Main code ----
main <- function(){
### Set path to variables
paths <- list(
    bg = './dgea_output/background_genes.csv',
    in.rds = './dgea_output/rds',
    in.genes = './dgea_output/venn_diagrams/FC_NA_pAdj_0.05/H_vs_DMSO/H_vs_DMSO_GeneLists.xlsx',
    out = './dgea_output/pea',
    gs = './dgea_output/gene_sets/raw',
    visualize.PEA = '/dgea_output/pea/visualisation'
)

# Create the output directory if it isn't already created
if(!dir.exists(paths$out)) dir.create(path = paths$out, showWarnings = FALSE, recursive = TRUE)
if(!dir.exists(paths$visualize.PEA)) dir.create(path = paths$visualize.PEA, showWarnings = FALSE, recursive = TRUE)

# Assuring all path exists
if(!all(file.exists(unlist(paths)))) stop("Not all file's path exist")

# Set other variables
cells <- c('LAPC4','VCaP')
count.cutoff <- 5
pval.cutoff <- 0.05

# Prepare datasets
gene.infos <- prep.genes.enricher(cells = cells, paths = paths)

# Modify information about the gene sets
gene.sets <- prep.gs.files(paths = paths)
gs <- split(gene.sets, stringr::str_extract(names(gene.sets), '(c[0-9]|h)'))

# Run Enricher ----
purrr::imap(gene.infos$gene.express, ~ {
    wb <- openxlsx::createWorkbook()
    genes <- .x$SYMBOL
    direction <- .y
    purrr::iwalk(gs, ~ {
        res <- purrr::map(.x, ~ clusterProfiler::enricher(gene = genes, TERM2GENE = .x)@result) %>%
            purrr::list_rbind() %>%
            dplyr::filter(pvalue <= pval.cutoff & Count > count.cutoff)

        # Add a new WorkBook sheet & Write data
        openxlsx::addWorksheet(wb, sheetName = .y)
        openxlsx::writeDataTable(wb, sheet = .y, x = res, xy = c(2, 2))

        # Adjust columns width for the data table
        openxlsx::setColWidths(wb, sheet = .y, cols = 2:(ncol(res) + 1), widths = "auto")

        # Generate directory to store visualisation tools
        output.visualisation <- file.path(paths$visualize.PEA, direction)
        if(!dir.exists(output.visualisation)) dir.create(output.visualisation)

        # Visualise data
        visualize_data.pathway_enrichment_analysis(
            df = res,
            outdir = output.visualisation,
            setId = .y
        )

    })

    # Save enriched data in a .xlsx file
    filename <- paste0(.y, paste('regulated', 'MSigDb','enrichR.xlsx', sep = '_' ))
    openxlsx::saveWorkbook(wb, file = file.path(paths$out, filename), overwrite = TRUE)

}, .progress = TRUE)
}
