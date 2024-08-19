library(openxlsx)
library(tidyverse)
library(AnnotationDbi)

source('./scripts/R/util.R')

prepare_ranking <- function(rds.files, intercept.files,
                            cells = c('LAPC4', 'VCaP', 'DU145')){
    set.seed(seed = 42)

    genes.11oxy <- process.data.11oxy(
        cells = cells,
        path.genes = intercept.files,
        path.rds = rds.files
        ) %>%
        dplyr::bind_rows() %>%
        # Group rows based on ensembl_gene (eg. ids)
        dplyr::group_by(ensembl_gene) %>%
        # Simplify the data.frame by taking the maximum FC for each genes
        dplyr::summarize(
            log2FoldChange = log2FoldChange[which.max(abs(log2FoldChange))],
            .groups = 'drop'
        ) %>%
        # Ordered the dataframe by descending order of the log2FoldChange column
        dplyr::arrange(desc(log2FoldChange)) %>%
        # Edit column name
        dplyr::rename(ENSEMBL = ensembl_gene)

    # Map the ENSEMBL gene IDs to the first ENTREZ Ids
    entrez.ids <- switchIds(genes.11oxy$ENSEMBL, to = 'SYMBOL')

    # Make the ranked gene list
    ranks <- genes.11oxy %>%
        # Set the column 'ENTREZID' as the rownames
        tibble::column_to_rownames('ENSEMBL') %>%
        # Order the vector based on the values of log2FoldChange
        dplyr::arrange(desc(log2FoldChange)) %>%
        # Generate the named vector based on the anterior computetions
        { setNames(.$log2FoldChange, entrez.ids[rownames(.)]) }

    # Return ranked list
    ranks
}
