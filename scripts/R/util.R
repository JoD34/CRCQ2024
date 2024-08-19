library(openxlsx)

write2xlsx <- function(df, sheetNames, out.dir, file.name){
    if (any(class(df) == 'data.frame')) write2xlsx.df(df, sheetNames, out.dir, file.name)
}

### Function: write2xlsx ----
###
### write2xlsx(df, sheetNames, filename)
###
##  Get log2 Fold Change for differentially expressed genes (DEG) for a given
##      hormones (R1881, 11KT, 11OHT) and a given cell lineage (VCaP, LAPC4)
##      for the comparaison Hormone vs DMSO
##
##  Arguments
##
##  df < data.frame > : Dataframe to write on .xslx sheet
##  cell < character > : cell lineage's name
##  filename < character > : directory path to files containing DEG
##
##  Valeur
##
##  NA - write a dataframe to a .xlsx file
##
##  Exemples
##
##  <expression(s) démontrant l'utilisation de la fonction>
##
write2xlsx.df <- function(df, sheetNames, out.dir, file.name){

    # Assigned values
    start.col = 2L
    start.row = 2L

    # Create a new workbook (eg new excel file)
    wb <- openxlsx::createWorkbook()

    # Create a new Sheet to work on
    sheet <- openxlsx::addWorksheet(wb, sheetName = sheetNames)

    # Write the dataframe to the appropriate sheet
    openxlsx::writeDataTable(wb, sheetNames, x = df, startCol = start.col, startRow = start.row)

    # Save workbook to the specified path
    openxlsx::saveWorkbook(wb, file = file.path(out.dir, file.name), overwrite = TRUE)

}

### Function: get.genes_infos ----
###
### get_fold_change(df, cell, filename)
###
##  Get log2 Fold Change for differentially expressed genes (DEG) for a given
##      hormones (R1881, 11KT, 11OHT) and a given cell lineage (VCaP, LAPC4)
##      for the comparaison Hormone vs DMSO
##
##  Arguments
##
##  df < list > : DEG for a given intersect of Venn diagram
##  cell < caracter > : cell lineage's name
##  filename < caracter > : directory path to files containing DEG
##
##  Valeur
##
##  < data.frame > log2 Fold associated to each DEG
##
##  Exemples
##
##
get.genes_infos <- function(df, file.infos, filename,
                            col.filter = c('ensembl_gene', 'log2FoldChange')){

    # Get gene list in character vector
    genes <- unlist(df, use.names = FALSE) %>% purrr::discard(is.na)

    # Separate the hormone's name if the set is an intersect of two hormones
    h <- if(grepl(pattern = '\\.', x = file.infos$hormone)){

        strsplit(file.infos$hormone, "\\.")[[1]]

    } else {

        file.infos$hormone
    }

    # Formatting the path to the directory
    file <- paste(rep(file.infos$cell, length(h)), h,'DMSO_de.rds', sep = '__')
    filename %<>%  file.path(file)

    # Retrieve the Log2FoldChange for the genes of interest
    dea <- filename %>%
        purrr::map(
            # Assigne the function readRDS to each element of filename
            ~ readRDS(.) %>%
                # Get specific rows (info on genes of interest)
                dplyr::filter(ensembl_gene %in% genes) %>%
                # Get the two specified columns
                dplyr::select(dplyr::all_of(col.filter))
        )

    # Apply modification if a gene has more than one fold change assigned
    if( length(dea) != 1 ){

        dea <- dplyr::bind_rows(dea)

    } else {

        dea <- dea[[1]]

    }

    # Return the information in a dataframe type value
    as.data.frame(dea)
}

### Function: switchIds ----
###
### switchIds(df, from, to)
###
##      Map IDs between two differents types. In case of a mapping 1:many, the
##          function will return the first mapping.
##
##  Arguments
##
##      ids < character > : IDs to get the mapping from
##      from <character> : KeyType of the IDs
##      to < character > : KeyType to map IDs to
##
##  Valeur
##
##      < data.frame > Mapping of IDs between type
##
##  Exemples
##
##      ensembl_ids <- c("ENSG00000139618", "ENSG00000284733", "ENSG00000141510",
##                       "ENSG00000157764", "ENSG00000155657")
##      exemple.map <- switchIds(ids = ensembl_ids)
##      exemple.map2 <- switchIds(ids = ensembl_ids, from = 'ENSEMBL', to = 'ENTREZID')
##
##      Get list of ID type:
##          library(org.Hs.eg.db)
##          library(AnnotationDbi)
##          AnnotationDbi::keytypes(org.Hs.eg.db)
##
switchIds <- function(ids, from = 'ENSEMBL', to = 'SYMBOL'){

    # Map the ENSEMBL gene IDs to the first ENTREZ Ids
    AnnotationDbi::mapIds(
        x = org.Hs.eg.db, # Database from which to get the IDs
        key = ids, # Input
        keytype = from, # Identifier type of the input
        column = to, # Identifier type of the output(s)
        multiVals = "first" # How to deal with mapping 1:many -> take the first
    )

}


### Function: process.data.11oxy ----
###
### process.data.11oxy(cells, path.genes, path.rds)
###
##      Retrieve gene list from .xlsx file and match it to a DESeq2 output file.
##      Usage: retrieve expression information of genes exclusive to a condition.
##
##  Arguments
##
##      cells < character > : Cell lineage for which to ger data (Used as sheet in .xlsx file)
##      path.genes <character> : Relative Paths to were gene of intersept are stored
##      path.rds < character > : Relative Paths to the DESeq2 output objects
##      col.filter < character > : Columns of the DESeq2 output to retrieve
##
##  Valeur
##
##      < data.frame > ENSEMBL IDs with the selected columns for genes in path.genes
##
##  Exemples
##
##      source('./scripts/R/exemples.R')
##      cells <- c("LAPC4")
##      path.rds <- create.rds_DESeq2_file()
##      path.genes <- create.intersect_file()
##      process.data.11oxy(cells, path.genes, path.rds)
##
##      Get all available column: exemple.col.filter()
##
process.data.11oxy <- function(cells, path.genes, path.rds,
                               col.filter = c('ensembl_gene', 'log2FoldChange')){

    # Iterate over all cell lineage()
    lapply(X = cells, function(sheet){

        # Load genes ENSEMBL IDs that are 11oxy specific
        genes <- openxlsx::read.xlsx(xlsxFile = path.genes, sheet = sheet) %>%
            dplyr::slice(-1) %>%
            dplyr::select(-matches('R1881|^X[0-9]+$'))

        # Get information from rds files
        lapply(X = names(genes), FUN = function(h){
            get.genes_infos(
                df = genes[[h]],
                file.infos = list(cell = sheet, hormone = h),
                filename = path.rds,
                col.filter = col.filter
            )
        }) %>% dplyr::bind_rows()
    })
}
