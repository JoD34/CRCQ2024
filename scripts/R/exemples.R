library(tidyverse)
library(UpSetR)
library(openxlsx)

create.rds_DESeq2_file <- function(){

    # Create path to store dummy RDS file
    path.output <- './scripts/exemples'
    if(!dir.exists(path.output)) dir.create(path.output, recursive = TRUE)

    # Forgo the creation of the file if it already exists
    filename <- file.path(path.output, 'dummy_DESeq2.rds')
    if(!file.exists(filename)) {

        # Generate a dummy deseq2 output file
        deseq2_dummy <- tibble(
            gene = c("ENSG00000139618", "ENSG00000141510", "ENSG00000157764", "ENSG00000155657", "ENSG00000146648"),
            baseMean = runif(5, 50, 1000),          # Random base mean values
            log2FoldChange = rnorm(5, 0, 2),        # Random log2 fold changes
            lfcSE = rnorm(5, 0.1, 0.5),             # Random standard error of log2 fold changes
            stat = rnorm(5, 0, 10),                 # Random test statistics
            pvalue = runif(5, 0, 1),                # Random p-values
            padj = p.adjust(runif(5, 0, 1), method = "BH") # Adjusted p-values using Benjamini-Hochberg
        )

        # Save RDS object to a specific location
        saveRDS(object = deseq2_dummy, file = filename)

    }

    # Return the relative path to the RDS file
    filename
}

create.intersect_file <- function(){

    # Create path to store dummy RDS file
    path.output <- './scripts/exemples'
    if(!dir.exists(path.output)) dir.create(path.output, recursive = TRUE)

    # Forgo the creation of the file if it already exists
    filename <- file.path(path.output, 'dummy_intersect.xlsx')

    if(!file.exists(filename)){

        genes_11KT <- c("ENSG00000139618", "ENSG00000141510", "ENSG00000157764")
        genes_11OHT <- c("ENSG00000139618", "ENSG00000155657", "ENSG00000146648")
        genes_R1881 <- c("ENSG00000139618", "ENSG00000157764", "ENSG00000146648")

        genes.list <- list('11KT' = genes_11KT, '11OHT' = genes_11OHT, 'R1881' = genes_R1881)
        write.gList2xlsx(gene.list = list(LAPC4 = genes.list), outdir = filename)
    }

    # Return filename
    filename
}

### Function: exemple.col.filter ----
###
### exemple.col.filter()
###
##      Get the column's name for a DESeq2 outputed RDS file
##
##  Valeur
##
##      < character >  Column of the DESeq2 object
##
exemple.col.filter <- function(){
    path <- './dgea_output/rds'
    files <- list.files(path = path, pattern = '.*_de.rds', full.names = TRUE)

    # Read RDS object
    rds.obj <- readRDS(files[1])

    # Return the column's name
    colnames(rds.obj)
}
