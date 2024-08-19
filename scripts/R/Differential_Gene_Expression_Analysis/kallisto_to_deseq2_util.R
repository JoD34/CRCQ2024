# Load requiered libraries ----
library(biomaRt)
library(tidyverse)
# make_annotation_file ----
###
### make_annotation_file(path)
###
##  Summary
##      Generate .csv file of transcripts' annotation
##
##  Args
##      path: (character vector) path to kallisto quantification file
##
##
##  Examples
##      make_annotation_file("data/")
##
make_annotation_file <- function(name, paths, path_out=getwd()){

    paths <- sub(pattern='.h5', replacement='.tsv', x=paths)

    # Get all unique transcipt's id
    trans <- lapply(paths, function(f){
        read.csv(f, sep='\t') %>%
            dplyr::select('target_id') %>%
            as.matrix %>%
            as.character
    }) %>%
        unlist() %>%
        unique()

    mart <- useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl")

    # Retrieve annotation of transcripts
    attr <- c(
        'ensembl_transcript_id_version',
        'ensembl_gene_id',
        'external_gene_name',
        'transcript_biotype'
    )

    anno <- getBM(
        attributes=attr,
        values=trans,
        filters = 'ensembl_transcript_id_version',
        mart=mart
    )

    names(anno) <- c('id', 'ensembl_gene', 'symbol', 'transcript_type')

    anno$entrez_id <- NA
    anno <- anno[ , c('id', 'ensembl_gene', 'symbol', 'entrez_id', 'transcript_type')]

    # Write file
    out_dir <- file.path(path_out, paste("annotation",name,".csv", sep='_'))
    write.csv(x=anno, file=out_dir, row.names = FALSE)
    out_dir
}


# modify_colnames ----
###
### modify_colnames(df, names)
###
##  Modify the dataframe abundance, counts and length of the tximport function
##
##  Arguments
##
##      df (list): output of the function tximport
##      names (charactor vector): new column names
##
##  Return
##
##      df (list): list with modified column names for the dataframe mentionned
##                 above
##
##  Exemples
##
##  modify_colnames(txi, new_column_names)
##
modify_colnames <- function(df, names){
    colnames(df$abundance) <- names
    colnames(df$counts) <- names
    colnames(df$length) <- names
    colnames(df$fpkm) <- names

    df
}


# make_design ----
###
### make_design(samples)
###
##      Generate the experimental design table
##
##  Arguments
##
##      Sample (character vector): Names given to the samples
##      cnames (character vector): Column's name
##
##  Return
##
##      df (dataframe) design table
##
##  Exemples
##
##  make_design(c(sample1, sample2))
##
make_design <- function(samples, cnames, condition){

    # Generate design table from scratch
    df <- stringr::str_extract_all(
        string = samples,
        pattern="[a-zA-Z0-9]+",
        simplify=TRUE
        ) %>%
        as.data.frame()

    # Format design dataframe
    if (any(grepl(pattern='ENZA', x=df))){
            colnames(df) <- cnames
        df$replicate <- ifelse(df$replicate == '', df$enza, df$replicate) %>%
            gsub(pattern='Rep', replacement='', x=.)
        df$enza <- as.numeric(df$enza == 'ENZA')
    } else {
        colnames(df) <- cnames[-which(cnames == 'enza')]
    }
    df$sample <- samples

    # Reformat columns and add factors
    df <- df %>%
        dplyr::relocate(sample, .before=everything()) %>%
        dplyr::relocate(replicate, .after=everything())

    df <- df %>% mutate(across(everything(), ~ as.factor(.x)))

    # Reassign the level for the column of importance in the DESeq2 analysis
    if(condition == "enza"){
        df$enza <- relevel(df$enza, ref='0')
    } else {
        pref <- c('DMSO', 'R1881', '11KT', '11OHT')
        ref <- sapply(X=pref, FUN=function(h) sum(df$hormone == h))
        df$hormone <- relevel(x=df$hormone, ref=pref[ref>0][1])
    }

    df
}
