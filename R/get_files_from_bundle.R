library(tidyverse)
library(xlsx)

###
### get_files(path, regex)
###
##  Summary
##      Get the path to all the needed files for further loading
##
##  Args
##      path (character vector): path to repository
##      regex (character vector): regex for the files I want to retrieve
##
##
##  Examples
##      get_files(path='\data', regex='VCaP_DMSO')
##      get_files(path='\data', regex='VCaP_DMSO|R1881')
##
get_files <- function(path, regex){
    # Separate regex
    regex <- regex %>%
        stringr::str_split_1(pattern='_')

    # Get all files
    work.files <- file.path(path, list.files(path))

    # Take what the regex asked for
    for(elem in regex){
        work.files <- work.files[grepl(pattern=elem, x=work.files, perl=TRUE)]
    }

    # Different return options
    if(length(regex) == 2) return(work.files)
    if (any(grepl(pattern='ENZA', x=regex))){
        return(work.files[grepl(pattern='ENZA', x= work.files)])
    }
    work.files[!grepl(pattern='ENZA', x= work.files)]
}

###
### get_samples_name(path, regex)
###
##  Summary
##      Get the path to all the needed files for further loading
##
##  Args
##      path (character vector): path to repository
##      regex (character vector): regex for the files I want to retrieve
##
##
##  Examples
##      get_files(path='\data', regex='VCaP_DMSO')
##      get_files(path='\data', regex='VCaP_DMSO|R1881')
##
get_samples_name <- function(bundle){

    # Extract the sample's name from the file's directory
    bundle <- bundle %>%
        stringr::str_extract_all(pattern='[a-zA-Z0-9_]+', simplify=TRUE) %>%
        '['(, 2)

    bundle
}

###
### make_comparison_file(path)
###
##  Summary
##      Get all combination of samples component.
##
##  Args
##      comp (list): components to generate combination from;
##          It is assumed that the list is in order of encompassment, making it
##          easier to group element.
##      path (character vector): directory to write combinations file
##
##
##  Examples
##      make_comparison_file()
##      make_comparison_file(path='./data')
##
make_comparison_file <- function(comp, path='.'){

    # First element is a segragation argument
    base <- unlist(comp[1])

    # Get all combination into a single 'word'
    comb <- expand.grid(comp) %>%
        apply(MARGIN=1, FUN=function(x) paste(x, collapse='_'))

    # Select pairwise comparison
    gp <- comb %>%
        split(f = base) %>%
        lapply(FUN = function(df) combn(df, 2, simplify = FALSE)) %>%
        unlist() %>%
        matrix(ncol = 2, byrow = TRUE) %>%
        data.frame()

    gp_enza <- gp
    gp_enza$X1 <- paste(gp_enza$X1, '_ENZA', sep = "")
    gp_enza$X2 <- paste(gp_enza$X2, '_ENZA', sep = "")


    # Generate Enza comparison
    enza <- comb %>%
        lapply(FUN=function(x) c(x, paste(x, 'ENZA', sep='_'))) %>%
        split(f=base) %>%
        unlist() %>%
        matrix(byrow=TRUE, ncol=2) %>%
        data.frame()

    w <- rbind(gp, gp_enza, enza)
    colnames(w) <- c('base', 'condition')

    # Write data into a .csv file
    readr::write_csv(
        x=data.frame(w),
        file=file.path(path, 'combinaison.csv')
        )
}

#create_directories ----
###
### create_directories(path)
###
##  Summary
##      Make directories for all output files
##
##  Examples
##      create_directories()
##
create_directories <- function(util){

    dir.create(path=util)
    dir.create(path='dgea_output')
    dir.create(path='dgea_output/volcanos')
    dir.create(path='dgea_output/data_visualization')
    dir.create(path='dgea_output/deg_list')
    dir.create(path='dgea_output/dds_output')
    dir.create(path='dgea_output/excel')
    dir.create(path='dgea_output/count_genes')
    dir.create(path='dgea_output/count_tmp')
}

###
### make_comparison_file(path)
###
##  Summary
##      Get all combination of samples component.
##
##  Args
##      comp (list): components to generate combination from;
##          It is assumed that the list is in order of encompassment, making it
##          easier to group element.
##      path (character vector): directory to write combinations file
##
##
##  Examples
##      make_comparison_file()
##      make_comparison_file(path='./data')
##
write_txi_file <- function(txi, comparison, path, namefile){

    # Get informations from txi object
    df_count <- get_raw_count_anno_df(txi) %>%
        dplyr::select(-c('ensembl_gene', 'entrez_id')) %>%
        get_stats(df=., celem=comparison)

    df_tmp <- get_tpm_anno_df(txi) %>%
        dplyr::select(-c('ensembl_gene', 'entrez_id')) %>%
        get_stats(df=., celem=comparison)


    # Change colnames of info data
    anno_col <- c('id', 'symbol', 'transcript_type')
    df_count <- remake_colNames(knames=anno_col, df=df_count, add_on='count')
    df_tmp <- remake_colNames(knames=anno_col, df=df_tmp, add_on='tmp')

    # Get count informations
    write_file <- df_count %>%
        inner_join(df_tmp, by = join_by(id, symbol, transcript_type) )

    readr::write_csv(
        x=write_file,
        file=file.path(path, paste0(namefile, '.csv'))
    )
}

###
### get_stats(path)
###
##  Summary
##      Get all combination of samples component.
##
##  Args
##      comp (list): components to generate combination from;
##          It is assumed that the list is in order of encompassment, making it
##          easier to group element.
##      path (character vector): directory to write combinations file
##
##
##  Examples
##      make_comparison_file()
##      make_comparison_file(path='./data')
##
get_stats <- function(df, celem){
    res <- df[,1:3]
    for( c in celem ){

        # Get column's informations
        v_loc <- grep(pattern=c, x=colnames(df))
        v <- df[,v_loc]

        # Generate a stat's subdataframe
        stat_df <- data.frame(mean=rowMeans(v), std=apply(v, 1, sd))
        colnames(stat_df) <- c(
            paste('mean', c, sep='_'),
            paste('std_dev', c, sep='_')
        )

        res <- cbind(res, stat_df)
    }
    res
}


###
### make_comparison_file(path)
###
##  Summary
##      Get all combination of samples component.
##
##  Args
##      comp (list): components to generate combination from;
##          It is assumed that the list is in order of encompassment, making it
##          easier to group element.
##      path (character vector): directory to write combinations file
##
##
##  Examples
##      make_comparison_file()
##      make_comparison_file(path='./data')
##
write_deg <- function(de, anno, path, namefile){

    # Add raw count information from tximport
    de$id <- rownames(de)

    add_on <- anno %>%
        dplyr::select(, -c('entrez_id', "ensembl_gene"))

    de <- as.data.frame(de) %>%
        inner_join(add_on, by = join_by(id))

    de <- de %>%
        dplyr::relocate(colnames(add_on), .before=everything())

    readr::write_csv(
        x=de,
        file=file.path(path, paste0(namefile, '.csv'))
    )
}

###
### make_comparison_file(path)
###
##  Summary
##      Get all combination of samples component.
##
##  Args
##      comp (list): components to generate combination from;
##          It is assumed that the list is in order of encompassment, making it
##          easier to group element.
##      path (character vector): directory to write combinations file
##
##
##  Examples
##      make_comparison_file()
##      make_comparison_file(path='./data')
##
remake_colNames <- function(knames, df, add_on){
    cnames <- colnames(df)
    change <- cnames[!(cnames %in% knames)]

    new_col <- sapply(X=change, FUN=function(x) paste(add_on, x, sep='_'))
    colnames(df) <- c(knames, new_col)
    df
}

###
### make_comparison_file(path)
###
##  Summary
##      Get all combination of samples component.
##
##  Args
##      comp (list): components to generate combination from;
##          It is assumed that the list is in order of encompassment, making it
##          easier to group element.
##      path (character vector): directory to write combinations file
##
##
##  Examples
##      make_comparison_file()
##      make_comparison_file(path='./data')
##
write_excel <- function(path_files){

}
