library(rlang)
library(furrr)
library(fs)
library(future)
library(magrittr)
library(rnaseq)
library(RColorBrewer)
library(tidyverse)

### generate_comparison_file(elem, path_write)
###
##  Summary
##      Make comparisons files used as instruction for rnaseq::batch_de
##
##  Args
##      elem (list): component for comparisons; first element of the list is
##                   used to separate by class
##      path_write (character vector): directory indicating where to write file
##
##  Examples
##      elem <- list(c('DU145', 'VCaP','LAPC4'),
##                   c('DMSO', 'R1881', '11KT', '11OHT'))
##      generate_comparison_file(elem = elem,path_write = './util')
##
generate_comparison_file <- function(elem, path_write){

    # Combined all elements by expanding the character vectors
    joinded <- apply(X = expand.grid(elem), 1, paste, collapse = '_')

    # Separated by specific classes
    sep <- split(joinded, elem[[1]])

    # Generate first dataframe of comparison, template for another one
    df <-  lapply(X = sep, combn, m = 2) %>%
        do.call(what = rbind) %>%
        matrix(ncol = 2, byrow = TRUE) %>%
        as.data.frame

    # Adding extra comparisons that doesn't match with previous comparisons
    extra <- lapply(X = sep, FUN = function(l) expand.grid(l[1], l[-1])) %>%
        do.call(what = rbind) %>%
        'colnames<-'(c('V1','V2'))
    extra$V2 <- paste(extra$V2, 'ENZA', sep='_')

    # Make final comparison dataframe, containing other comparisons
    res <- rbind(
        df,
        sapply(X=df, paste, 'ENZA', sep = '_'),
        unlist(sep) %>% data.frame(V1 = paste(., 'ENZA', sep='_'), V2 = .),
        extra
        ) %>%
        'row.names<-'(1:nrow(.)) %>%
        'colnames<-'(c('contrast_2', 'contrast_1'))

    id_mod <- stringr::str_split(string=res$contrast_2, pattern = '_', simplify=TRUE)[,-1] %>%
        apply(MARGIN = 1, FUN = paste, collapse = '_') %>%
        sub(pattern = '_$', replacement = '')

    # Construct dataframe used in rnaseq::batch_de
    c1 <- stringr::str_split(res$contrast_1, '_', simplify=TRUE)
    c <- apply(c1[, 2:3], 1, paste, collapse = '_') %>%
        sub(pattern = '_$', replacement = '') %>%
        paste(c1[, 1], . , sep = '__')

    id_de <-  paste(c, id_mod, sep='__')
    de_infos <- data.frame(id_de =id_de) %>%
        dplyr::mutate(
            group = 'group', formula = '~group',
            filter = 4, count_matrix = NA) %>%
        cbind(., res)

    # Write the dataframe to a file without rownames
    readr::write_csv(x = de_infos, file = path_write)

    # Return the comparison id for directory storage
    id_de
}

### generate_pca_files(samples, path_metadata, path_pca_infos)
###
##  Summary
##      Make comparisons files used as instruction for rnaseq::batch_de
##
##  Args
##      elem (list): component for comparisons; first element of the list is
##                   used to separate by class
##      path_write (character vector): directory indicating where to write file
##
##  Examples
##      elem <- list(c('DU145', 'VCaP','LAPC4'),
##                   c('DMSO', 'R1881', '11KT', '11OHT'))
##      generate_comparison_file(elem = elem,path_write = './util')
##
generate_pca_files <- function(samples, paths){
    generate_pca_metadata(id = samples, output = paths$pca_meta)
    generate_pca_infos(output = paths$pca_infos)
}

### generate_pca_metadata(id, path_metadata)
###
##  Summary
##      Make metadata files for pca ultimately used in  rnaseq::batch_pca
##
##  Args
##      id (character vector): samples names
##      path_metadata (character vector): directory for where to write file
##
##  Examples
##      id <- c('DU145_DMSO_Rep1','DU145_DMSO_Rep1','DU145_DMSO_Rep1')
##      generate_pca_metadata(id, path_metadata='./util)
##
generate_pca_metadata <- function(id, output){

    group <- gsub(pattern = '_Rep[A-z0-9]', replacement = '', x = id)
    shape <- grepl(pattern = 'ENZA', x = group) + 21

    # Associate colors corresponding to the sample's hormone
    color_sets <- stringr::str_split(string = group, pattern = '_', simplify = TRUE)[, 2]
    unic_elem <- unique(color_sets)
    dict <- setNames(
        object = brewer.pal(n = length(unic_elem), name = 'Set1'),
        nm = unic_elem
    )
    color <- dict[color_sets]

    # Reformat groups
    new_group <- as.character(
        grepl('ENZA', group) +
            grepl('DU145', group) * 1 +
            grepl('LAPC4', group) * 3 +
            grepl('VCaP', group) * 5
    )
    df <- data.frame(id = id, shape = shape, color = color, group = paste0('group', new_group))

    # Write metadata used for pca graph's generation
    readr::write_csv(x = df, file = output)
}


### generate_pca_metadata(id, path_metadata)
###
##  Summary
##      Make metadata files for pca ultimately used in  rnaseq::batch_pca
##
##  Args
##      id (character vector): samples names
##      path_metadata (character vector): directory for where to write file
##
##  Examples
##      id <- c('DU145_DMSO_Rep1','DU145_DMSO_Rep1','DU145_DMSO_Rep1')
##      generate_pca_metadata(id, path_metadata='./util)
##
generate_pca_infos <- function(output){
    title = c(
        'PCA analysis:  lignee DU145, sans enzalutamide',
        'PCA analysis:  lignee DU145, avec enzalutamide',
        'PCA analysis:  lignee LAPC4, sans enzalutamide',
        'PCA analysis:  lignee LAPC4, avec enzalutamide',
        'PCA analysis:  lignee VCaP, sans enzalutamide',
        'PCA analysis:  lignee VCaP, avec enzalutamide'
    )
    infos <- data.frame(id_plot=c('DU145', 'DU145_ENZA','LAPC4', 'LAPC4_ENZA',
                                  'VCaP', 'VCaP_ENZA')) %>%
        dplyr::mutate(
            group='group',
                group_val = paste0('group', 1:6),
            use_normalisation = 'none',
            use_ruv = FALSE,
            use_combat = FALSE,
            min_counts = 5,
            id_metadata = 'id',
            size = 2,
            title = title,
            legend.position = 'bottom',
            legend.box = 'vertical',
            show_names = TRUE
    )

    # Write pca infos files
    readr::write_csv(x = infos, file = output)
}


### generate_pca_metadata(id, path_metadata)
###
##  Summary
##      Make metadata files for pca ultimately used in  rnaseq::batch_pca
##
##  Args
##      id (character vector): samples names
##      path_metadata (character vector): directory for where to write file
##
##  Examples
##      id <- c('DU145_DMSO_Rep1','DU145_DMSO_Rep1','DU145_DMSO_Rep1')
##      generate_pca_metadata(id, path_metadata='./util)
##
order_de_files <- function(path, directories, file_type, volcano){

    idx <- ifelse(volcano, 5, 4)
    d <- stringr::str_split(string = directories, pattern = '__', simplify = TRUE)[,-1] %>%
        apply(MARGIN = 1, FUN = paste, collapse='__') %>%
        unique

    lapply(X = d, FUN = function(c) {
        p <- file.path(path, c)
        if (!dir.exists(p)){dir.create(path = p)}

        mv_files <- list.files(
            path = path,
            pattern = paste0(c, file_type),
            full.names = TRUE,
            recursive = TRUE
            )

        names <- stringr::str_split(
            string = mv_files,
            pattern = '/',
            simplify = TRUE)[, idx]

        file.rename(from = mv_files, to = file.path(p, names))
        })
}


### generate_pca_metadata(id, path_metadata)
###
##  Summary
##      Make metadata files for pca ultimately used in  rnaseq::batch_pca
##
##  Args
##      id (character vector): samples names
##      path_metadata (character vector): directory for where to write file
##
##  Examples
##      id <- c('DU145_DMSO_Rep1','DU145_DMSO_Rep1','DU145_DMSO_Rep1')
##      generate_pca_metadata(id, path_metadata='./util)
##
generate_info_txi <- function(txi, path, cells){
    subfunc <- function(FUN){

        # Extract information fromt txi based on a specific function
        elem <- FUN(txi = txi) %>% 'rownames<-'(.$ensembl_gene)

        # Remove columns related to annotation and change column names
        subelem <- elem %>% dplyr::select(-colnames(txi$anno))
        new.colnames <- colnames(subelem) %>%
            stringr::str_remove(pattern = '_Rep[0-9]') %>%
            unique()

        # Split data by unique column names and compute the mean and standard deviation
        stat <- map(new.colnames, ~ {

            # Split the previous data.frame by unique columns names
            set <- subelem %>%
                dplyr::select(dplyr::matches(.x)) %>%
                as.matrix()

            # Generate new data.frame
            dplyr::tibble(
                # Generate the column names and assign the corresponding values
                !!paste(.x, '_mean') := rowMeans(set),
                !!paste(.x, '_std') := apply(set, 1, sd)
            )

        }) %>%
            purrr::list_cbind()

        # Return the newly generated data.frame
        stat
    }

    # Function to rename the column names appropriately
    format_col <- function(df, suffix){
        df %>%
            # On its own, rename_with replace the column names
            dplyr::rename_with( ~ {
                # Make a custom function to generate the new columns name
                stringr::str_replace(., '.mean', paste0('_mean_', suffix))
                }) %>%
            dplyr::rename_with( ~ {
                # Make a custom function to generate the new columns name
                stringr::str_replace(., 'std', paste0('_std_', suffix))
            }) #%>%
            # Integrate the row names within the data.frame
           # dplyr::mutate(id = rownames(df))
    }

    # Generate 'count' and 'tpm' values
    count <- subfunc(FUN = get_raw_count_anno_df) %>% format_col(suffix = 'counts')
    tpm <- subfunc(FUN = get_tpm_anno_df) %>% format_col(suffix = 'tpm')

    # Full join between the two generated data.frame
    df <- dplyr::full_join(x = count, y = tpm, by = 'id')

    walk(cells, ~ {
        sub.df <- df %>%
            dplyr::select(matches(.x)) %>%
            dplyr::rename_with( ~ { stringr::str_replace(., paste0('^', .x, '_'), '') }) %>%
            dplyr::mutate(id = df$id) %>%
            dplyr::left_join(txi$anno, by = 'id') %>%
            dplyr::select(ensembl_gene, symbol, entrez_id, everything(), -id, -transcript_type)

        # Write data.frame to a .csv file
        readr::write_csv(x = sub.df, file = fs::path(path, paste0(.x,'.csv')))
    })
}


### generate_pca_metadata(id, path_metadata)
###
##  Summary
##      Make metadata files for pca ultimately used in  rnaseq::batch_pca
##
##  Args
##      id (character vector): samples names
##      path_metadata (character vector): directory for where to write file
##
##  Examples
##      id <- c('DU145_DMSO_Rep1','DU145_DMSO_Rep1','DU145_DMSO_Rep1')
##      generate_pca_metadata(id, path_metadata='./util)
##
order_volcano_files <- function(path, directories){
    d <- stringr::str_split(string = directories, pattern = '__', simplify = TRUE)[,-1] %>%
        apply(MARGIN = 1, FUN = paste, collapse='__') %>%
        unique
}


### generate_pca_metadata(id, path_metadata)
###
##  Summary
##      Make metadata files for pca ultimately used in  rnaseq::batch_pca
##
##  Args
##      id (character vector): samples names
##      path_metadata (character vector): directory for where to write file
##
##  Examples
##      id <- c('DU145_DMSO_Rep1','DU145_DMSO_Rep1','DU145_DMSO_Rep1')
##      generate_pca_metadata(id, path_metadata='./util)
##
generate_batch_volcano <- function(path_de, output){
    # List files
    de_files <- list.files(path_de, recursive = TRUE, full.names = TRUE)

    # Get differential expression files' id
    de <- fs::path_file(de_files) %>% stringr::str_remove(., pattern = '.csv')

    de_results <- map(de_files, read.csv) %>% setNames(de)

    # Generate titles for each volcano plot to be generated
    elem <- stringr::str_split(string = de, pattern ='__', simplify = TRUE)
    title <- paste('Differential Expression:', elem[, 2], 'vs',
                   elem[, 3], 'on', elem[, 1], 'cells') %>%
        gsub(pattern = '_', replacement = ' with ', x = .)

    # Generate the dataframe for the volcano information
    volcano_infos <- data.frame(id_plot = de, id_de = de) %>%
        dplyr::mutate(
            y_axis = 'padj',
            p_threshold = 0.05,
            fc_threshold = 1.3,
            title = title,
            show_signif_counts = TRUE,
            show_signif_lines = 'vertical',
            show_signif_color = TRUE,
            col_up = '#b00c0c',
            col_down = '#09078a',
            size = 2
            )

    # Write the volcano plot information in the util directory
    readr::write_csv(x = volcano_infos, file = output)

    de_results
}


### generate_pca_metadata(id, path_metadata)
###
##  Summary
##      Make metadata files for pca ultimately used in  rnaseq::batch_pca
##
##  Args
##      id (character vector): samples names
##      path_metadata (character vector): directory for where to write file
##
##  Examples
##      id <- c('DU145_DMSO_Rep1','DU145_DMSO_Rep1','DU145_DMSO_Rep1')
##      generate_pca_metadata(id, path_metadata='./util)
##
modify_batch_de_output <- function(path) {

    f <- list.files(path, full.names = TRUE)

    walk(f, ~ {
        df <- read.csv(.x) %>% dplyr::select(-c('transcript_type', 'id'))
        fold_change <- 2^df$log2FoldChange
        df <- df %>% dplyr::mutate(fold_change, .after = log2FoldChange)
        df <- add_tmp_counts(df = df, txi = txi, path = .x)
        readr::write_csv(x = df, file = .x)
    })
}


### generate_pca_metadata(id, path_metadata)
###
##  Summary
##      Make metadata files for pca ultimately used in  rnaseq::batch_pca
##
##  Args
##      id (character vector): samples names
##      path_metadata (character vector): directory for where to write file
##
##  Examples
##      id <- c('DU145_DMSO_Rep1','DU145_DMSO_Rep1','DU145_DMSO_Rep1')
##      generate_pca_metadata(id, path_metadata='./util)
##
add_tmp_counts <- function(df, txi, path = p){
    rep <- paste0('Rep', 1:3)
    comparison <- stringr::str_split(string=path, pattern='/', simplify=TRUE)[,4] %>%
        sub(pattern='\\.csv', replacement='',x=.) %>%
        stringr::str_split(string=., pattern ='__', simplify=TRUE)

    cols <- list(
        paste(comparison[1], comparison[3], rep, sep='__'),
        paste(comparison[1], comparison[2], rep, sep='__')
    )

    # Add counts and tpm information
    df <- get_tpm_stats(txi = txi, df = df , cols = cols) %>%
        get_counts_stats(counts = txi$counts, df = ., cols = cols)
    df[order(df$padj),]
}


### generate_pca_metadata(id, path_metadata)
###
##  Summary
##      Make metadata files for pca ultimately used in  rnaseq::batch_pca
##
##  Args
##      id (character vector): samples names
##      path_metadata (character vector): directory for where to write file
##
##  Examples
##      id <- c('DU145_DMSO_Rep1','DU145_DMSO_Rep1','DU145_DMSO_Rep1')
##      generate_pca_metadata(id, path_metadata='./util)
##
get_counts_stats <- function(counts, df, cols){
    stats <- lapply(X = cols, FUN = function(x){
        title <- stringr::str_split(string = x[1], pattern = '__', simplify = TRUE) %>%
            '['(, -c(1,ncol(.)))

        x <- gsub(pattern = '__', replacement = '_', x = x)

        sub_counts <- as.data.frame(counts) %>% dplyr::select(all_of(x))
        stats <- data.frame(rowMeans(sub_counts), apply(sub_counts, 1, sd))
        names(stats) <- c(paste0(title, '_mean_count'), paste0(title, '_std_count'))
        stats
    }) %>%
        do.call(what = cbind)

    cbind(df, stats)
}


### generate_pca_metadata(id, path_metadata)
###
##  Summary
##      Make metadata files for pca ultimately used in  rnaseq::batch_pca
##
##  Args
##      id (character vector): samples names
##      path_metadata (character vector): directory for where to write file
##
##  Examples
##      id <- c('DU145_DMSO_Rep1','DU145_DMSO_Rep1','DU145_DMSO_Rep1')
##      generate_pca_metadata(id, path_metadata='./util)
##
get_tpm_stats <- function(txi, df, cols){

    tpm <- get_tpm_anno_df(txi)

    stats <- map(cols, ~ {

    }) %>% purrr::list_cbind()

    stats <- lapply(X = cols, FUN = function(x){
        title <- stringr::str_split(string = x[1], pattern = '__', simplify = TRUE) %>%
            '['(, -c(1,ncol(.)))
        x <- gsub(pattern = '__', replacement = '_', x = x)

        sub_tpm <- tpm %>% dplyr::select(all_of(x)) %>% 'row.names<-'(tpm$id)
        stats <- data.frame(rowMeans(sub_tpm), apply(sub_tpm, 1, sd))
        names(stats) <- c(paste0(title, '_mean_tpm'), paste0(title, '_std_tpm'))
        stats
    }) %>%
        do.call(what = cbind)

    cbind(df, stats)
}

### get_deg_list(outfile, infile, fc.seuil, p.adj)
###
##  Summary
##      Generate files for differentially expressed gene list depending on fold change
##      and p-value
##
##  Args
##      outfile (character): directory to write the differentially expressed gene list
##      infile (character): directory to get gene list
##      fc.seuil (numeric): threshold for the fold change
##      p.adj (numeric): threshold for the adjusted p-value
##
##  Examples
##      get_deg_list(outifle = './deg', infile = './gene_list', fc.seuil = 1.5, p.adj = 0.05)
##
get_deg_list <- function(outfile, infile, fc.seuil = NA, p.adj = NA){

    # Get all files (search is done recursively) from directory
    files <- fs::dir_ls(infile, recurse = TRUE, type = 'file')

    # Generate output directory and create the output file if not already done
    out.dir <- file.path(outfile, paste('FC',fc.seuil,'pAdj', p.adj, sep = '_'))
    if(!dir.exists(outfile)) dir.create(outfile, recursive = TRUE)

    # Set up parallel processing
    future::plan(multisession)

    # Process each files in parallel for the multicores machines
    furrr::future_walk(files, ~ {

        # Read and filter data
        id <- read_csv(file = .x, show_col_types = FALSE) %>%
            # Remove NA values
            tidyr::drop_na() %>%
            # Filter data for fold change and pValue given in the parameters
            dplyr::filter(
                (is.na(fc.seuil) | abs(log2FoldChange) > log2(fc.seuil)) &
                    (is.na(p.adj) | padj < p.adj)
            ) %>%
            # Retrieve associated ENSEMBL IDs
            dplyr::select(ensembl_gene)

        # Write .csv corresponding to selected data
        readr::write_csv(x = id, file = fs::path(out.dir, fs::path_file(.x)), col_names = FALSE)

        # Ensure reproductibility of the parallelize jobs
    }, .options = furrr_options(seed = TRUE))
}
