library(tidyverse)
library(RColorBrewer)
library(rnaseq)

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
generate_pca_files <- function(samples, out_meta, out_infos){
    generate_pca_metadata(id = samples, output = out_meta)
    generate_pca_infos(output = out_infos)
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
generate_info_txi <- function(txi, path){
    subfunc <- function(FUN){

        elem <- FUN(txi = txi)
        rownames(elem) <- elem$ensembl_gene
        subelem <- elem %>% dplyr::select(-colnames(txi$anno))

        new.coln <- colnames(subelem) %>%
            sub(pattern = '_Rep[0-9]', replacement = '', x = .) %>% unique

        extracted <- stringr::str_extract(string = colnames(subelem), pattern = paste(new.coln, collapse = '|'))
        elem.split <- lapply(X = new.coln, FUN = function(coln) {
            subset(subelem, select = (extracted == coln))
            }) %>% setNames(new.coln)

        stat <- lapply(X = elem.split, FUN = function(set) {
            data.frame(rowMeans(set), apply(set, 1, sd)) %>%
            'colnames<-'(c('mean', 'std'))
        }) %>% do.call(what = cbind)

        stat
    }

    format_prep <- function(df, count){
        add_str <- ifelse(count, 'count', 'tpm')
        colnames(df) <- sub(
            pattern = '.mean',
            replacement = paste('_mean', add_str, sep = '_'),
            x = colnames(df)
        )
        colnames(df) <- sub(
            pattern = '.std',
            replacement = paste('_std', add_str, sep = '_'),
            x = colnames(df)
        )
        count <- df %>% dplyr::mutate(id = rownames(df))
    }

    cells <- c('DU145', 'LAPC4','VCaP')

    count <- subfunc(FUN = get_raw_count_anno_df) %>%
        format_prep(., TRUE)

    tpm <- subfunc(FUN = get_tpm_anno_df) %>%
        format_prep(., FALSE)

    df <- plyr::join(x = count, y = tpm, by = 'id', type = "full")

    sets <- lapply(X = cells, FUN = function(c){
        sub.df <- subset(df, select = grepl(pattern = c, colnames(df))) %>%
            dplyr::mutate(id = df$id)
        colnames(sub.df) <- sub(pattern = paste0(c,'_'), replacement = '', x = colnames(sub.df))
        sub.df <- sub.df %>%
            plyr::join(x = ., y = txi$anno, by = 'id', type = 'left') %>%
            dplyr::select(-c('id', 'transcript_type')) %>%
            dplyr::relocate(c('ensembl_gene', 'symbol', 'entrez_id'), .before = everything())
        readr::write_csv(x = sub.df, file = file.path(path, paste0(c,'.csv')))
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
    de_files <- list.files(path_de, recursive = TRUE, full.names = TRUE)

    # Get differential expression files' id
    de <- stringr::str_split(
        string = de_files,
        pattern = '/',
        simplify = TRUE)[, 5] %>%
        sub(pattern = '.csv', replacement = '', x = .)

    de_results <- lapply(X = de_files, FUN = read.csv)
    names(de_results) <- de

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
    lapply(X = list.files(path, full.names = TRUE), FUN = function(p){
        df <- read.csv(p) %>% dplyr::select(-c('transcript_type', 'id'))
        fold_change <- 2^df$log2FoldChange
        df <- df %>% dplyr::mutate(fold_change, .after = log2FoldChange)
        df <- add_tmp_counts(df = df, txi = txi, path = p)
        readr::write_csv(x = df, file = p)
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

    files <- list.files(infile, full.names = TRUE, recursive = TRUE)

    lapply(X = files , FUN = function(p) {
        filename <- stringr::str_split(string = p, pattern = '/')[[1]][5]

        id <- readr::read_csv(file = p) %>% na.omit

        if (!is.na(fc.seuil)){id <- subset(id, abs(id$fold_change) > fc.seuil)}
        if (!is.na(p.adj)){id <- subset(id, id$padj < p.adj)}

        file <- file.path(outfile, filename)
        id <- subset(id, select = 'ensembl_gene')
        readr::write_csv(x = id, file = file, col_names = FALSE)
        })
}
