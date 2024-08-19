
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
