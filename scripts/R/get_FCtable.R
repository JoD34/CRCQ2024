
#' Compute number of significant genes based on fold change.
#'
#' @description Reads gene expression data from CSV files located in a specified directory,
#'              filters the data based on fold change values for different conditions,
#'              and computes the number of unique genes that meet the filtering criteria for each condition.
#'
#' @param in_path Character string specifying the directory path where the CSV files are located.
#' @param fold.changes Numeric vector specifying fold change values used as thresholds for filtering genes.
#'
#' @return A data frame where each row corresponds to a fold change value in \code{fold.changes}.
#'         Each column corresponds to a condition, defined by combinations of cell types ('DU145', 'LAPC4', 'VCaP')
#'         and treatments ('11KT__11KT', 'R1881__R1881', '11OHT__11OHT'). Each cell contains the count of unique genes
#'         meeting the filtering criteria for the respective condition.
#'
#' @examples
#' # Define example inputs
#' in_path <- "/path/to/your/files/"
#' fold_changes <- c(1.5, 2.0, 3.0)
#'
#' # Load required libraries
#' #' @importFrom tidyverse %>%
#' #' @importFrom stringr %>%
#' #' @importFrom readr %>%
#' #' @importFrom purrr %>%
#' #' @importFrom dplyr %>%
#' library(tidyverse)
#' library(stringr)
#' library(readr)
#' library(purrr)
#' library(dplyr)
#'
#' # Compute fold changes
#' results <- compute_foldchange(in_path, fold_changes)
#'
#' @export
compute_foldchange <- function(in_path, fold.changes){
    cells <- c('DU145', 'LAPC4', 'VCaP')
    comp <- c('11KT__11KT', 'R1881__R1881', '11OHT__11OHT')

    to.do <- c(
        function(df, fold.changes) process_file(df, fold.changes, 'up'),
        function(df, fold.changes) process_file(df, fold.changes, 'down')
        #function(df, fold.changes) process_file(df, fold.changes, 'total')
    )

    all.files <- list.files(path = in_path ,recursive = TRUE, full.names = TRUE)
    sub.set <- subset_groups(df = all.files, comp = comp)
    my.files <- lapply(sub.set, function(FUN){
        FUN(all.files) %>%
        split(x = ., f = stringr::str_extract(string = ., pattern = paste(cells, collapse = '|')))
    })

    use.files <- unlist(my.files) %>%
        lapply(readr::read_csv) %>%
        setNames(unlist(my.files))

    result <- lapply(X = my.files, function(file.set){

        lapply(to.do, function(FUN){

            sapply(X = file.set , FUN = function(group.files) { # Select a cell type

                group.results <- lapply(use.files[group.files], FUN, fold.changes = fold.changes)

                apply(simplify2array(group.results), MARGIN = 1, FUN = function(row) {
                    unlist(row) %>% unique %>% length
                })

            }) %>%
                as.data.frame

        }) %>%
            setNames(c('Up_regulated', 'Down_regulated'))

    }) %>%
        setNames(c('H_vs_DMSO', 'HE_vs_DMSO', 'HE_VS_DMSOE', 'H_VS_HE'))

    df <- lapply(result, function(i){
        lapply(cells, function(c){

            data.frame(
                subset(i$Down_regulated, select = c),
                subset(i$Up_regulated, select = c)
                ) %>% 'colnames<-'(paste(c, c('Down','Up'), sep = '.'))

            }) %>%
            do.call(what = cbind) %>%
            dplyr::mutate('Fold Change' = ifelse(is.na(fold.changes), 'NA', fold.changes)) %>%
            dplyr::relocate('Fold Change', .before = dplyr::everything())
        })

    # Write Excel file
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb = wb, sheetName = 'Data')
    excel <- write_excel(df = df, wb = wb, sheetName = 'Data')
    openxlsx::saveWorkbook(wb = excel, file = "Quantification_FC.xlsx", overwrite = TRUE)
}

#' Process a CSV file to filter genes based on fold changes.
#'
#' @description Reads a CSV file containing gene data, filters the genes based
#'              on specified fold changes and type of filtering ('up', 'down', 'total'),
#'              and returns a list of gene names that meet the filtering criteria
#'              for each fold change value.
#'
#' @param file_paths Character vector of file paths to CSV files containing gene data.
#' @param fold.changes Numeric vector of fold change values used as thresholds
#'                     for filtering genes.
#' @param type Character string indicating the type of filtering ('up', 'down', 'total').
#'
#' @return A list where each element corresponds to a fold change value in
#'         \code{fold.changes}. Each element is itself a character vector of gene names
#'         that meet the filtering criteria.
#'
#' @examples
#' # Define example inputs
#' files <- c("file1.csv", "file2.csv")
#' fold_changes <- c(1.5, 2.0, 3.0)
#' type <- "up"
#'
#' # Load required libraries
#' #' @importFrom readr %>%
#' #' @importFrom dplyr %>%
#' library(readr)
#' library(dplyr)
#'
#' # Process files
#' results <- process_file(files, fold_changes, type)
#'
#' @export
process_file <- function(df, fold.changes, type){

    purrr::map(fold.changes, ~ {

        fc <- switch(
            type,
            'up'    = ifelse(is.na(.x), 0, log2(.x)),
            'down'  = ifelse(is.na(.x), 0, -log2(.x)),
            'total' = ifelse(is.na(.x), 0, abs(log2(.x)))
        )

        df %>%
            dplyr::filter(
                case_when(

                    type == 'up'    ~ log2FoldChange > fc,
                    type == 'down'  ~ log2FoldChange < fc,
                    type == 'total' ~ abs(log2FoldChange) > fc

                    ) & padj < 0.05
                ) %>%
            dplyr::pull(ensembl_gene)
    })
}

#' Write Data Frames to an Excel Workbook
#'
#' @description This function creates an Excel workbook, writes multiple data frames to it on a single sheet,
#'              adds titles for each data frame, and formats these titles in bold.
#'
#' @param df A named list of data frames. Each element in the list should be a data frame,
#'           and each name of the list will be used as the title for the corresponding data frame in the Excel sheet.
#'
#' @return This function does not return a value. It creates and saves an Excel file named "Quantification_FC.xlsx"
#'         in the working directory.
#'
#' @examples
#' # Load required library
#' library(openxlsx)
#'
#' # Define example data frames
#' df1 <- data.frame(A = 1:5, B = letters[1:5])
#' df2 <- data.frame(X = 5:1, Y = LETTERS[1:5])
#'
#' # Combine into a named list
#' df_list <- list("First DataFrame" = df1, "Second DataFrame" = df2)
#'
#' # Write data frames to Excel
#' write_excel(df_list)
#'
#' @importFrom openxlsx createWorkbook addWorksheet writeData addStyle writeDataTable saveWorkbook createStyle
#' @export
write_excel <- function(df, wb, sheetName){

    rnum <- 1L
    ctitle <- 5L
    tab.col <- 2L
    bold <- openxlsx::createStyle(textDecoration = "bold")

    for (i in names(df)){
        openxlsx::writeData(wb = wb, sheet = sheetName, x = i, startCol = ctitle, startRow = rnum)
        openxlsx::addStyle(wb, sheet = sheetName, style = bold, rows = rnum, cols = ctitle, gridExpand = TRUE)
        rnum <- rnum + 1

        openxlsx::writeDataTable(
            wb = wb, sheet = sheetName, x = df[[i]], startCol = tab.col, startRow = rnum
            )
        rnum <- rnum + nrow(df[[i]]) +2
        }

    wb
}

#' @title Filter data frame based on two patterns
#' @description Filters a character vector containing a first pattern and
#'              optionally excludes rows containing pattern_2.
#'
#' @param df Data frame to be filtered.
#' @param pattern_1 Character vector specifying the pattern(s) to include.
#' @param pattern_2 Character vector specifying the pattern(s) to exclude
#' @param invert Logical value indicating if rows containing pattern_2 should be
#' kept (FALSE, default) or excluded (TRUE).
#'
#' @return A data frame containing rows matching the filtering criteria.
#'
#' @examples
#' # Example data frame
#' df <- data.frame(x = c("apple", "banana", "cherry", "grape"),
#'                   y = c(1, 2, 3, 4))
#'
#' # Include rows with "apple" and exclude rows with "3"
#' filtered_df <- filter_data(df, pattern_1 = "apple", pattern_2 = "3", invert = TRUE)
#' print(filtered_df)
#'
#' # Include rows with any digit and exclude rows with "apple"
#' filtered_df <- filter_data(df, pattern_1 = "\\d+", pattern_2 = "apple", invert = TRUE)
#' print(filtered_df)
filter_data <- function(df, pattern_1, pattern_2, invert){
    df[grep(pattern = pattern_1, x = df)] %>%
        '['(grep(pattern = pattern_2, x = ., invert = invert))
}

#' @title Subset data frame by group based on file names
#' @description Filter the character vector based on specific groups
#'             defined by file name patterns.
#'
#' The groups are defined by the following logic:
#'  * Group 1: Files finishing with "DMSO.csv" AND excluding "ENZA"
#'  * Group 2: Files finishing with "DMSO.csv" AND "ENZA"
#'  * Group 3: Files finishing "DMSO_ENZA.csv" AND excluding  "__DMSO__"
#'  * Group 4: Files containing any combination of characters specified in the
#'  `comp` argument (separated by '|')
#'
#' @param df Character vector containing file paths
#' @param comp Character vector specifying patterns to include in Group 4
#'
#' @return A list of functions. Each function can be used to filter the
#'  data frame for a specific group.
#'
#' @examples
#' # Example data frame (assuming 'filepath' column contains file paths)
#' df <- data.frame(filepath = c("data/DMSO.csv", "data/control.txt",
#'                               "data/DMSO_ENZA.csv", "data/ENZA_data.csv"))
#'
#' # Get filtering functions for each group
#' filter_funcs <- subset_groups(df)
#'
#' # Filter data for Group 1 (DMSO.csv but not ENZA)
#' group1_data <- filter_funcs[[1]](df)
#' print(group1_data)
#'
#' # Filter data for Group 4 (any combination from 'comp' argument)
#' comp <- c("treatment", "control")  # Example patterns
#' filter_funcs <- subset_groups(df, comp = comp)
#' group4_data <- filter_funcs[[4]](df)
#' print(group4_cast)  # Typo corrected: cast -> data
subset_groups <- function(df, comp){
    list(
       H_vs_DMSO = function(df) filter_data(df = df, pattern_1 = '.*DMSO.csv$', pattern_2 = 'ENZA', invert = TRUE),
       HE_vs_DMSO = function(df) filter_data(df = df, pattern_1 = '.*DMSO.csv$', pattern_2 = 'ENZA', invert = FALSE),
       HE_VS_DMSOE = function(df) filter_data(df = df, pattern_1 = '.*DMSO_ENZA.csv$', pattern_2 = '__DMSO__', invert = TRUE),
       H_VS_HE = function(df) filter_data(df = df, pattern_1 = paste(comp, collapse = '|'), pattern_2 = '', invert = FALSE)
    )
}

#' quant_express_distr
#'
#' Function to process and analyze quantitative expression data from multiple files.
#'
#' This function processes files located in a specified directory, separates them based on specified comparisons,
#' calculates fold changes for up- and down-regulated genes, and outputs the results into an Excel file.
#'
#' @param in_path Character string specifying the path to the directory containing input files.
#' @param fold.changes Numeric value specifying the fold change threshold for gene expression changes.
#' @return No direct return; outputs an Excel file named 'quant_UpDown_hormone.xlsx' containing processed dataframes.
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' quant_express_distr(in_path = "/path/to/your/files", fold.changes = 2)
#' }
#' @references
#' Adapted from original function by [Author Name] in [Publication/Source].
#' @seealso
#' \code{\link{process_file}}, \code{\link{subset_groups}}
#' @keywords data processing
quant_express_distr <- function(in_path, fold.changes){
    cells <- c('DU145', 'LAPC4', 'VCaP')
    hormone <- c('11KT', '11OHT', 'R1881')
    comp <- sapply(hormone, function(h) paste(h, h, sep = '__'))

    to.look <- c(
        up = function(df, fold.changes) process_file(df, fold.changes, 'up'),
        down = function(df, fold.changes) process_file(df, fold.changes, 'down')
    )

    all.files <- list.files(path = in_path ,recursive = TRUE, full.names = TRUE)
    sub.set <- subset_groups(df = all.files, comp = comp) # Get functions to separate files
    my.files <- lapply(sub.set, function(FUN){ # Separate by comparison and cells types
        FUN(all.files) %>%
            split(x = ., f = stringr::str_extract(string = ., pattern = paste(cells, collapse = '|')))
    }) %>%
        setNames(c('H_vs_DMSO', 'HE_vs_DMSO', 'HE_VS_DMSOE', 'H_VS_HE'))

    use.files <- unlist(my.files) %>%
        lapply(readr::read_csv) %>%
        setNames(unlist(my.files))

    res <- lapply(X = my.files, function(file.set){

        lapply(X = file.set , FUN = function(group.files) { # Select a cell type

            lapply(use.files[group.files], function(df){
                up <- to.look$up(df = df, fold.changes = fold.changes) %>% sapply(FUN = length)
                down <- to.look$down(df =df , fold.changes = fold.changes) %>% sapply(FUN = length)
                data.frame(Down_reg = down, Up_reg = up)
            }) %>%
                setNames(stringr::str_extract(string = names(.), pattern = paste(hormone, collapse = '|'))) %>%
                do.call(what = cbind) %>%
                dplyr::mutate('Fold Change' = ifelse(is.na(fold.changes), 'NA', fold.changes)) %>%
                dplyr::relocate('Fold Change', .before = everything())
        })
    })

    # Write Excel file with dataframes
    wb <- openxlsx::createWorkbook()
    lapply(names(res), FUN = function(comparisons){
        openxlsx::addWorksheet(wb = wb, sheetName = comparisons)
        write_excel(df = res[[comparisons]], wb = wb, sheetName = comparisons)
        })
    openxlsx::saveWorkbook(wb = wb, file = 'quant_UpDown_hormone.xlsx', overwrite = TRUE)
}
