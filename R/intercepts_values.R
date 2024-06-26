library(tidyverse)
library(UpSetR)

### make_intercepts_list(anno, files_dir)
###
##  Summary
##      Compute list of genes subsetted for various fold changes
##
##  Args
##      anno (data.frame): annotation data.frame from the tximport function of the 'rnaseq' pipeline
##      files_dir (character): path where genes lists are stored
##
##  Examples
##      compute_foldchange(anno = txi$anno, files_dir = './data')
##
make_intercepts_list <- function(anno, files_dir){

    #### Hardcoded components ----
    cells <- c('DU145', 'VCaP', 'LAPC4')
    venn_file_path <- './dgea_output/venn_diagrams'
    venn_list_path <- file.path(venn_file_path, 'lists')
    venn_upset_path <- file.path(venn_file_path, 'upset')
    queried_list <- './genesMichele/genes_michele.csv'

    if(!dir.exists(venn_file_path)){dir.create(venn_file_path)}
    if(!dir.exists(venn_list_path)){dir.create(venn_list_path)}
    if(!dir.exists(venn_upset_path)){dir.create(venn_upset_path)}

    elem <- list(
        c("11OHT__DMSO", "R1881__DMSO", "11KT__DMSO"),
        c("11KT__DMSO",  "11OHT__DMSO"),
        c("11KT__DMSO",  "R1881__DMSO"),
        c("11OHT__DMSO", "R1881__DMSO"),
        c("11KT__DMSO"),
        c("11OHT__DMSO"),
        c("R1881__DMSO")
    )

    queried_genes <- subset(anno$id, anno$symbol %in% (
        readr::read_csv(file = queried_list, col_names = FALSE) %>%
            unlist %>% as.character)
        )

    #### Read all differentially expressed gene lists ----
    set_names <- list.files(files_dir) %>%
        sub(pattern = '.csv', replacement = '', x = .) %>%
        subset(x = ., subset = !grepl(pattern = 'ENZA', x = .) & grepl(pattern = 'DMSO', x = .))

    f <- stringr::str_extract(string = set_names, pattern = paste(cells, collapse = '|'))

    file_names <- file.path(venn_list_path, unique(f))
    for(fn in file_names){if(!dir.exists(fn)){dir.create(fn)}}

    files <- list.files(path = files_dir, pattern = paste(set_names, collapse = '|'), full.names = TRUE)

    lists <- lapply(X = files, FUN = readr::read_csv, col_names = FALSE) %>%
        setNames(set_names) %>%
        split(x = ., f = f)

    #### Generate intersection of list ----
    for(n in cells) {
        list <- lists[[n]]
        genes <- unlist(list) %>% unique

        m <- regexpr('__', names(list))
        names(list) <- regmatches(x = names(list), m = m, invert = TRUE) %>%
            sapply(FUN = '[[', 2)

        upset_list <- sapply(X = list, FUN = function(x) unlist(x) %>% as.character)
        upset_data <- fromList(input = upset_list) %>% 'rownames<-'(genes)
        #make_upsetPlot(upset_data = upset_data, path = venn_upset_path, name = n)

        for(i in seq_along(elem)){
            i <- elem[[i]]
            subname <- paste(i, collapse = '_')
            path <- file.path(venn_list_path, n, paste0(subname, '.csv'))
            extract_genes(positif = i, df = upset_data, anno = anno, path = path)
            }
    }
}

### extract_genes(positif, df, anno, path)
###
##  Summary
##      write files for various gene intersects
##
##  Args
##      positif:
##      df:
##      anno:
##      path:
##
##  Examples
##      extract_genes(out_path = './output', in_path = './data', fold.changes = c(1.5, 2))
##
extract_genes <- function(positif, df, anno, path){
    cond <- rowSums(df) == length(positif) &
        rowSums(sapply(positif, function(c) df[c] == 1)) == length(positif)

    goi <- subset(x = df, subset = cond) %>% rownames
    symbols <- anno[anno$id %in% goi, ]$symbol
    readr::write_csv(x = as.data.frame(symbols), file = path, col_names = FALSE)
}

### make_upsetPlot(upset_data)
###
##  Summary
##      Make upset plot for datasets
##
##  Args
##      upset_data (data.frame): representation of genes sets
##
##  Examples
##      data <- data.frame(11KT = C(1,1,1), 11OHT = c(1,0,0), R1881 = c(0,0,0)) %>%
##              'rownames<-'(c('ENSG00000103222', 'ENSG00000143554', 'ENSG00000198791'))
##      make_upsetPlot(upset_data = data)
##
make_upsetPlot <- function(upset_data, path, name){

    upset_path <- file.path(path, paste0(name, '.pdf'))
    pdf(file = upset_path)
    UpSetR::upset(
        data = upset_data,
        order.by = 'freq',
        nsets = length(list),
        empty.intersections = "on"
    )

    dev.off()

    my.subset <- subset(upset_data, genes %in% queried_genes)

    if(nrow(my.subset)){

        upset_path <- file.path(path, paste0(name,'_queried', '.pdf'))
        pdf(file = upset_path)
        UpSetR::upset(
            data = my.subset,
            order.by = 'freq',
            nsets = length(list),
            empty.intersections = "on"
        )

        dev.off()

    }
}
