library(tidyverse)
library(openxlsx)

genes.path <- './dataClementine/DEG_list.csv'
genes <- readr::read_csv(file = genes.path, col_names = FALSE) %>% as.matrix %>% as.character

# Get files, H vs DMSO, H vs HE
files.dir <- './dgea_output/de'
hormones <- c('R1881', '11KT', '11OHT', 'DMSO')
h.he <- sapply(hormones, function(h) paste(h,h, sep = '__'), USE.NAMES = FALSE)

deg.files <- list.files(path = files.dir, recursive = TRUE, full.names = TRUE)

deg.files.list <- list(
    deg.files %>% subset(grepl(pattern = 'DMSO.csv$', x = .) & !grepl(pattern = 'ENZA', x = .)),
    deg.files %>% subset(grepl(pattern = 'DMSO.csv$', x = .) & grepl(pattern = 'ENZA', x = .)),
    deg.files %>% subset(grepl(pattern = paste(h.he, collapse = '|'), x = .))
)
sheetName <- c('Hormone vs DMSO', 'Hormone + ENZA vs DMSO', 'Hormone vs Hormone + ENZA')

wb <- openxlsx::createWorkbook()
for (i in seq_along(sheetName)){

    openxlsx::addWorksheet(wb = wb , sheetName = sheetName[i])
    n.row <- 2
    n.col <- 2

    for(f in deg.files.list[[i]]){
        deg.infos <- readr::read_csv(file = f) %>% na.omit
        sub.deg <- subset(deg.infos, deg.infos$symbol %in% genes)
        comp.name <- stringr::str_split(string = f, pattern = '/')[[1]][5] %>%
            gsub(pattern = '.csv', replacement = '', x = .)

        openxlsx::writeData(wb = wb, sheet = sheetName[i],  x = comp.name, startCol = n.col, startRow = n.row)
        openxlsx::addStyle(wb = wb, sheet = sheetName[i], styl = openxlsx::createStyle(textDecoration = "bold"), rows = n.row, cols = n.col)
        n.row <- n.row + 1
        openxlsx::writeDataTable(wb = wb, sheet = sheetName[i], x = sub.deg, startCol = n.col, startRow = n.row)
        n.row <- n.row + nrow(sub.deg) + 3
    }
}

openxlsx::saveWorkbook(wb = wb, file = './dataClementine/concat_deg.xlsx', overwrite = TRUE)
