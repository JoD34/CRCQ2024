library(tidyverse)
library(ggplot2)
library(ggVennDiagram)
library(VennDiagram)
library(ggpubr)
library(RColorBrewer)
library(UpSetR)
library(openxlsx)
library(org.Hs.eg.db)
library(rlist)
library(magick)
library(pdftools)
library(magrittr)
library(rje)

source('./scripts/R/get_FCtable.R')
source('./scripts/R/util.R')

darker_color <- function(color, reduction){

    # Hex to rgb colors and darkening the color
    rgb <- grDevices::col2rgb(color) * reduction

    # Reverse to Hex colors
    grDevices::rgb(red = t(rgb), maxColorValue = 255)
}

### FUNCTION: triple_set_venn ----
#' @title Generate Venn Diagrams for Three Gene Sets

#' @description Generates Venn diagrams to compare three sets of genes identified
#'              from differentially expressed gene (DEG) analysis. It reads gene
#'              lists (CSV format) and creates separate Venn diagrams for each
#'              pairwise comparison between the sets. Finally, it combines the individual
#'              Venn diagrams into a single plot for visualization.

#' @return Creates a PDF file named "venn_hormone_DMSO.pdf" containing the combined Venn diagrams.

#' @examples
# (Assuming 'dgea_output/deg_list' directory exists with gene lists)
triple_set_venn <- function(output.dir, input.dir, cell, hormones, byCell = TRUE){

    if(!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)

    # Set Regex
    regex.cell <- paste(cell, collapse = '|')
    regex.hormones <- paste(hormones, collapse = '|')
    regex.files <- paste0('.*(', regex.cell,')__(',regex.hormones, ').*')

    # Generate comparisons for hormone vs hormones
    comp <- paste(hormones, hormones, sep = '__')

    # Get list of all differential expression analysis
    all.files <- list.files(path = input.dir, full.names = TRUE)

    # Get functions for subdividing sets
    sub.set <- subset_groups(df = all.files, comp = comp)

    # Set pattern for splitting files according to 'by' parameter
    group <- ifelse(byCell, '\\1', '\\2')

    # Subdivided all files into comparison groups
    my.files <- lapply(sub.set, function(FUN){
        files <- FUN(all.files)
        files <- subset(x = files, grepl(pattern = regex.files, x = files))
        split(x = files, f = sub('.*/(.*)__(.*)__.*\\.csv', group, files))
    }) %>% setNames(names(sub.set))


    # Set colors for sub groups of Venn (Hormone type)
    nCol <- ifelse(!byCell, length(cell), length(hormones))
    myCol <- brewer.pal(max(length(cell), length(hormones)), "Set1")[1:nCol]

    # Get a darker color for the subset's name
    darken.col <- darker_color(color = myCol, reduction = 0.7)

    # Iterate over all group comparison (ex. Hormone vs DMSO, Hormone vs Hormone & ENZA)
    lapply(X = names(my.files), FUN = function(group.name){

        tmp.outdir <- file.path(output.dir, group.name)
        if(!dir.exists(tmp.outdir)){dir.create(tmp.outdir)}

        group <- my.files[[group.name]]

        # Iterate over each cell lineage (eg DU145, LAPC4, VCaP)
        gene.list <- lapply(X = names(group), FUN = function(venn.group.name){

            # Get datasets of cell lineage
            venn.group <- group[[venn.group.name]]

            # Ordered the other way thant the earlier files
            group <- ifelse(byCell, '\\2', '\\1')
            venn.names <- sub('.*/(.*)__(.*)__.*', group, x = venn.group)

            # Load gene list for current comparison
            venn.list <- lapply(X = venn.group, FUN = function(f){
                readr::read_csv(file = f, col_names = FALSE) %>% unlist(use.names=FALSE)
            }) %>% setNames(venn.names)

            # Compute number of genes per hormones and per cell lineage
            num.genes <- paste0('(n = ', sapply(venn.list, length),')')
            genes.ligne <- unlist(venn.list) %>% unique %>% length

            # Generate the Venn diagrams
            VennDiagram::venn.diagram(
                x = venn.list,
                category.names = paste(names(venn.list), num.genes, sep = '\n'),
                filename = file.path(tmp.outdir,paste0(venn.group.name,'.png')),
                output = TRUE,

                # Output features
                imagetype = 'png',
                disable.logging = TRUE,
                main = paste0(venn.group.name,' (n = ',genes.ligne,')'),
                main.cex = 2,
                force.unique = FALSE,

                # Circle
                col = myCol,
                lwd = 2,
                lty = 1,
                fill = sapply(myCol, function(col) alpha(col, 0.3)),
                reverse = TRUE,

                # Numbers
                cex = 1.3,
                fontface = "bold",
                fontfamily = "sans",

                # Set names
                cat.default.pos = "outer",
                cat.fontfamily = "sans",
                cat.col = darken.col,
                cat.cex = 1.5
            )

            # Return list of gene to make tables
            venn.list
    }) %>%
            setNames(names(group))

        # Generate differentially expressed gene lists
        write.gList2xlsx(
            gene.list = gene.list,
            outdir = file.path(tmp.outdir, paste0(group.name, '_GeneLists.xlsx'))
            )

        # Concatenate Venn Diagrams for a given comparison (e.g. Hormone vs DMSO)
        prepare_venn_plots(dir = tmp.outdir, name = group.name, title)

    })
}

### write.gList2xlsx ----
#' Write Differentially Expressed Genes to Excel
#'
#' This function takes a list of differentially expressed genes and writes them to an Excel file. Each element in the list represents a sheet in the Excel workbook. Within each sheet, gene sets are organized into tables with merged cells for set names.
#'
#' @param gene.list A named list where each element is a list of gene sets. Each element of the outer list will be written to a separate sheet in the Excel workbook.
#' @param outdir The output directory path and filename for the Excel file.
#'
#' @return This function does not return a value. It writes an Excel file to the specified output directory.
#'
#' @details
#' The function creates a new Excel workbook and iterates over the elements of `gene.list`. For each element, it creates a new sheet in the workbook. Within each sheet, it iterates over the gene sets, writes the set names into merged cells, and writes the corresponding gene data into tables. The set names are center-aligned in the merged cells. The workbook is saved to the specified output directory at the end.
#'
#' @examples
#' \dontrun{
#' gene.list <- list(
#'   "Comparison1" = list(
#'     "Set1" = data.frame(Gene = c("GeneA", "GeneB"), FoldChange = c(2.0, -1.5)),
#'     "Set2" = data.frame(Gene = c("GeneC", "GeneD"), FoldChange = c(1.2, -2.3))
#'   ),
#'   "Comparison2" = list(
#'     "Set1" = data.frame(Gene = c("GeneE", "GeneF"), FoldChange = c(1.5, -2.1)),
#'     "Set2" = data.frame(Gene = c("GeneG", "GeneH"), FoldChange = c(2.4, -1.8))
#'   )
#' )
#'
#' write.gList2xlsx(gene.list, "output/degs.xlsx")
#' }
#'
#' @import openxlsx
#' @import dplyr
write.gList2xlsx <- function(gene.list, outdir){

    # Create Excel Workbook
    wb <- openxlsx::createWorkbook()
    sheets <- names(gene.list)

    lapply(X = sheets, FUN = function(sheet){

        # Create a new sheet for each cell lineage
        openxlsx::addWorksheet(wb = wb, sheetName = sheet)
        col.start <- 2

        # Get data corresponding to the cell lineage
        genes <- gene.list[[sheet]]

        # Write intersect into gene lists
        my.sets <- make_sets(lists = genes)

        for(i in names(my.sets)){

            df <- my.sets[[i]]

            # Manage excel formating for column's title
            openxlsx::mergeCells(wb, sheet, col.start:(col.start + 1), rows = 2)
            openxlsx::writeData(wb, sheet, i, xy = c(col.start, 2))
            openxlsx::addStyle(wb, sheet, openxlsx::createStyle(halign = 'center'), col = col.start, row = 2)
            openxlsx::writeDataTable(wb, sheet, df, startCol = col.start, startRow = 3)

            # Iterate the column so not to override one
            col.start <- col.start + ncol(df)
        }
    })

    openxlsx::saveWorkbook(wb = wb, file = outdir, overwrite = TRUE)
}

### make_sets ----
#' make_sets
#'
#' Function to create sets of genes based on input lists and retrieve corresponding gene symbols.
#'
#' This function takes a list of gene lists, generates binary matrices representing presence of genes,
#' creates sets of gene combinations, and retrieves gene symbols using org.Hs.eg.db.
#'
#' @param gene.lists A list of character vectors where each vector contains gene names.
#' @return A list where each element corresponds to a set of gene combinations with associated gene symbols.
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' gene_lists <- list(c("GeneA", "GeneB"), c("GeneB", "GeneC"))
#' sets <- make_sets(gene_lists)
#' }
#'
#' @references
#' Adapted from original function by [Author Name] in [Publication/Source].
#' @seealso
#' \code{\link{org.Hs.eg.db}}, \code{\link[stringr]{str_split}}, \code{\link[dplyr]{select}}
#' @keywords data manipulation
make_sets <- function(lists){

    # List all genes
    unique.genes <- lists %>% unlist(use.names = FALSE) %>% unique

    # Boolean dataframe presenting the set in which a gene is present
    df <- UpSetR::fromList(lists) %>% 'rownames<-'(unique.genes)

    # Generate every interesect's name
    sets <- rje::powerSet(x = colnames(df), m = ncol(df))[-1]
    sets.names <- sapply(X = sets, paste, collapse = '.')

    # Generate list of ENSEMBL IDs and corresponding Symbol for each intersects
    lapply(X = sets, FUN = function(check.cols){

        # Subset ENSEMBL IDs of genes that are part of the given intersect
        RowSums.equals <- df[rowSums(df[check.cols]) == rowSums(df), ]
        gene.int. <- rownames(RowSums.equals[rowSums(RowSums.equals) == length(check.cols), ])

        # Get ENSEMBL IDs and Symbols data.frame
        if(length(gene.int.) > 0) ensembl.symbols <- switchIds(ids = gene.int.)
        else {
            gene.int. <- NA
            ensembl.symbols <- NA
        }

        data.frame(ENSEMBL = gene.int., SYMBOLS = ensembl.symbols, row.names = NULL)

    }) %>% setNames(sets.names)
}

### FUNCTION : prepare_venn_plots ----
#' Prepare Venn Diagrams and Combine into a PDF
#'
#' This function reads PNG files from a directory, prepares Venn diagrams by converting them
#' to `ggplot` objects, and saves them as a combined PDF file with an appropriate title.
#'
#' @param dir A character string specifying the directory containing PNG files.
#' @param name A character string used to create the title of the combined PDF file.
#' The name is split by underscores, and abbreviations are replaced to form proper nouns.
#'
#' @return This function does not return a value. It creates and saves a PDF file
#' in the specified directory.
#'
#' @details This function performs the following steps:
#' \itemize{
#'   \item Reads PNG files from the specified directory.
#'   \item Converts the images to `ggplot` objects.
#'   \item Constructs a title by replacing abbreviations in the `name` parameter.
#'   \item Uses the `write_venns` function to save the combined images as a PDF.
#' }
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' prepare_venn_plots("path/to/directory", "H_E_combined")
#' }
#'
#' @import magick
#' @import ggplot2
#' @import ggpubr
#' @import stringr
#' @import grid
#'
prepare_venn_plots <- function(dir, name, title){

    # List files with the '.png" extension
    png.files <- list.files(path = dir, pattern = '\\.png$', full.names = TRUE)

    # Load images
    images <- lapply(X = png.files, FUN = image_read)

    # Change the title to proper nouns of abbreviations
    title <- stringr::str_split(string = name, pattern = '_', simplify = TRUE) %>%
        sapply(FUN = gsub, pattern = 'H', replacement = 'Hormone') %>%
        sapply(FUN = gsub, pattern = 'E', replacement = ' avec Enzalutamide') %>%
        paste(collapse = ' ')

    # Switch magick objects to grob object so ggarrange can be used to join graphs
    venns_ggplots <- lapply(X = images, FUN = function(img) {
        ggplot() +
            annotation_custom(grob = rasterGrob(image = as.raster(x = img), interpolate = TRUE)) +
            theme_void() +
            theme(panel.background = element_rect(fill = "transparent", colour = NA))
    })

    # Send graph list to print pdf file
    write_venns(
        venns = venns_ggplots,
        title = bquote(atop(bold(.(title)), phantom(0) * "p.adj < 0.05 ; |FC| >" ~ log[2](1.3))),
        filename = file.path(dir, paste(name, "combined.pdf", sep = '_'))
    )
}

### FUNCTION: double_set_venn ----
#' @title Create Venn Diagrams for Differential Expression Analysis

#' @description
#' This function processes differential expression (DE) results files
#' (expected to be named "DMSO.csv" and contain a column named "ensembl_gene")
#' to generate Venn diagrams comparing gene expression between hormone-treated
#' and control (DMSO) samples for different cell lines.

#' @return This function invisibly returns a list containing the created
#' Venn diagrams, but it saves the resulting plots as PNG images in the specified output directory.

#' @examples
#' # Assuming your DE results files are located as expected by the function
#' double_set_venn()
double_set_venn <- function(){
    cells <-  c('DU145', 'LAPC4', 'VCaP')
    hormones <- c('R1881', '11KT', '11OHT')

    outdir <- './dgea_output/venn_diagrams/HvsDMSO_venn_HEvsDMSO'
    if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
    sub.outidir <- file.path(outdir, cells) %>% setNames(cells)
    for (i in sub.outidir) if(!dir.exists(i)) dir.create(i)

    file_set <- list.files(path = './dgea_output/de', pattern = 'DMSO.csv$', recursive = TRUE, full.names = TRUE) %>%
        split(., f = stringr::str_extract(., pattern = paste(cells, collapse = '|'))) %>%
        map(~ split(., f = stringr::str_extract(., pattern = paste(hormones, collapse = '|'))))

    lapply(X = names(file_set), FUN = function(name_set){

        set <- file_set[[name_set]]
        lapply(X = set, function(f) {
            file.name <-stringr::str_extract(string = f, pattern = '[A-z0-9]+__[A-z0-9]+')
            loaded <- lapply(f, function(df) readr::read_csv(df) %>% na.omit %>% dplyr::pull(ensembl_gene))

            output <- paste(file.name, collapse = '__') %>%
                paste('Venn',paste(file.name, collapse = '__'), sep = '__') %>% paste0(., '.png') %>%
                file.path(outdir, name_set,.)

            venn.diagram(
                x = loaded,
                category.names = file.name,
                filename = output,
                cat.pos = 0,
                output = TRUE,
                rotation= 1)
            })
    })

}

### FUNCTION: Make_2sets_ARp ----
#' @title Create Two-way Venn Diagrams for AR Positive Cell Lines

#' @description Processes differential expression (DE) results files (expected
#'              to be named "DMSO.csv" and contain a column named "ensembl_gene")
#'              to generate two-way Venn diagrams comparing AR positive cell
#'              lines for each hormone treatment.

#' @return Returns a list containing the created Venn diagrams, but it saves the
#'         resulting plots as PNG images in the specified output directory.

#' @examples
#' # Assuming your DE results files are located as expected by the function
#' Make_2ways_ARp()
Make_2sets_ARp <- function(){

    cells <- c('LAPC4', 'VCaP')
    hormones <- c('R1881', '11KT', '11OHT')
    indir <- c('./dgea_output/deg_list')

    outdir <- c('./dgea_output/venn_diagrams/Hormones_ARp')
    if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

    used.files <- list.files(path = indir, pattern = 'DMSO.csv', recursive = TRUE, full.name = TRUE) %>%
        '['(grep(pattern = 'DU145|ENZA', x = ., invert = TRUE))

    split.files <- split(
        x = used.files,
        f = stringr::str_extract(string = used.files, pattern = paste(hormones, collapse = '|'))
    )

    wb <- createWorkbook()
    myCol <- brewer.pal(3, "Set1")[-1]

    for(group.name in names(split.files)){

        group <- split.files[[group.name]]
        loaded <- lapply(X = group, FUN = function(file) {
            readr::read_csv(file = file, col_names = FALSE) %>% na.omit %>%
                as.matrix %>% as.character
        }) %>%
            setNames(stringr::str_extract(string = group, pattern = paste(cells, collapse = '|')))

        total.pop <- unlist(loaded) %>% as.character %>% unique %>%  length
        num.genes <- paste0(' (n = ', sapply(loaded, length, USE.NAMES = FALSE), ')')

        # Make Venn Diagrams
        VennDiagram::venn.diagram(
            x = loaded,
            category.names = paste(names(loaded), num.genes, sep = '\n'),
            filename = file.path(outdir, paste0(group.name, '.png')),
            output = TRUE,

            # Output features
            imagetype = 'png',
            disable.logging = TRUE,
            main = paste0(group.name, ' (n = ',total.pop,')'),
            main.cex = 2,
            sub.cex = 1.5,
            force.unique = FALSE,

            # Circle
            col = myCol,
            lwd = 2,
            lty = 1,
            fill = sapply(myCol, function(col) alpha(col, 0.3)),
            reverse = TRUE,

            # Numbers
            cex = 1.3,
            fontface = "bold",
            fontfamily = "sans",

            # Set names
            cat.default.pos = "outer",
            cat.fontfamily = "sans",
            cat.col = myCol,           # Color of category name (matching the circle's)
            cat.pos = 200 ,            # Position of category name (degree arount circle)
            cat.dist = 0.05,           # Distance of category name from circl
            cat.cex = 1.5,             # Size of category name writting

            # Perform hypergeometric test for 2set Venn Diagram
            hyper.test = TRUE,
            total.population = total.pop
        )

        append_workbook_ensembl(wb = wb, gene.list = loaded, group.name = group.name)
    }
    openxlsx::saveWorkbook(wb = wb, file = file.path(outdir,'gene_list_PerHormone_ARp.xlsx'), overwrite = TRUE)
    prepare_2set_venn_plots(dir = outdir, filename = file.path(outdir, 'concat.pdf'))
}


prepare_2set_venn_plots <- function(dir, filename){
    png.files <- list.files(path = dir, pattern = '\\.png$', full.names = TRUE)
    images <- lapply(X = png.files, FUN = image_read)

    # Switch magick objects to grob object so ggarrange can be used to join graphs
    venns_ggplots <- lapply(X = images, FUN = function(img) {
        ggplot() +
            annotation_custom(grob = rasterGrob(image = as.raster(x = img), interpolate = TRUE)) +
            theme_void() +
            theme(panel.background = element_rect(fill = "transparent", colour = NA))
    })

    # Send graph list to print pdf file
    write_venns(
        venns = venns_ggplots,
        title = "Comparaison: lignées AR+",
        filename = filename
    )
}


### FUNCTION: write_venns ----
#' @title Create a combined Venn diagram plot and save it to a PDF

#' @description
#' This function takes a list of Venn diagrams (`venns`) and a title (`title`) as input.
#' It arranges the Venn diagrams horizontally with a slight spacing between them and saves the resulting plot as a PDF.

#' @param venns A list containing individual Venn diagrams created using packages like `venn` or `ggvenn`.

#' @param title The title to be displayed at the top of the combined plot.

#' @return This function invisibly returns the arranged plot object (`v`) but saves the plot as a PDF file.

#' @examples
#' library(venn) # Or any other library you used to create the Venn diagrams in 'venns'

#' # Assuming you have created Venn diagrams and stored them in a list named 'venns'
#' write_venns(venns = venns, title = "My Venn Diagram Comparison")
write_venns <- function(venns, title, filename){

    null.list <- rep(list(NULL), length(venns))
    venns_list <- mapply(list, venns, null.list)
    dist <- c(0.75, 0.05)
    my.widths <- rep(dist, length.out = length(venns) * length(dist) - 1)

    v <- ggarrange(plotlist = venns_list, ncol = 5, nrow = 1, widths = my.widths, heights = 2)
    annotate_figure(p = v, top = text_grob(title, face = "bold", size = 14))

    ggsave(
        filename = filename,
        height = unit(x = 5, units = 'inch'),
        width = unit(x = 11, units = 'inch')
    )
}

### Function: append_workbook_ensembl ----
#' @title Append Worksheet to Existing Workbook for Upset Analysis Results

#' @description
#' This function takes an existing Openxlsx workbook (`wb`), a gene list (`gene.list`),
#' and a group name (`group.name`) as input. It processes the gene list to identify genes
#' specific to each cell line (LAPC4 or VCaP) or common to both, and then appends a new worksheet
#' to the workbook containing these results.

#' @param wb An existing Openxlsx workbook object.

#' @param gene.list A list containing gene expression data for each cell line (LAPC4 and VCaP).
#' The list elements should be named "LAPC4" and "VCaP", each containing character vectors of ENSEMBL IDs.

#' @param group.name A character string specifying the name for the new worksheet in the workbook.

#' @return This function invisibly returns the modified workbook (`wb`) with the appended worksheet.

#' @examples
#' library(openxlsx)
#' # Assuming you have a workbook 'my_workbook.xlsx' and a gene list 'my_genes'
#' wb <- loadWorkbook("my_workbook.xlsx")
#' append_workbook(wb = wb, gene.list = my_genes, group.name = "Upset_Results")
#' saveWorkbook(wb, "my_workbook_with_upset.xlsx")
append_workbook_ensembl <- function(wb, gene.list, group.name){

    # Organize data ----
    ensembl.id <- unlist(gene.list) %>% unique
    upset.data <- data.frame(
        LAPC4 = as.numeric(ensembl.id %in% gene.list$LAPC4),
        VCaP = as.numeric(ensembl.id %in% gene.list$VCaP)
    ) %>%
        'rownames<-'(ensembl.id)

    # Subset data depending on set ----
    my.lists <- list(
        subset(upset.data, rowSums(upset.data) == 1 & upset.data$LAPC4 == 1),
        subset(upset.data, rowSums(upset.data) == 1 & upset.data$VCaP == 1),
        subset(upset.data, rowSums(upset.data) == 2)
    ) %>%
        setNames(c('LAPC4_Only', 'VCaP_Only', 'Both'))

    # Write Excel file ----
    start.col <- 2
    openxlsx::addWorksheet(wb = wb, sheetName = group.name)
    for(lname in names(my.lists)){

        # Organize data
        l <- my.lists[[lname]]
        info <- select(org.Hs.eg.db, keys = rownames(l), keytype = 'ENSEMBL', columns = c('SYMBOL'))
        x <- data.frame(info) %>% setNames(paste0(lname, c('.ENSEMBL','.SYMBOL')))

        # Write data tables
        openxlsx::writeDataTable(wb = wb , sheet = group.name, x = x, startCol = start.col, startRow = 2)
        start.col <- start.col + 3
    }
}
