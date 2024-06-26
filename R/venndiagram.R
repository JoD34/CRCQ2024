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
source('./R/get_FCtable.R')

### FUNCTION: make_venn ----
#' @title Create a Venn Diagram

#' @description This function creates a Venn diagram to visualize the overlap between
#' multiple sets of data. It expects a list named `plot` where each element
#' represents a set of data (e.g., character vector of gene IDs).

#' @param plot A list containing character vectors, where each element represents a set of data
#' to be included in the Venn diagram.

#' @param title A character string specifying the title for the Venn diagram.

#' @return A ggplot object representing the Venn diagram.

#' @examples
#   Example plot data (replace with your actual data)
#'      set1 <- c("A", "B", "C", "D")
#'      set2 <- c("B", "C", "E", "F")
#'      set3 <- c("C", "D", "G", "H")
#'      plot_data <- list(set1 = set1, set2 = set2, set3 = set3)

#   Generate Venn diagram
#'      venn_plot <- make_venn(plot = plot_data, title = "Gene Overlap")
#'      print(venn_plot)
make_venn <- function(plot, title){
    tot <- unlist(plot) %>% unique() %>% length()
    categories <- stringr::str_split(string = names(plot), pattern = '_', simplify = TRUE)[,3] %>%
        paste0(.,c('\n','\n',''))
    nb <- sapply(plot, function(x) length(x))

    ggVennDiagram(x=plot, category.names = paste0(categories, 'n=', nb),
                  label_alpha = 0, label='count', label_size = 4,
                  edge_size=1, set_color = c("#FFD700","#228B22","#99000d"),) +
        scale_fill_distiller(palette = "Blues", direction = 1) +
            labs(title = paste0(title, ' (n = ', tot, ')')) +
        theme(legend.position = "none")
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
triple_set_venn <- function(){

    cell <- c('DU145', 'VCaP', 'LAPC4')
    hormones <- c('R1881', '11KT', '11OHT')
    comp <- sapply(hormones, function(h) paste(h, h, sep = '__'), USE.NAMES = FALSE)

    outdir <- c('./dgea_output/venn_diagrams/all_venn')
    if(!dir.exists(outdir)){dir.create(outdir)}

    # Get list of all differential expression analysis
    all.files <- list.files(path = './dgea_output/deg_list', full.names = TRUE)

    # Get functions for subdividing sets
    sub.set <- subset_groups(df = all.files, comp = comp)

    # Subdivided all files into comparison groups
    my.files <- lapply(sub.set, function(FUN){
        FUN(all.files) %>%
            split(x = ., f = stringr::str_extract(string = ., pattern = paste(cell, collapse = '|')))
    }) %>%
        setNames(c('H_vs_DMSO', 'HE_vs_DMSO', 'HE_VS_DMSOE', 'H_VS_HE'))

    # Set colors for sub groups of Venn (Hormone type)
    myCol <- brewer.pal(3, "Set1")

    # Iterate over all group comparison (ex. Hormone vs DMSO, Hormone vs Hormone & ENZA)
    lapply(X = names(my.files), FUN = function(group.name){

        tmp.outdir <- file.path(outdir, group.name)
        if(!dir.exists(tmp.outdir)){dir.create(tmp.outdir)}

        group <- my.files[[group.name]]

        # Iterate over each hormone group
        gene.list <- lapply(X = names(group), FUN = function(venn.group.name){
            venn.group <- group[[venn.group.name]]
            venn.names <- stringr::str_extract(string = venn.group, pattern = paste(hormones, collapse = '|')) %>% unlist

            # Load gene list for current comparison
            venn.list <- lapply(X = venn.group, FUN = function(f){
                readr::read_csv(file = f, col_names = FALSE) %>% as.matrix %>% as.character
            }) %>%
                setNames(venn.names)

            # Generate the Venn diagrams
            VennDiagram::venn.diagram(
                x = venn.list,
                filename = file.path(tmp.outdir,paste0(venn.group.name,'.png')),
                output = TRUE,

                # Output features
                imagetype = 'png',
                disable.logging = TRUE,
                main = paste('Lignée ', venn.group.name),
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
                cat.col = myCol,
                cat.cex = 1.5
            )

            # Return list of gene to make tables
            venn.list
    }) %>%
            setNames(names(group)) %>%
            write_deg2excel(     # Generate differentially expressed gene lists
                gene.list = .,
                outdir = file.path(tmp.outdir, paste(group.name, 'GeneLists.xlsx', sep = '_'))
            )

        # Concatenate Venn Diagrams for a given comparison (e.g. Hormone vs DMSO)
        prepare_venn_plots(dir = tmp.outdir, name = group.name)

    })
}

### write_deg2excel ----
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
#' write_deg2excel(gene.list, "output/degs.xlsx")
#' }
#'
#' @import openxlsx
#' @import dplyr
write_deg2excel <- function(gene.list, outdir){
    # Create Workbook for Excel file
    wb <- openxlsx::createWorkbook()

    lapply(X = names(gene.list), FUN = function(sheetName){

        openxlsx::addWorksheet(wb = wb, sheetName = sheetName)
        col.start <- 2

        genes <- gene.list[[sheetName]]
        my.sets <- make_sets(gene.lists = genes)

        for(i in names(my.sets)){

            df <- my.sets[[i]]

            openxlsx::mergeCells(wb = wb, sheet = sheetName, cols = col.start:(col.start + 1), rows = 2)
            openxlsx::writeData(wb = wb, sheet = sheetName, x = i, xy = c(col.start, 2))
            openxlsx::addStyle(wb = wb, sheet = sheetName, style = openxlsx::createStyle(halign = 'center'), col = col.start, row = 2)

            openxlsx::writeDataTable(
                wb = wb,
                sheet = sheetName,
                x = df,
                startCol = col.start,
                startRow = 3)

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
make_sets <- function(gene.lists){
    unique.genes <- unlist(gene.lists) %>% unique %>% as.character
    gene.df <- sapply(X = gene.lists, FUN = function(x) as.numeric(unique.genes %in% x)) %>%
        'rownames<-'(unique.genes) %>%
        as.data.frame

    set.names <- lapply(X = ncol(gene.df):1, FUN = function(len) {
        combn(x = colnames(gene.df), m = len, FUN = paste, collapse = '.')
    }) %>% unlist

    lapply(X = set.names, FUN = function(df.col.name){
            check.cols <- stringr::str_split(string = df.col.name, pattern = '\\.')[[1]]
            subset.rows <- subset(
                x = gene.df,
                subset = (rowSums(gene.df) == length(check.cols)) &
                    (rowSums(gene.df %>% dplyr::select(all_of(check.cols))) == length(check.cols))
            )

            info <- select(
                x = org.Hs.eg.db,
                keys = rownames(subset.rows),
                keytype = c('ENSEMBL'),
                columns = c('SYMBOL')
            )
    }) %>%
        setNames(set.names)
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
prepare_venn_plots <- function(dir, name){
    png.files <- list.files(path = dir, pattern = '\\.png$', full.names = TRUE)
    images <- lapply(X = png.files, FUN = image_read)

    # Change the title to proper nouns of abbreviations
    title <- stringr::str_split(string = name, pattern = '_', simplify = TRUE) %>%
        sapply(FUN = gsub, pattern = 'H', replacement = 'Hormone') %>%
        sapply(FUN = gsub, pattern = 'E', replacement = ' avec Enzalutamide', USE.NAMES = FALSE) %>%
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
    indir <- c('./dgea_output/de')

    outdir <- c('./dgea_output/venn_diagrams/Hormones_ARp')
    if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

    used.files <- list.files(path = indir, pattern = 'DMSO.csv', recursive = TRUE, full.name = TRUE) %>%
        '['(grep(pattern = 'DU145|ENZA', x = ., invert = TRUE))

    split.files <- split(
        x = used.files,
        f = stringr::str_extract(string = used.files, pattern = paste(hormone, collapse = '|'))
    )

    wb <- createWorkbook()

    for(group.name in names(split.files)){

        group <- split.files[[group.name]]
        tmp_outdir <- file.path(outdir, group.name)
        if(!dir.exists(tmp_outdir)) dir.create(tmp_outdir)

        loaded <- lapply(X = group, FUN = function(g) {
            readr::read_csv(g) %>% na.omit %>%
                dplyr::select(ensembl_gene) %>%
                as.matrix %>% as.character
        }) %>%
            setNames(stringr::str_extract(string = group, pattern = paste(cells, collapse = '|')))

        # Make Venn Diagrams
        venn.diagram(
            x = loaded,
            category.names = names(loaded),
            filename = file.path(tmp_outdir, paste0(group.name, '.png')))

        append_workbook_ensembl(wb = wb, gene.list = loaded, group.name = group.name)
    }
    openxlsx::saveWorkbook(wb = wb, file = file.path(outdir,'gene_list_PerHormone_ARp.xlsx'), overwrite = TRUE)
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

    venns_list <- list()
    my.widths <- c()

    for(element in venns){
        venns_list <- list.append(venns_list, element, NULL)
        my.widths <- c(my.widths, 0.75, 0.05)
    }

    venns_list[length(venns_list)] <- NULL
    my.widths <- my.widths[-length(my.widths)]

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
