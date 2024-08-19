library(tidyverse)
library(pheatmap)
library(ComplexHeatmap)
library(DESeq2)
library(ClassDiscovery)
library(org.Hs.eg.db)
library(grDevices)
library(RColorBrewer)
library(magrittr)


#' Convert Ensembl IDs to Gene Symbols
#'
#' This function converts Ensembl gene IDs to gene symbols for a given data frame.
#' The conversion is done using the `org.Hs.eg.db` package which provides
#' the necessary mappings.
#'
#' @param df A data frame with Ensembl gene IDs as row names. Each row should
#'   correspond to a gene, and the row names should be Ensembl IDs.
#'
#' @return A data frame with gene symbols as row names. The content of the data frame
#'   remains unchanged, but the row names are updated to gene symbols.
#'
#' @importFrom org.Hs.eg.db select
#' @importFrom AnnotationDbi keys
#' @importFrom AnnotationDbi columns
#' @import org.Hs.eg.db
#'
#' @examples
#' \dontrun{
#' # Example usage
#' # Assuming you have a data frame 'df' with Ensembl IDs as row names
#' df <- data.frame(matrix(rnorm(100), nrow=10))
#' rownames(df) <- c("ENSG00000141510", "ENSG00000284733", "ENSG00000238009",
#'                   "ENSG00000141837", "ENSG00000139618", "ENSG00000182185",
#'                   "ENSG00000138022", "ENSG00000182157", "ENSG00000163654",
#'                   "ENSG00000185102")
#' df_converted <- ensembl_to_symbol(df)
#' }
#' @export
ensembl_to_symbol <- function(df){

    # Get column of ensembl_id and remove it from df
    genes <- df$ensembl_gene
    df <- df %>% dplyr::select(-ensembl_gene)

    # Get associated symbols and assigned as row names
    symbols <- AnnotationDbi::select(org.Hs.eg.db, keys = genes, keytype = "ENSEMBL", columns = "SYMBOL")

    symbols <- symbols %>%
        group_by(ENSEMBL) %>%
        summarise(SYMBOL = paste(unique(SYMBOL), collapse = "/")) %>%
        ungroup()

    rownames(df) <- symbols$SYMBOL

    df
}

#' Create and Save a Heatmap with Clustering
#'
#' This function generates a heatmap with hierarchical clustering for both genes and samples.
#' The heatmap is saved to a file specified by the user. Clustering is performed using
#' Pearson correlation distance and average linkage.
#'
#' @param plot.data A numeric matrix where rows represent genes and columns represent samples.
#'   The values in the matrix should be suitable for clustering and heatmap visualization.
#'
#' @param filename A character string specifying the file path where the heatmap image will be saved.
#'   The file format is inferred from the file extension (e.g., ".png", ".pdf").
#'
#' @return The function does not return any value. It saves a heatmap image to the specified file.
#'
#' @importFrom pheatmap pheatmap
#' @importFrom stats hclust dist
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' # Example usage
#' # Create a sample data matrix
#' set.seed(123)
#' sample_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
#' rownames(sample_data) <- paste0("Gene", 1:10)
#' colnames(sample_data) <- paste0("Sample", 1:10)
#'
#' # Generate and save the heatmap
#' Make_HeatMap(plot.data = sample_data, filename = "heatmap.png")
#' }
#' @export
Make_HeatMap <- function(plot.data, filename, title, scale = TRUE, rcluster = TRUE,
                         ccluster = TRUE, color = NA, breaks = NA, dendogram = TRUE,
                         shorten = FALSE, display_numbers = FALSE, show_rownames = TRUE){

    # Set seed for reproductibility
    set.seed(seed = 42)

    if(rcluster){

        # Manually generate the clustering for genes
        row.cluster <- plot.data %>% t %>%
            distanceMatrix(metric="pearson") %>%
            hclust(method="average")

    } else { row.cluster <- FALSE }

    if(ccluster){

        # Manually generate the clustering for both samples
        col.cluster <- plot.data %>%
            distanceMatrix(metric="pearson") %>%
            hclust(method="average")

    } else{ col.cluster <- FALSE }

    if(scale){

        # Z - score
        means <- rowMeans(plot.data)
        var <- apply(X = plot.data, MARGIN = 1, FUN = var)
        plot.data <- (plot.data - means)/ sqrt(var)

    }

    if(any(is.na(color))){
        color <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100)
    }

    # Make and save de heatmap
    pheatmap::pheatmap(
        mat = plot.data,            # Data to be plotted

        # Cell esthetic
        border_color = NA,                     # Border
        cellwidth = ifelse(shorten, NA, 20),   # Width
        cellheight = ifelse(shorten, NA, 4),   # Height
        display_numbers = display_numbers,     # log2 Fold Change or Z-score
        color = color,                         # Color of ploted data

        # Legend
        breaks = breaks ,
        legends_breaks = c('upregulated', 'downregulated'),

        # Fontsizes
        fontsize = 5,
        fontsize_row = 4,           # Gene name
        fontsize_col = 4,           # Sample name
        fontsize_number = 3,        # log2 Fold Changes

        # Clustering
        cluster_rows = row.cluster, # Row (gene) clustering
        cluster_cols = col.cluster, # Column (sample) clustering
        treeheight_row = ifelse(dendogram, 50, 0) ,
        treeheight_col = ifelse(dendogram, 50, 0),
        lwd = 0.5,

        # Display names
        show_rownames = ifelse(shorten,FALSE,TRUE),       # Show gene names
        show_colnames = TRUE,       # Show sample names
        main = title,

        # Output
        filename = filename,        # Directory to save the heatmap

    )
}

#' Create and Save a Subset Heatmap Based on Column Pattern
#'
#' This function generates a heatmap from a subset of the columns in the input data.
#' The subset is defined by a specified pattern present in the column names. The
#' resulting heatmap is saved to a file.
#'
#' @param plot.data A data frame or matrix where rows represent genes (or features)
#'   and columns represent samples. The data should be suitable for clustering and
#'   visualization in a heatmap.
#' @param pattern A character string containing a regular expression pattern.
#'   Columns whose names match this pattern will be selected for the heatmap.
#'
#' @return This function does not return a value. It saves a heatmap image to a file.
#'   The filename is constructed using the specified pattern and some predefined
#'   variables (`file.output` and `name.lineage`).
#'
#' @importFrom pheatmap pheatmap
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' # Example usage
#' # Assuming you have a data frame 'plot.data' with appropriate data
#' subset_heatmap(plot.data, pattern = '11OHT')
#' }
#' @export
subset_heatmap <- function(plot.data, pattern, title, beg.filename, scale = TRUE,
                           unscale = FALSE, rcluster = TRUE, ccluster = TRUE){

    # Select column for sub HeatMap
    pattern_col_select <- grep(pattern = pattern, x = colnames(plot.data), value = TRUE)
    plot.data <- subset(plot.data, select = pattern_col_select)

    if(grepl(pattern = '11KT|11OHT', x = pattern)){

        infos <- correct.color_scaling(df=plot.data)
        color <- infos[[1]]
        breaks <- infos[[2]]

        } else {

        color <- NA
        breaks <- NA

        }

    if(scale){

        # Make scaled Sub HeatMap
        filename <- paste(beg.filename, pattern,'exclusif_scaled.pdf', sep = '_')
        Make_HeatMap(plot.data = plot.data, filename = filename, title = title,
                     rcluster = rcluster, ccluster= ccluster, color = color,
                     breaks = breaks)
    }

    if(unscale){

        # Make unscaled Sub HeatMap
        filename <- paste(beg.filename, pattern,'exclusif_unscaled.pdf', sep = '_')
        Make_HeatMap(plot.data = plot.data, filename = filename, title = title,
                     scale = FALSE, rcluster = rcluster, ccluster = ccluster,
                     color = color)
    }
}

#'
#' This function reads an Excel file and extracts specific columns based on the provided patterns.
#' It returns a unique list of genes after filtering out unwanted entries.
#'
#' @param excel.file A character string specifying the path to the Excel file.
#' @param sheet A character string or numeric value indicating the sheet name or index to read from the Excel file.
#' @param col.wanted A character vector specifying patterns to select desired columns. Columns matching these patterns will be included in the extraction.
#'
#' @return A character vector of unique gene names or identifiers after filtering out entries containing 'R1881' and any empty strings.
#'
#' @importFrom openxlsx read.xlsx
#' @importFrom dplyr %>%
#'
#' @examples
#' \dontrun{
#' # Example usage
#' # Assuming you have an Excel file with gene data
#' genes <- get_genes_excel_table("path/to/your/file.xlsx", sheet = 1, cols = c("11OHT", "Gene"))
#' }
#' @export
get_genes_excel_table <- function(xlsx.file, sheet, cols){

    map(xlsx.file, ~ {
        # Load file
        wb <- openxlsx::read.xlsx(.x, sheet = sheet)[-1,] %>%
            # Get columns related to 11oxy sets
            dplyr::select(all_of(cols))

        # Return a character vector
        unlist(x =  wb, use.names = FALSE) %>% na.omit %>% unique
    })
}

#' Subset and Plot Heatmap
#'
#' This function subsets the top 100 rows of the input data based on the values of a specified column
#' and generates a heatmap using the `subset_heatmap` function. The heatmap file is saved with a filename
#' that includes information about the ordering direction.
#'
#' @param plot.data A data frame containing the data to be plotted in the heatmap.
#' @param file.dir A character string specifying the directory where the heatmap file will be saved.
#' @param pattern A character string indicating the column name by which the data will be ordered and subset.
#' @param color A character string or vector specifying the color(s) to be used in the heatmap.
#' @param decreasing A logical value indicating whether the data should be sorted in decreasing order.
#'        Default is `FALSE` (ascending order).
#' @param scale A logical value indicating whether the data should be scaled before plotting. Default is `FALSE`.
#' @param unscale A logical value indicating whether the data should be unscaled after plotting. Default is `TRUE`.
#'
#' @details The function first orders the data based on the values in the specified column (`pattern`),
#'          selects the top 100 rows, and then generates a heatmap. The filename for the saved heatmap
#'          includes an indication of the ordering direction ('up' for ascending, 'down' for descending).
#'
#' @return This function does not return a value. It generates and saves a heatmap plot.
#'
#' @seealso `subset_heatmap`
#'
#' @examples
#' \dontrun{
#' # Example usage
#' data <- data.frame(A = rnorm(200), B = rnorm(200))
#' subset_heatmap_100(data, "output_directory", "A", color = "red", decreasing = TRUE)
#' }
#'
#' @export
subset_heatmap_100 <- function(plot.data, file.dir, pattern, color,
                               decreasing = FALSE, scale = FALSE, unscale = TRUE){

    # Subset data according to hormones
    plot.data <- plot.data[order(plot.data[[pattern]], decreasing = !decreasing),][1:100,]

    # Add info to filename
    direction <- ifelse(decreasing, 'down', 'up')
    filename <- paste(file.dir, paste0('100', direction), sep = '_')

    # Generate heatmap
    subset_heatmap(
        plot.data = plot.data,
        pattern = pattern,
        beg.filename = filename,
        scale = scale,
        unscale = unscale,
        rcluster = FALSE,
        ccluster = FALSE
    )
}

#' Correct Color Scaling for Heatmaps
#'
#' This function calculates the appropriate color scaling for a data frame containing both negative and positive values,
#' ensuring that the color gradient is centered around zero. It returns a list containing the color palette and corresponding breaks.
#'
#' @param df A numeric data frame. It is expected that the data frame contains both negative and positive values.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{color}: A character vector of color codes representing the gradient.
#'   \item \code{breaks}: A numeric vector representing the breaks corresponding to the color gradient.
#' }
#'
#' @details The function computes the proportion of negative values in the data to ensure that the color scale is balanced.
#'          It uses blue for negative values and red for positive values, with white in the middle representing zero.
#'          The function is useful for visualizing heatmaps where the zero value (neutral) should be visually distinguished.
#'
#' @examples
#' # Example usage
#' data <- matrix(rnorm(100, mean = 0, sd = 3), nrow = 10)
#' df <- as.data.frame(data)
#' color_info <- correct.color_scaling(df)
#' print(color_info)
#'
#' @importFrom grDevices colorRampPalette
#' @export
correct.color_scaling <- function(df){
    min.percent <- abs(min(df)) / (max(df) + abs(min(df)))

    # Number of colors for negative and positive values
    neg.colors <- round(100 * min.percent, digits = 0) # For negative values
    pos.colors <- 100 - neg.colors # For positive values

    # Define color palettes for negative and positive values
    neg_palette <- colorRampPalette(c('#2000e6', '#FFFFFF'))(neg.colors)
    pos_palette <- colorRampPalette(c('#FFFFFF', '#E60000'))(pos.colors+1)[-2]

    # Combine palettes, placing white in the middle for zero
    color <- c(neg_palette[-neg.colors], pos_palette[-1])

    # Calculate breaks for color mapping
    breaks <- c(seq(min(df), 0, length.out = neg.colors),
                seq(0, max(df), length.out = pos.colors + 1)[-1])

    list(color, breaks)
}



#' Generate Heatmap for AR+ Cell Lineages
#'
#' This function generates a heatmap for androgen receptor positive (AR+) cell lineages by processing and visualizing
#' gene expression data. It reads RDS files, extracts log2 fold changes, and matches gene symbols from an Excel file.
#'
#' @param pattern.files A character string specifying the pattern to match RDS files containing differential expression data.
#' @param rds_files A character string specifying the directory containing the RDS files.
#' @param dir.commun A character string specifying the directory containing the Excel file with common genes.
#' @param xlsx.sheetname A character vector specifying the sheet names in the Excel file to be read.
#' @param xlsx.col A character vector specifying the column names to include from the Excel sheets.
#' @param filename.heatmap A character string specifying the filename for the output heatmap image.
#' @param title.heatmap A character string specifying the title for the heatmap.
#'
#' @details The function performs the following steps:
#' \itemize{
#'   \item Reads RDS files containing differential expression data and extracts the log2 fold changes.
#'   \item Retrieves a list of genes from an Excel file.
#'   \item Merges the gene expression data based on common ensembl_gene identifiers.
#'   \item Filters and restructures the data to include only relevant genes.
#'   \item Converts ensembl_gene IDs to gene symbols.
#'   \item Generates and saves the heatmap with the specified filename and title.
#' }
#'
#' @return None. The function generates a heatmap and saves it to a file.
#'
#' @examples
#' # Example usage:
#' generate_whole_heatmap(
#'     pattern.files = "AR+",
#'     rds_files = "/path/to/rds/files",
#'     dir.commun = "/path/to/excel/files",
#'     xlsx.sheetname = c("Sheet1", "Sheet2"),
#'     xlsx.col = c("GeneName", "log2FoldChange"),
#'     filename.heatmap = "heatmap.png",
#'     title.heatmap = "AR+ Cell Lineages Heatmap"
#' )
#'
#' @importFrom dplyr select matches
#' @importFrom base readRDS
#' @export
generate_whole_heatmap <- function(
        pattern.files, rds_files, dir.commun,
        xlsx.sheetname, xlsx.col, filename.heatmap, title.heatmap,
        scale = TRUE
)
{
    # Generate heatmap for AR+ cell lineages (e.g. VCaP & LAPC4) ----
    r_object <-list.files(path = rds_files, pattern = pattern.files)

    sets <- sub(pattern='__DMSO_de.rds', replacement='', x=r_object)

    # Get exclusive gene's log2 Fold Change ----
    dir.r_objects <- file.path(rds_files, r_object)  # Get directory
    files.r_objects <- lapply(X = dir.r_objects, FUN = function(path){
        base::readRDS(file = path) %>%
            dplyr::select(ensembl_gene,log2FoldChange)
    }) %>% setNames(sets)

    # Get Excel file absolute path to directories ----
    file.commun <- list.files(path = dir.commun, pattern = '^[^~].*\\.xlsx$', full.names = TRUE)

    # Get list of genes exclusive to 11oxy ----
    genes.11oxy <- lapply(X = xlsx.sheetname, function(sheet){
        get_genes_excel_table(xlsx.file = file.commun, sheet = sheet , col.include = xlsx.col)
    }) %>% unlist %>% unique

    # Join gene expression into a single dataset on ensemble_id
    df.heatmap <- Reduce(function(x, y) merge(x, y, by = "ensembl_gene", all = TRUE), files.r_objects)
    colnames(df.heatmap) <- c('ensembl_gene', sets)

    # Restructure heatmap
    df.heatmap <- df.heatmap %>%
        dplyr::select(-matches("11OHT"), matches('11OHT')) %>%
        dplyr::select(-matches("ENZA"), matches('ENZA')) %>%
        subset(x=., subset = ensembl_gene %in% genes.11oxy) %>%
        na.omit

    # Switch ID from ensembl ID to symbol
    df.heatmap <- ensembl_to_symbol(df = df.heatmap)

    # Build output directory
    Make_HeatMap(plot.data = df.heatmap, filename = filename.heatmap,
                 ccluster = FALSE, dendogram = FALSE, title = title.heatmap,
                 shorten = FALSE, display_numbers = TRUE, scale = scale)

}
