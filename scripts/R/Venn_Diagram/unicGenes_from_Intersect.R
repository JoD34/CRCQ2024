library(tidyverse)
library(openxlsx)
library(magrittr)


intersect2unic <- function(input.dir, output.dir, sheets, cols, intersect = FALSE){

    # Check for validity of inputed data
    if(!file.exists(input.dir)) stop("The input file doesn't exist")
    if(!all(sheets %in% openxlsx::getSheetNames(input.dir))) {
        stop("Not all specified sheet aren't in the .xlsx file")
    }

    # Include intersect if asked for in the parameters
    if(intersect) cols <- c(cols, paste(cols, collapse = '.'))

    df <- lapply(X = sheets, FUN = function(s){

        # Load excel workbook
        doc <- openxlsx::read.xlsx(xlsxFile = input.dir, sheet=s)

        # Check if all wanted column exist
        if(!all(cols %in% colnames(doc))) stop("Some columns you're looking for don't exist")

        # Select subset of columns (correspond to 11oxy) and remove 'ENSEMBL'
        subset(x = doc, select = cols) %>%

        # Flatten the dataframe to get a character vector
        unlist(x=., use.names=FALSE)

    }) %>%

        # Explode list to get a character vector
        unlist(x=., use.names=FALSE) %>%

        # Extract unique ENSEMBL gene from AR+ lineage
        unique %>%

        # Remove all that isn't an ENSEMBL IDs
        subset(x = ., subset = grepl(pattern='^ENSG.*', x=.)) %>%

        # Place the genes in ordre of number
        '['(order(.))

    # Create output directory if doesn't exist already
    if(!dir.exists(output.dir)) dir.create(output.dir)

    # Write ENSEMBLE IDs to file
    readr::write_csv(
        x = data.frame(ENSEMBL_ID = df),
        file = file.path(output.dir, 'gene_list.csv')
    )
}
