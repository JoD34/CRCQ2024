### loading libraries ----
library(tidyverse)
library(fs)
library(msigdbr)
library(clusterProfiler)
library(openxlsx)
library(org.Hs.eg.db)
library(stats)
library(enrichplot)
library(ggridges)
library(rnaseq)
library(ggplotify)
library(pheatmap)
library(RColorBrewer)
library(pathview)
library(ReactomePA)

source('./scripts/R/util.R')

### Visualization plots ----
vis.dotplot <- function(obj, db, filename, path){

    # Make dotplot ----
    enrichplot::dotplot(
        object = obj,
        color = 'p.adjust',
        font.size = 10,
        title = paste(db, "Dotplot of enriched pathways", sep = ' - ')
    )

    # Save dotPlot
    ggsave(filename = filename, path = path )
}

vis.ridgeplot <- function(obj, db, filename, path){

    # Generate Ridge plot
    enrichplot::ridgeplot(
        x = obj,
        label_format = 30
    ) +
        # Reduce the character size for the axis informations
        theme(
            axis.text.x = element_text(size = 8),  # Adjust the x-axis text size
            axis.text.y = element_text(size = 8)   # Adjust the y-axis text size
        ) +
        labs(
            title = "Ridgeplot of enriched pathways",
            subtitle = paste('Database used : MSigDb', db, sep = ' - '),
            caption = Sys.time())

    # Save Ridge Plot
    ggsave(filename = filename, path = path )
}

vis.cnetplot <- function(obj, ranked.list, db, path, filename){
    # Useful to visualize the relationship between enriched pathways
    enrichplot::cnetplot(
        obj,
        categorySize = "p.adjust",
        showCategory = nrow(obj@result),
        node_label = "category",
        foldChange = ranked.list,
        colorEdge = FALSE,
        cex_label_category = 0.5
    ) +
        # Personalize the color scale for the Fold Change
        scale_color_gradient2(name='log2(FC)', low='darkblue', high='firebrick') +
        # Add labels to the general plot
        labs(
            title = 'Graph representation of shared expressed gene',
            subtitle = paste('Database : MSigDB -', db, 'set'),
            caption = Sys.time())

    # Save Concept Network Plot
    ggsave(filename = filename, path = path)
}

vis.gseaplot <- function(obj, path, txi.info){
    # Extract title information
    title <- obj@result$Description
    obj@result$Description <- ''
    # Extracted genes to subset the txi object for the heatmap
    enriched_genes <- obj@result$core_enrichment %>%
        stringr::str_split(., '/')
    # Generate the output directory
    out.gsea <- fs::path(path, 'plotGSEA')

    # Iterate over all enriched pathways
    for (i in seq_len(nrow(obj@result))) {
        # Generate GSEAplot
        gsea.plot <- enrichplot::gseaplot2(
            x = obj,
            geneSetID = i,
            color = "#86B875",
            title = paste(title[i],'(',txi.info$cell,'-',txi.info$hormone,')'),
            pvalue_table = TRUE,
            ES_geom = 'line'
        )

        # Susbet heatmap information
        heat.data <- txi$counts %>%
            as.data.frame %>%
            dplyr::select(matches(txi.info$cell)) %>%
            dplyr::select(matches(txi.info$hormone), matches('DMSO')) %>%
            dplyr::filter(rownames(.) %in% enriched_genes[[i]]) %>%
            setNames( sub('.*_(.*)_.*','\\1', colnames(.)))

        heat.plot <- pheatmap::pheatmap(
            heat.data,
            color = grDevices::colorRampPalette(
                rev(RColorBrewer::brewer.pal(n = 7, name ="RdBu")))(100),
            scale = 'row',
            cluster_cols = FALSE,
            treeheight_row = 0,
            cellwidth = 10,
            show_rownames = FALSE,
            na_col = 'darkgrey'
        )
        # Change obejct type from their respective type to ggplot
        gg.heat <- ggplotify::as.ggplot(heat.plot)
        gg.gsea <- ggplotify::as.ggplot(gsea.plot)

        # Concatenate plot unto a single pdf file
        gg.gsea + gg.heat + patchwork::plot_layout(widths = c(5,1))
        # Save plot
        ggsave(
            filename = paste(title[i] ,'gseaplot_heatmap.pdf', sep = '_'),
            path = out.gsea,
            height = 4.34, width = 8.31, units = 'in'
        )
    }
}

excel.table <- function(df){
    df$core_enrichment <- df$core_enrichment %>%
        map(~ {
            paste(switchIds(stringr::str_split_1(string = .x, pattern = '/')),
                  collapse = '/')
        })
    df
}

## Define paths ----
paths <- list(
    rds = './dgea_output/rds',
    output = './output_deseq2_pval0.01/gsea/clusterProfiler',
    kallisto = './deseq2_files/',
    anno = './files_util/anno'
)

paths$rds.output <- fs::path(paths$output, 'rds')
if(!dir.exists(paths$rds.output))dir.create(paths$rds.output, recursive = TRUE)

### Generate txi information ----
# Get kallisto files
k.files <- fs::dir_ls(paths$kallisto, regexp = '.h5$', recurse = TRUE) %>%
    '['(!grepl('DU145|ENZA',.))
samples <- fs::path_dir(k.files) %>% fs::path_file(.)
# Get annotation informations
anno.file <- fs::dir_ls(paths$anno, regexp = 'protein_coding.csv$')

# Import full kallisto informations
txi <-  rnaseq::import_kallisto(
    filenames = k.files,
    anno = anno.file,
    ignoreTxVersion = TRUE
)
colnames(txi$counts) <-  samples
colnames(txi$abundance) <- samples
colnames(txi$length) <- samples

# Set Long name for the abbreviations of the gene sets' name
categories <- c(
    H = 'Hallmark Gene Set',
    C1 = 'Positional Gene Set',
    C2 = 'Curated Gene Set',
    C3 = 'Regulatory targe Gene Set',
    C4 = 'Computational Gene Set',
    C5 = 'Ontology Gene Set',
    C6 = 'Oncogenic Signature Gene Set',
    C7 = 'Immunologic Signature Gene Set',
    C8 = 'Cell Type Signature Gene Set'
)

# Generate named vectors
deg.files <- fs::dir_ls(paths$rds, regexp = '.*/(VCaP|LAPC4)__.*[^ENZA]__DMSO_de.rds$') %>%
    purrr::map(~ {
        # Read DESeq2 output
        readRDS(.x) %>%
            # Extract specific columns from DESeq2 output
            dplyr::select('ensembl_gene','log2FoldChange', 'padj') %>%
            # Keep only significant sets
            dplyr::filter(padj < 0.05) %>%
            # Order base on log2FC
            '['(order(.$log2FoldChange, decreasing = TRUE), ) %>%
            # Get rid of NAs
            na.omit() %>%
            # Generate the ranked named vectors
            { setNames(.$log2FoldChange, .$ensembl_gene) }
    })

# GSEA ----
iwalk(deg.files, ~ {

    # Extract ranked gene list
    geneList <- .x

    ## Extract metadata from path ----
    # Extract cell lineage information
    f <- fs::path_file(.y)
    rds.f <- sub('(.*)__.*','\\1', f)
    dir.out <- sub('^(.*?)__(.*?)__.*$', '\\1\\\\\\2',  f)

    # Generate output pathways
    cell.out <- fs::path(paths$output, dir.out)
    if(!dir.exists(cell.out)) dir.create(cell.out, recursive = TRUE)

    # Create a new workbook to stock enriched pathways' information
    wb <- openxlsx::createWorkbook()

    # Iterate over gene set from MSigDb ----
    purrr::walk(names(categories), ~ {

        # Download MSigDb gene sets for a given categories
        term2gene <- msigdbr::msigdbr(species = 'Homo sapiens', category = .x) %>%
            dplyr::select(gs_name, ensembl_gene)

        # Make GSEAnalysis against the gene set of all MSigSB Categories
        res <- clusterProfiler::GSEA(
            geneList = geneList,
            TERM2GENE = term2gene,
            eps = 1e-300,
            pvalueCutoff = 0.05,
            seed = 42
        )

        # Save GSEAresult object for comparison analysis
        res.path <- fs::path(paths$rds.output, paste0(paste(rds.f, .x, sep = '_')), ext = 'rds')
        saveRDS(object = res, file = res.path)

        # Write enriched information to an .xlsx file
        openxlsx::addWorksheet(wb, sheetName = .x)

        # Reformat output table
        table <- excel.table(df = res@result)
        openxlsx::writeDataTable(wb, sheet = .x, x = table, xy = c(2, 2))

        # Adjust columns width for the data table
        openxlsx::setColWidths(wb, sheet = .x, cols = 2:(ncol(res@result) + 1), widths = "auto")

        if((.x == 'H') & (nrow(table) > 0)){
            cat.title <- categories[.x]
            tmp.dir <- fs::path(cell.out, .x)
            if(!dir.exists(tmp.dir)) dir.create(tmp.dir)
            res@result$Description <- str_remove(res@result$ID, 'HALLMARK_') %>%
                gsub('_', ' ', .)

            # Make dotplot ----
            vis.dotplot(
                obj = res,
                db = cat.title,
                path = tmp.dir,
                filename = paste0(.x, '_dotplot.pdf'))

            # Make RidgePlot ----
            vis.ridgeplot(
                obj = res,
                db = cat.title,
                filename = paste0(.x, '_ridgeplot.pdf'),
                path = tmp.dir )

            # Make concept network plot ----
            vis.cnetplot(
                obj = res,
                ranked.list = geneList,
                db = cat.title,
                path = tmp.dir,
                filename = paste0(.x,'_cnetplot.pdf'))

            # Plot GSEA plot with heatmap informations -----
            meta <- sub('(.*)_.*','\\1', fs::path_file(res.path)) %>%
                str_split_1(.,'__')

            vis.gseaplot(
                obj = res,
                path = tmp.dir,
                txi.info = list(cell = meta[1], hormone = meta[2])
            )

        }
    })

    # Save enriched pathways informations from GSEA in an .xlsx file
    openxlsx::saveWorkbook(wb, file = paste0(cell.out,'_MSigDb.xlsx'), overwrite = TRUE)
})


# TODO -> Implement the other source of GSEA

#entrezID.list <- setNames(geneList, switchIds(ids = names(geneList), to = 'ENTREZID'))

## Kegg GSEA
#kegg.res <- clusterProfiler::gseKEGG(geneList = entrezID.list,eps = 1e-300)

## REACTOME GSEA
#reactome.res <- ReactomePA::gsePathway(geneList = entrezID.list, eps = 1e-300)

### KEGG ----
kegg.analysis <- function(){
    pathview(
        gene.data = geneList,
        organism = "hsa",
        gene.idtype = "ENSEMBL",
        pathway.id = 'hsa05215',
        species = 'hsa',
        high = list(gene = "firebrick", cpd = "yellow"),
        low = list(gene = "darkgreen", cpd = "blue")
    )

}

### GO ----
GO.GSEA <- function(){
    # Do GSEA on GO terms
    gseGO_res <- clusterProfiler::gseGO(
        geneList = geneList,
        ont = 'BP',
        keyType = 'ENSEMBL',
        OrgDb = 'org.Hs.eg.db',
        eps = 1e-300,
        pAdjustMethod = "BH",
        verbose = FALSE,
        seed = 42
    )

    title <- gseGO_res@result$Description
    gseGO_res@result$Description <- ''

    p <- gseaplot2(
        x = gseGO_res,
        geneSetID = 1,
        title = title[1],
        color = c("#E495A5", "#86B875", "#7DB0DD"),
        pvalue_table = TRUE,
        ES_geom = 'line'
    )

    ggsave(filename = filename, path = paths$output)
}
