# Load libraries ----
library(RColorBrewer)
library(DESeq2)
library(tidyverse)

# make_annotation_file ----
###
### make_dist_matrix(dds)
###
##  Summary
##      Visualize the distance between samples
##
##  Args
##      dds: (DESeqDataSet) Object generated from the DESeq2 package
##
##
##  Examples
##      txi <- get_demo_txi()
##      design <- get_demo_design()
##      dds <- deseq2_analysis(txi, design, ~ group)
##      make_dist_matrix(dds)
##
make_dist_matrix <- function(dds){

    # Data stabilization for distance computation
    vsd <- DESeq2::vst(dds, blind=F)

    # Compute sample-to-sample distance
    sampleDists <- dist(t(assay(vsd)))
    sampleDistMatrix <- as.matrix(sampleDists)

    # Generate heatmap of the distance matrix
    rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$sample, sep=":")
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette(rev(brewer.pal(name="BuPu", n=9)))(255)
    pheatmap(
        sampleDistMatrix,
        clustering_distance_rows=sampleDists,
        clustering_distance_cols=sampleDists,
        col=colors,
        main="Distance matrix")
}

# make_volcano ----
###
### make_volcano(dds)
###
##  Summary
##      Represent the fold change of gene expression in accordance with
##          the significancy of the change
##
##  Args
##      de: (DESeqResults) Object generated from the DESeq2 package,
##          result from the comparison
##
##
##  Examples
##      txi <- get_demo_txi()
##      design <- get_demo_design()
##      dds <- deseq2_analysis(txi, design, ~ group)
##      make_expr_heatmap(dds)
##
make_volcano <- function(de, fc.seuil, p.seuil, title, path, filename){

    # Generate information for further usage in volcano plot
    titre <- paste('Expression différentielle:',
                   title[1], '(Défaut)', 'vs', title[2], sep=' ')
    subtitle <- paste('seuil pValue:', p.seuil, ';',
                      'seuil fold-change:', fc.seuil, sep=' ')
    caption <- paste("Produced on", Sys.time())

    # Generate and save volcano plot
    vol <- produce_volcano(
        de_res=de,
        fc_threshold=fc.seuil,
        p_threshold=p.seuil,
        size=1,
        title=titre,
        subtitle=subtitle,
        caption=caption
    )
    ggsave(filename=file.path(path, paste0(filename, '.pdf')))
}

# make_expr_heatmap ----
###
### make_expr_heatmap(dds)
###
##  Summary
##      Visualize the distance between samples
##
##  Args
##      dds: (DESeqDataSet) Object generated from the DESeq2 package
##
##
##  Examples
##      txi <- get_demo_txi()
##      design <- get_demo_design()
##      dds <- deseq2_analysis(txi, design, ~ group)
##      make_expr_heatmap(dds)
##
make_expr_heatmap <- function(dds){

    # Data stabilization for further computation
    vsd <- vst(dds, blind=F)

    # Get counts of dds object (e.g. DESeq2 analysis)
    mine <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:200]
    mine <- assay(vsd)[mine,]
    # Format data
    df <- as.data.frame(colData(dds))
    anno_heat <- df %>%
        dplyr::select('lignee', 'hormone', 'replicate')

    # Generate heatmap
    pheatmap(
        mat=mine,
        show_colnames=FALSE,
        show_rownames=FALSE,
        annotation_col=anno_heat
        )

}

# make_PCA ----
###
### make_PCA(dds)
###
##  Summary
##      Visualize principal component analysis
##
##  Args
##      dds: (DESeqDataSet) Object generated from the DESeq2 package
##
##
##  Examples
##      txi <- get_demo_txi()
##      design <- get_demo_design()
##      dds <- deseq2_analysis(txi, design, ~ group)
##      make_dist_matrix(dds)
##
make_PCA <- function(dds, intgroup, titre){

    # Data stabilization for further computation
    vsd <- vst(dds, blind=F)

    # PCA analysis
    pcaData <- plotPCA(vsd, intgroup=intgroup, returnData=T, ntop=50000)
    pcaData$rep <- stringr::str_extract(string=pcaData$name, pattern='Rep[0-9]')
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    namesPCA <- attr(pcaData,"names")
    graph_title <- paste("Principal Component Analysis (PCA):",
                         paste(titre, collapse=' '),
                         sep=' ')

    # plotPCA
    ggplot(pcaData) +
        aes(PC1, PC2, color=group, shape=rep) +
        geom_point(size=3) +
        labs(title = graph_title,
             x = paste0(namesPCA[1], ": ", percentVar[1], " %"),
             y = paste0(namesPCA[2], ": ", percentVar[2], " %"),
             caption = paste("Produced on", Sys.time())) +
        coord_fixed() +
        theme_bw()
}
