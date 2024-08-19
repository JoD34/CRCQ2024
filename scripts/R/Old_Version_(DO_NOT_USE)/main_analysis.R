## Load libraries ----
library(rnaseq)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

# dgea ----
###
### dgea(dds)
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
dgea <- function(comparison, kallisto_path=getwd(), util){

    # relative paths
    path_deg <- './dgea_output/dds_output'
    path_deg_list <- './dgea_output/deg_list'
    path_txi <- 'dgea_output/count_tmp'

    namefile <- paste(comparison, collapse='_vs_')
    regX <- sapply(comparison, function(x) paste(x, 'Rep?', sep='_')) %>%
        paste(collapse='|')

    myDir <- list.files(path=kallisto_path, pattern=regX, full.names=TRUE)
    file <- file.path(myDir, "abundance.h5")

    # Generate the anntotation file
    out_dir <- make_annotation_file(
        name = namefile,
        paths = myDir,
        path_out = util
        )

    # Import Kallisto data
    txi <- import_kallisto(filename =file, anno=out_dir)

    mySamples <- get_samples_name(bundle=myDir)
    txi <- modify_colnames(df=txi, names=mySamples)

    # Write data generated to file
    write_txi_file(
        txi=txi,
        comparison=comparison,
        path=path_txi,
        namefile=namefile
        )

    # Make PCA
    pca_df <- produce_pca_df(txi=txi)
    plot_pca(res_pca=pca_df)

    # Make formula used in DESeq2 analysis
    condition <- ifelse(
        test=sum(grepl(pattern='ENZA', x=comparison)) == 1,
        yes='enza',
        no='hormone'
        )
    f <- formula(paste("~",condition))

    # Make design table
    design <- make_design(
        samples=mySamples,
        cnames=c('lignee', 'hormone', 'enza', 'replicate'),
        condition=condition
    )

    # Differencial analysis using DESeq2
    dds <- deseq2_analysis(txi=txi, design=design, formula=f)

    # Get results of differencial analysis
    deg <- DESeq2::resultsNames(dds) %>%
        '['(!grepl(pattern="Intercept",x=.)) %>%
        DESeq2::results(object=dds, name=.)

    # Generate graphs
    fc.seuil=0.01
    p.seuil=0.05
    make_volcano(
        deg,
        fc.seuil=fc.seuil,
        p.seuil=p.seuil,
        title=comparison,
        path='./dgea_output/volcanos',
        filename = namefile
        )
    plots <- list(
        make_expr_heatmap(dds),
        make_dist_matrix(dds)
        )
    pca <- make_PCA(
        dds=dds,
        intgroup=condition,
        titre=unlist(strsplit(namefile, split='_'))
    )
    pdf(file=file.path(('./dgea_output/data_visualization/'),
                       paste0(namefile, '.pdf')))
    invisible(lapply(plots, function(p){
        grid::grid.newpage()
        grid::grid.draw(p$gtable)
        }))

    print(pca)
    dev.off()
    dev.off()

    # Get significant genes (p_adj < 0.05)
    seuils <- deg$padj < p.seuil
    deg <- deg[!is.na(deg$padj) & seuils, ]

    # Write DEG files
    write_deg(de=deg, anno=txi$anno, path=path_deg, namefile=namefile)
    readr::write_csv(
        x=data.frame(rownames(deg)),
        file=file.path(path_deg_list, paste0(namefile, '.csv')),
        col_names = FALSE
    )
}
