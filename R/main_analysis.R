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
dgea <- function(comparison, kallisto_path=getwd(), util, outfile_deg){
    namefile <- paste(comparison, collapse='_vs_')
    regX <- sapply(comparison, function(x) paste(x, 'Rep?', sep='_')) %>%
        paste(collapse='|')

    myDir <- list.files(path=path_kallistos, pattern=regX, full.names=TRUE)
    file <- file.path(myDir, "abundance.h5")

    # Generate the anntotation file
    out_dir <- make_annotation_file(
        name = namefile,
        path = myDir,
        path_out = util
        )

    # Import Kallisto data
    txi <- import_kallisto(filename =file, anno=out_dir)

    mySamples <- get_samples_name(bundle=myDir)
    txi <- modify_colnames(df=txi, names=mySamples)

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
    make_volcano(
        deg,
        fc.seuil=1.5,
        p.seuil=0.01,
        title=comparison,
        path='./dgea_output/volcanos',
        filename = namefile
        )

    pdf(file=file.path(('./dgea_output/data_visualization/'),
                       paste0(namefile, '.pdf')))
    heat <- make_expr_heatmap(dds)
    dist.mat <- make_dist_matrix(dds=dds)
    pca <- make_PCA(
        dds=dds,
        intgroup=condition,
        titre=unlist(strsplit(namefile, split='_'))
        )
    print(pca)
    dev.off()

    write_de_files(
        de=deg,
        txi=txi$anno,
        path=outfile_deg
        )
}
