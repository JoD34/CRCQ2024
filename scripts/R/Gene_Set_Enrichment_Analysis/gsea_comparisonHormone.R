library(clusterProfiler)
library(tidyverse)
library(fs)
library(msigdbr)
library(org.Hs.eg.db)

source('./scripts/R/util.R')

### Set variables
paths = list(
    rds.input = './dgea_output/gsea/clusterProfiler/MSigDb_padj0.05/rds/',
    de.rds.input = './dgea_output/rds/'
)

### Read rds output files of clusterProfiler::GSEA
f <- fs::dir_ls(path = paths$rds.input, regexp = '.*_H.rds$') %>%
    split(., sub('(VCaP|LAPC4).*','\\1', fs::path_file(.)))

f2 <- fs::dir_ls(
    path = paths$de.rds.input,
    recurse = TRUE,
    regexp = '.*(VCaP|LAPC4)__.*([^ENZA])__DMSO_de.rds$'
)

g.cell <- sub('(.*)__.*__.*','\\1', fs::path_file(f2))
g.hormone <- sub('.*__(.*)__.*', '\\1', fs::path_file(f2))

deg <- map(f2, ~ readRDS(.x) %>%
               dplyr::select(all_of(c('id', 'log2FoldChange', 'padj'))) %>%
               dplyr::filter(padj < 0.05) %>%
               dplyr::mutate(
                   cell = sub('(.*)__.*__.*','\\1', fs::path_file(.x)),
                   hormone = sub('.*__(.*)__.*', '\\1', fs::path_file(.x))
                   ) %>%
               dplyr::arrange(desc(log2FoldChange)) %>%
               {setNames(.$log2FoldChange, .$id)}
            ) %>%
    setNames(sub('.*/(.*)__DMSO_de.rds','\\1',f2))

# Load MSigDb Hallmark gene set
H <- msigdbr::msigdbr(category = 'H') %>%
    dplyr::select(gs_name, ensembl_gene) %>%
    dplyr::mutate_all(~ sub('^HALLMARK_', '', .))

res <- compareCluster(
    geneClusters = deg,
    fun = 'GSEA',
    TERM2GENE = H,
    eps = 1e-300
)

# Add metadata to the output table
res@compareClusterResult$Hormone <- sub('.*__(.*)','\\1', res@compareClusterResult$Cluster)
res@compareClusterResult$cell <- sub('(.*)__.*','\\1', res@compareClusterResult$Cluster)
res@compareClusterResult$Cluster <- sub('__','.', res@compareClusterResult$Cluster)
res@compareClusterResult$Description <- gsub('_',' ', res@compareClusterResult$Description)

# Graph data using various presentation plot
dotplot(res, x = 'Hormone') +
    facet_grid(~ cell)

enrichplot::cnetplot(
    res,
    node_label = 'category',
    circular = FALSE,
    layout = 'kk',
    cex.params = list(category_node = 1)
    )

plots <- map2(f, meta, function(obj, txi.info){

    txi.info <- list(cell = txi.info[1], hormone = txi.info[2])

    # Extract title information
    title <- 'ANDROGEN RESPONSE'
    obj@result$Description <- ''

    # Extracted genes to subset the txi object for the heatmap
    enriched_genes <- obj@result$core_enrichment[[1]] %>% stringr::str_split(., '/') %>% '[['(1)

    # Generate GSEAplot
    gsea.plot <- enrichplot::gseaplot2(
        x = obj,
        geneSetID = grep('ANDROGEN', obj@result$ID),
        color = "#86B875",
        title = paste(title,'(',txi.info$cell,'-',txi.info$hormone,')'),
        pvalue_table = TRUE,
        ES_geom = 'line'
    )

    # Susbet heatmap information
    heat.data <- txi$counts %>%
        as.data.frame %>%
        dplyr::select(matches(txi.info$cell)) %>%
        dplyr::select(matches(txi.info$hormone), matches('DMSO')) %>%
        dplyr::filter(rownames(.) %in% enriched_genes) %>%
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
    p <- gg.gsea + gg.heat + patchwork::plot_layout(widths = c(5,1))
    p
})

p <- cowplot::plot_grid(plotlist = plots, align = 'hv', ncol = 1, labels = LETTERS[1:2])
ggsave(filename = 'Rplot01.pdf', path = './dgea_output/gsea/clusterProfiler/MSigDb_padj0.05/',
       width = 8.5, height = 8, unit = 'in')
