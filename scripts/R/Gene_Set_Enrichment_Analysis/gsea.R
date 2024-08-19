### Clean environement ----
rm(list = ls(all.names = TRUE))
gc()
options(max.print = .Machine$integer.max, stringsAsFactors = FALSE, dplyr.summarise.inform = FALSE)

### Libraries ----
# Analyse d'enrichissement
library(gprofiler2) # Tools for accessing the GO enrichment results using g:Profiler webtool
library(GSVA) # Gene Set Variation Analysis
library(clusterProfiler) # Provides a suite of tools for functional expression analysis; Make sure R-Forge isn't a repo when installing
library(enrichplot) # Great for making the standard GSEA enrichment plote;  Make sure R-Forge isn't a repo when installing
# GSEA
library(GSEABase) # Functions and methods for Gene Set Enrichment analysis
library(Biobase) # Base function for bioconductor; required by GSEAbase
library(fgsea) # Fast Gene Set Enrichment Analysis
library(BiocParallel) # Reimplementation of certain functions
library(DOSE)
# Biological annotations
library(org.Hs.eg.db)
library(AnnotationDbi)
library(biomaRt)
# databases
library(msigdbr) # Molecular Signature Database (for R)
library(ReactomePA) # Reactome Pathways
library(rWikiPathways)
library(meshes)
library(MeSHDbi)
# Manipulation de données
library(DT) # Interactive ans searchable tables of our GSEA results
library(htmltools)
library(stats)
library(openxlsx) # Set of tool to write an .xlsx (Excel) file
library(tidyverse)
# graphical package
library(ggplot2) # Graph it baby
library(ggplotify)
library(ggdag) # Used for GO graph (eg nodes and edges)
# Custom code
source('./scripts/R/util.R')
source('./scripts/R/Gene_Set_Enrichment_Analysis/gsea.prep_rank.functions.R')
source('./scripts/R/Pathway_Enrichment_Analysis/enricher.R')

## Function: Main ----

# Manage directories
intercept.files <- './dgea_output/venn_diagrams/FC_NA_pAdj_0.05/H_vs_DMSO/H_vs_DMSO_GeneLists.xlsx'
rds.files <- './dgea_output/rds'
path.gsea.output <- './dgea_output/gsea'
if(!dir.exists(path.gsea.output)) dir.create(path.gsea.output)

cells <- c('LAPC4', 'VCaP')

# Prepare ranking of gene expression - Custom function
ranks <- prepare_ranking(
    rds.files = rds.files,
    intercept.files = intercept.files,
    cells = cells
)

# Prepare gene set to check againt for significancy - Custom function
gs <- prepare.gene_sets()

# Gene Set Enrichment Analysis - Fast GSEA ----
# Using fgseaMultilevel is more accurate than the simple fgsea function
# Used mainly for small gene sets or small dataset (fewer than a coupl thousand genes)
fgseaRes <- fgsea::fgseaMultilevel(
    pathways = gs$msigdb, # Gene Set used
    stats = ranks,        # Ranked genes
    minSize = 5,          # Minimal size of a gene set to test. Reduce noise.
    nPermSimple = 10000    # Number of permutation for sampling to do
)

sig.set <- fgseaRes[fgseaRes$pval <= 0.05] %>%
    dplyr::arrange(pval)

# Visulise enriched pathways with interactive table
DT::datatable(data = sig.set, options = list(pageLength = 10))

## Switch ENSEMBL GENE IDs to SYMBOLs
# Get unique ENSEMBL IDs
ensembl <- sig.set$leadingEdge %>% unlist %>% unique

# Get corresponding SYMBOL
symbols <- switch_ids(input = ensembl, in.type = 'ENSEMBL', out.type = 'SYMBOL')

# Generate a named character vector to switch faster
ensembl_to_symbols <- setNames(symbols$SYMBOL, symbols$ENSEMBL)

# Remplace the listed ENSEMBL IDS with SYMBOLS
sig.set$leadingEdge <- lapply(sig.set$leadingEdge, FUN = function(s) ensembl_to_symbols[s])

## Write dataframe to .xlsx files
write2xlsx(df = sig.set, sheetNames = 'GSEA - MSigbDb', out.dir = path.gsea.output, file.name = 'fgseaMultilevel.xlsx')

## Go Enrichment ----
go.res <- gprofiler2::gost(query = names(ranks), organism = 'hsapiens', correction_method = 'fdr')
go.plot <- gprofiler2::gostplot(gostres = go.res, interactive = TRUE, capped = TRUE)

# Save GO enrichment plot
gprofiler2::publish_gostplot(
    p = go.plot,
    filename = file.path(path.gsea.output, 'gprofiler2_GOenrich_plot.pdf')
)

# Save GO enrichment table
write.go.enrich <- go.res$result %>%
    dplyr::select(c('source', 'term_id','term_name',  'term_size', 'intersection_size', 'p_value')) %>%
    dplyr::filter(term_size > 5)

write2xlsx(df = write.go.enrich, sheetNames = 'GOenrichment', out.dir = path.gsea.output, file.name = 'gprofiler2_GOenrich_plot.xlsx')

prepare.gene_sets <- function(){

    ### Categories of MSigDb explained ----
    #
    #   Link to details: https://www.gsea-msigdb.org/gsea/msigdb
    #
    #   'H' - Hallmark gene sets : coherently expressed signatures derived by
    #               aggregating many MSigDB gene sets to represent well-defined
    #               biological states or processes.
    #
    #   'C1' - Positional gene sets : human chromosome cytogenetic bands.
    #
    #   'C2' - Curated gene sets : online pathway db, publications in PubMed,
    #               and knowledge of domain experts.
    #
    #   'C3' - Regulatory target gene sets : gene target predictions for microRNA
    #               seed sequences and predicted transcription factor binding sites.
    #
    #   'C4' - Computational gene sets : defined by mining large collections of
    #               cancer-oriented expression data.
    #
    #   'C5' - Ontology gene sets : genes annotated by the same ontology term.
    #
    #   'C6' - Oncogenic signature gene sets : defined directly from microarray
    #               gene expression data from cancer gene perturbations.
    #
    #   'C7' - Immunologic signature gene sets : cell states and perturbations
    #               within the immune system.
    #
    #   'C8' - Cell type signature gene sets : curated from cluster markers identified
    #               in single-cell sequencing studies of human tissue.
    # ----

    ### Generating gene sets with MSigDb ----
    categories <- c('H', 'C2', 'C3', 'C4', 'C6')

    # Load the gene sets from MSigDb
    gs.msigdb <- msigdbr::msigdbr(species = 'Homo sapiens') %>%
        # Subset the database to the selected categories
        dplyr::filter(gs_cat %in% categories) %>%
        # Split the dataframe by the values of the gs_name (group by gene sets)
        dplyr::group_by(gs_name) %>%
        # list all ensembl_gene ids from a given gene sets
        dplyr::summarise(ensembl_gene = list(ensembl_gene), .groups = 'drop') %>%
        # Restrict the number of gene sets for those having more than 5 genes
        dplyr::filter(lengths(ensembl_gene) > 5)

    # Set list name (Assign gene set name to its corresponding gene)
    gene.sets.msigdb <- setNames(gs.msigdb$ensembl_gene, gs.msigdb$gs_name)

    ### Generating gene sets with rWikiPathways
    path.db <- "./files_util/databases"
    if(!dir.exists(path.db)) dir.create(path = path.db, recursive = TRUE)
    path.file.db <- file.path(path.db, list.files(path.db, '.gmt$'))

    # Download WikiPathway archive if the file isn't already downloaded
    if(!file.exists(path.file.db)){
        rWikiPathways::downloadPathwayArchive(
            organism = 'Homo sapiens',
            format = 'gmt',
            destpath = path.db
        )
    }

    # Load db File
    w.pathway <- rWikiPathways::readPathwayGMT(file = path.file.db)%>%
        dplyr::group_by(name) %>%
        dplyr::summarise(entrezid = list(gene), .groups = 'drop') %>%
        dplyr::filter(lengths(entrezid) > 5)

    gene.sets.wikipathway <- setNames(w.pathway$entrezid, w.pathway$name)


    # Return the gene sets
    list(
        msigdb = gene.sets.msigdb,
        wikipathway = gene.sets.wikipathway
    )
}

# Gene Set Enrichment Analysis - Cluster Profiler ----
GSEA.CP <- function(ranks){

    # gseKEGG - GSEA using KEGG
    gseKEGG.res <- clusterProfiler::gseKEGG(geneList = ranks, seed = 3)

    # gseGO - GSEA using GO term
    lapply(X = c('BP', 'MF', 'CC'), FUN = function(ont){

        # Make Gene Set Enrichment with GO terms
        gseGO.res <- clusterProfiler::gseGO(geneList = ranks, ont = ont, OrgDb = org.Hs.eg.db, seed = 1)

        # If pathways are found to be enriched
        if(nrow(gseGO.res@result != 0)){
            gseGO.plot <- clusterProfiler::plotGOgraph(x = gseGO.res)
            plot(gseGO.plot$complete.dag)
            }
        })

    # gseMKEGG - GSEA using MKEGG
    gseMKEGG.res <- clusterProfiler::gseMKEGG(geneList = ranks, seed = 2)

    # GSEA using MSigDB - Molecular Signature Databases
    files <- list.files('./dgea_output/gene_sets/modified/', pattern = '.*\\.rds$', full.names = TRUE)
    gs <- split(
        x = lapply(X = files, FUN = readRDS),
        f = sub(pattern = '.*/(c[0-9]+|h).*', replacement = '\\1', x = files)
    )

    purrr::imap(ranks, ~ {
        wb <- openxlsx::createWorkbook()

        purrr::iwalk(gs, ~ {
            res <- purrr::map_dfr(.x, ~ clusterProfiler::GSEA(geneList = ranks, TERM2GENE = .x)@result)

            # Add a new WorkBook sheet & Write data
            openxlsx::addWorksheet(wb, sheetName = .y)
            openxlsx::writeDataTable(wb, sheet = .y, x = res, xy = c(2, 2))

            # Adjust columns width for the data table
            openxlsx::setColWidths(wb, sheet = .y, cols = 2:(ncol(res) + 1), widths = "auto")

        })

        # Save enriched data in a .xlsx file
        filename <- paste0(.y, paste('regulated', 'MSigDb','enrichR.xlsx', sep = '_' ))
        openxlsx::saveWorkbook(wb, file = file.path(path.gsea.output, filename), overwrite = TRUE)

    })

}

# Gene Set Enrichment Analysis - Reactome ----
GSEA.REACTOME <- function(ranks){
    ReactomePA::gsePathway(geneList = ranks, seed = 4)
}

# Gene Set Enrichment Analysis - Medical Subjects Headings ----
GSEA.MESH <- function(){
    meshes::gseMeSH(geneList = ranks, MeSHDb = )
}

# Gene Set Enrichment Analysis - Disease Ontology Semantic Analysis ----
GSEA.DOSE <- function(){

}
