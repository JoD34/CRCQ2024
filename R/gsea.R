library(fgsea)
library(data.table)
library(ggplot2)
library(gprofiler2)
library(GSVA)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)

pathways <- msigdbr::msigdbr(species = 'Homo sapiens')



fgseaRes <- fgsea(pathways = pathways,
                  stats    = exampleRanks,
                  minSize  = 15,
                  maxSize  = 500)

head(fgseaRes[order(pval), ])
