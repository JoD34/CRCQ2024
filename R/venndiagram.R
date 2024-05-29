library(tidyverse)
library(VennDiagram)


make_names <- function(){
    # Loaded gene files
    cell <- c('DU145', 'VCaP', 'LAPC4')
    categories <- c('R1881', '11KT', '11OHT')
    base <- 'DMSO'
    files <- data.frame(
        cell1=rep(cell, each=length(categories)),
        dmso=rep(base,length(cell)*length(categories)),
        vs=rep('vs',length(cell)*length(categories)),
        cell2=cell1,
        condition=rep(categories,3)
    )

# Size new graphical window for output
par(mfrow=c(1,3))

venn.diagram(
    x=du145,
    category.names = ,

)
