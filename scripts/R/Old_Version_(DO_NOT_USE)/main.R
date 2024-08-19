## Load all R scripts ----
s_dir <- './scripts/R'
for (script in file.path(s_dir, list.files(path = s_dir, recursive = TRUE))){
    if (script != './scripts/R/main.R') source(script)
}
rm(s_dir, script)
source("../rnaseq/R/volcano.R")

# Get data and generate required data structures ----
path_kallistos <- './deseq2_files'
util_files <- 'files_util'
component <- list(
    cell=c('DU145', 'VCaP','LAPC4'),
    hormone=c('DMSO', 'R1881', '11KT', '11OHT')
)

create_directories(util=util_files)
make_comparison_file(path=util_files, comp=component)
comp <- read.csv(file.path(util_files,'combinaison.csv'))
rm(component)

for (i in 40:nrow(comp)) {
    print(paste('Row:', i, sep=' '))
    comparison <- as.character(comp[i, ])
    dgea(
        comparison=comparison,
        kallisto_path=path_kallistos,
        util=util_files
        )
}




