library(tidyverse)
library(fs)
library(ggplot2)
library(openxlsx)
library(gridExtra)
library(ggpubr)
# Set paths
paths <- list(
    files = 'dgea_output/rds',
    output = 'Intersection_analysis'
)

# Select files
f <- fs::dir_ls(
    path = paths$files,
    regex = '(VCaP|LAPC4)__(R1881|11KT|11OHT).*_DMSO_de.rds'
    )

# Split files
motifs <- sub("^(.*?__.*?)(?:__.*)?\\.rds$", "\\1", fs::path_file(f)) %>%
    sub('_ENZA','',.)


# Load files
data <- map(f, ~ readRDS(.x) %>% dplyr::select('id','log2FoldChange')) %>%
    split(x = ., f = motifs)

# Set plot title
t <- c('Relation des DEG entre les traitements hormonaux\navec et sans enzalutamide ')

# Generate Graph
metadata <- imap(data, ~ {

    # set output directory
    tmp.out <- fs::path(paths$output, .y)
    if(!dir.exists(tmp.out)){dir.create(tmp.out)}

    # Generate ploting data
    p.data <- Reduce(function(x, y) dplyr::inner_join(x,y, by = 'id'), .x) %>%
        na.omit() %>%
        as_tibble() %>%
        column_to_rownames('id') %>%
        'colnames<-'(c('Enza', 'Normal'))

    # Get coloring information
    d <- apply(p.data, 1,  function(row) diff(row))

    # Compute normality test
    gg.qq <- ggpubr::ggqqplot(
        d,
        fill = 'lightgrey',
        title = 'QQ plot des différences entre traitement hormonaux',
        ggtheme = theme_bw(),
        color = "navyblue",
        size = 0.5
        )
    gg.den <- ggpubr::ggdensity(
        d,
        fill = 'lightgrey',
        title = "Density distribution des différences entre traitement hormonaux",,
        ylab = "Density",
        ggtheme = theme_bw()
        )

    # Get PP plot
    empirical_prob <- ppoints(length(d))
    theoretical_prob <- pnorm(sort(d), mean = mean(d), sd = sd(d))
    pp.df <- data.frame(theoretical_prob, empirical_prob)
    gg.pp <- ggplot(data = pp.df)+
        geom_point(mapping = aes(x = theoretical_prob, y = empirical_prob), color = "navyblue", size = 0.5) +
        geom_abline(intercept = 0, slope = 1, color = 'black')+
        labs(
            title = "P-P Plot",
            x = "Theoretical Probabilities",
            y = "Empirical Probabilities"
        )


    pdf(file = fs::path(tmp.out, paste(.y,'normality_plots.pdf', sep = '_')))
    grid.arrange(gg.qq, gg.den, gg.pp, ncol = 1)
    dev.off()

    # Compute normality test with
    d_jittered <- jitter(d)
    ks <- ks.test(d_jittered, "pnorm", mean = mean(d_jittered), sd = sd(d_jittered))
    anderson <- nortest::ad.test(d)
    skew <- moments::skewness(d)
    kurt <- moments::kurtosis(d)

    # Compute difference
    dist <- scale(d)
    p.data$pvalue <- rep('Not.Sig.', length = length(dist))
    p.data$pvalue[abs(dist) >= 1.96] = '<= 0.05'
    p.data$pvalue[abs(dist) >= 2.58] = '<= 0.01'

    # Set subtitle information
    subt.elem <- str_split(.y, '__') %>% unlist()
    sub.t <- paste(
        'Cell:', subt.elem[1],
        '; Hormone:', subt.elem[2],
        'vs DMSO; Mesure: log2(FoldChange)'
    )

    # Plot data
    p <- ggplot(data = p.data) +
        geom_point(mapping = aes(x = Normal, y = Enza, color = pvalue), size = 0.75) +
        scale_color_manual(values = c(
            'Not.Sig.' = 'gray60',
            '<= 0.05' = "firebrick",
            '<= 0.01' = 'navyblue')) +
        labs(title = t, subtitle = sub.t, caption = Sys.time()) +
        xlab(paste(subt.elem[2], 'without Enza')) +
        ylab(paste(subt.elem[2], 'with Enza')) +
        theme_bw() +
        theme(legend.position = "bottom") +
        guides(color = guide_legend(override.aes=list(shape = 15, size = 4)))

    ggsave(filename = paste(.y, 'foldChange_relation.pdf', sep = '_'), path = tmp.out)

    # Get distribution data
    data.frame(
        id = .y,
        Kurtosis = kurt,
        Skewness = skew,
        "KS test stat" = ks$statistic,
        "KS test pvalue" = ks$p.value,
        "AD test stat" = anderson$statistic,
        "AD test pvalue" = anderson$p.value
    )

    }) %>% purrr::list_rbind() %>%
    as_tibble(.) %>%
    tibble::column_to_rownames("id")

# Write stat data to excel file
wb <- openxlsx::createWorkbook()
sheet <- "Metadata"
openxlsx::addWorksheet(wb, sheetName = sheet)
openxlsx::writeDataTable(
    wb,
    sheet = sheet,
    x = metadata,
    xy = c(2,2),
    rowNames = TRUE,
    tableName = "Tests_statistiques_pour_la_distribution_des_différences_expression"
)

openxlsx::saveWorkbook(
    wb,
    file = fs::path(paths$output, 'statistic_data.xlsx'),
    overwrite = TRUE
)
