suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggridges))
library(ggplot2)

## Load method colors as a named vector
col   <- unlist(yaml::yaml.load_file('code/helpers/colors.yml'))

## Unlist correlations/ convert to data frames for ggplot2/ggridges
list2df <- function(est.list) {
  est.list<-est.full_sub 
  reshape2::melt(est.list) %>%
    dplyr::rename(dataset=L1) %>%
    bind_cols(
      stringr::str_split_fixed(.$dataset, "_|\\.", 4) %>%
        as.data.frame(stringsAsFactors=FALSE) %>%
        mutate(rep=as.numeric(V1),
               n=as.numeric(V2),
               method=V3) %>%
        select(rep, n, method)
    )
}

mainord <- yaml::yaml.load_file('code/helpers/data_order.yml')[['main']]
suppord <- yaml::yaml.load_file('code/helpers/data_order.yml')[['supp']]

plotDensity <- function(df, main="", levels=order) {
  ggplot(aes(x=value, y=factor(method, levels=rev(levels)), fill=method), data=df %>% filter(method %in% levels)) +
    facet_wrap(~group) + #,scales = "free_x") +
    geom_density_ridges(scale=1.3) +
    scale_fill_manual(values=col) +
    guides(fill=FALSE) + ylab("") +
    scale_x_continuous(limits=c(-0.25,0.25), breaks=c(-0.2,-0.1,0,0.1,0.2)) + # comment out for free scales- see full tails
    ggtitle(main) +
    theme_ridges(font_size = 16,center_axis_labels = TRUE) + theme_bw() + 
    theme(axis.title.x=element_blank()) + theme(axis.text.y = element_text(size = 15))
}

