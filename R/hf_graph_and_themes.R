## _________________________________________________
##
## Helper Functions: Graphs and Themes
##
## Aim: Include in one single file the functions and variables needed to bring a
## consistent theme to the graphs used in the analysis.
##
## Author: Mr. Guillermo Rocamora Pérez
##
## Date Created: 31/50/2023
##
## Copyright (c) Guillermo Rocamora Pérez, 2023
##
## Email: guillermorocamora@gmail.com
## _________________________________________________
##
## Notes:
##
## Many themes will be further modified per graph as required by the situation.
## This script includes only a guideline on the themes to use and other
## variables required to properly run the graphs.
##
## Please contact guillermorocamora@gmail.com for further assistance.
## _________________________________________________

custom_gg_theme <- theme(plot.title = element_text(size = 14, face = "bold"),
                         panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
                         axis.text.x = ggplot2::element_text(color = "black", size = 9, angle = 0, hjust = 0.5),
                         axis.text.y = ggplot2::element_text(color = "black", size = 9),
                         axis.title.x = ggplot2::element_text(face = "bold", size = 11, margin=margin(5,0,0,0)),
                         axis.title.y = ggplot2::element_text(face = "bold", size = 11, margin=margin(0,10,0,0)),
                         panel.grid.minor = element_line(color = "#444444", linewidth  = 0.05, linetype = 2),
                         panel.grid.major.y = element_line(color = "#444444", linewidth  = 0.05, linetype = 2),
                         panel.grid.major.x = element_line(color = "#444444", linewidth  = 0.1),
                         panel.background = element_rect(fill = "#FBFBFB"),
                         legend.title = element_text(size=12),
                         legend.text = element_text(size=10),
                         legend.position = "top",
                         legend.key = element_rect(color="black", linewidth = 0.2),
                         legend.key.size = unit(1, 'lines'),
                         strip.text.x = element_text(color = "black", face = "bold", size = 9),
                         strip.background = element_rect(color="black", linewidth=1, linetype="solid"),
                         strip.text.y = element_text(color = "black", face = "bold", size = 9),
                         plot.margin = margin(0.5, 0.5, 0.2, 0.5, "cm"))

custom_gg_theme_subtitle <- custom_gg_theme + 
  theme(plot.subtitle = element_text(size=12, vjust = 1, color="black", margin = margin(b = -20, t = -5)),
        plot.title = element_text(size = 15, face = "bold"),
        legend.margin=margin(b = -5),
        legend.justification = "right",
        legend.title = element_text(size=11),
        legend.text = element_text(size=10))

novel_labels <- c("novel_acceptor" = "Novel Acceptor", "novel_donor" = "Novel Donor")
novel_labels <- c("novel_donor" = "Novel Donor", "novel_acceptor" = "Novel Acceptor")
novel_mini_labels <- c("donor" = "Novel Donor", "acceptor" = "Novel Acceptor")
sequence_labels <- c("acceptor exon" = "Acceptor exon",
                     "acceptor intron" = "Acceptor intron",
                     "donor exon" = "Donor exon",
                     "donor intron" = "Donor intron")
delta_MES_labels <- c("delta MES 5'ss" = "Delta MES 5'ss", "delta MES 3'ss" = "Delta MES 3'ss")


plotMetadataSubsample <- function(metadata_subsample, level, ref_group = NULL, output_file = "", ratio = 1, dpi = 300){
  p <- ggplot(metadata_subsample, aes(x = !!sym(level), y = RIN)) +
    geom_boxplot() +
    geom_dotplot(aes(fill = Diagnosis), stackratio = 1.1, binaxis = "y", 
                 color = "black", stackdir = "center", stackgroups = T, binpositions = "all", dotsize = 0.7) +
    ggpubr::geom_pwc(label = " p = {p.adj.format}", p.adjust.method = "bonferroni", p.adjust.by = "group",
                     ref.group = ref_group) +
    custom_gg_theme
  
  if(output_file != "") ggsave(output_file, plot = p, width = 183, height = 183*ratio, units = "mm", dpi = 300)
  return(p)
}


addArrows <- function(x_arrow, x_text_offset, y_start, y_end, facet_info, clusters){
  data_arrow <- tibble(x_start = x_arrow, 
                       x_end = x_arrow, 
                       y_start = c(y_start, -y_start), 
                       y_end = c(y_end, -y_end), 
                       facet_info)
  
  data_text <- tibble(x = x_arrow-x_text_offset,
                      y = c((y_start+y_end)/2*0.95,
                            -(y_start+y_end)/2*0.95),
                      label = c(parse(text = paste0('"MSR[', clusters[1], ']', '>', 'MSR[', clusters[2], ']"')),
                                parse(text = paste0('"MSR[', clusters[1], ']', '<', 'MSR[', clusters[2], ']"'))))
  
  list(geom_segment(aes(x = x_start, xend = x_end, y = y_start, yend = y_end), 
                    arrow = arrow(type = "closed", length = unit(0.1, "npc")), linewidth = 0.2,
                    data = data_arrow),
       geom_text(aes(x = x, y = y, label = label),
                 parse = T, size = 3, data = data_text))
}