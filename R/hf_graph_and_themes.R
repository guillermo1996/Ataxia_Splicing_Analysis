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

