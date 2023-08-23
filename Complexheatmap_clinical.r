library(tidyverse)
library(data.table) 
library(ComplexHeatmap)


###### Fig2A
# df = div_trd or div_trg

hp <- df %>% select(shannon_entropy) %>% t() %>% as.matrix()
colnames(hp) = df %>% select(Sample) %>% unlist


col_fun = colorRamp2(c(0,8),c("white","red"))
col_fun1 = colorRamp2(c(0,0.8),c("white","#3385c8"))



column_ha1 = HeatmapAnnotation(
        # empty = anno_empty(border = FALSE),
        foo = anno_block(gp = gpar(fill = c("#F06869","#7775CF","#3C97CD","#4AB79E","#B6B6B9"))
        ),
        # Type = df$Type,
        Sex = df$Sex,
        Survival = df$Surv,
        `Distant metastasis` = df$`Distant metastasis`,
        `Shannon entropy` = df$shannon_entropy_G,
        Clonality = df$Clonality_G,
        ##bar
        `Unique TCRs` = anno_barplot(df$unique.clonotypes_G,bar_width = 0.6,
                            axis_param = list(side = "right"),beside = TRUE,attach = TRUE),
        Fraction = anno_barplot(df %>% select(`T10000__G`,`T1000__G`,`T100__G`,`T10__G`), 
            gp = gpar(fill = c("#4DBBD5FF","#00A087FF","#3C5488FF","#6f99ad"), 
                      col = c("#4DBBD5FF","#00A087FF","#3C5488FF","#6f99ad"),alpha=0.9), 
            bar_width = 0.6,axis_param = list(side = "right",at = c(0,0.25,0.50,0.75,1))),
        annotation_height = unit(c(0.2,0.3,0.3,0.3,0.3,0.3,1,4),'cm'),
        col = list(
          Sex = c("male"="#2068bbc3", "female"="#ee3017c9"),
          `Shannon entropy` = col_fun,
          Clonality = col_fun1,
          Survival = c("1"="#256613", "0"="#46AF58", "NA"="#AFD093"),
          `Distant metastasis` = c("Present"="#46ADEC", "Absent"="#A6DCE5", "/"="#e9dede")
          ),
        show_legend = FALSE,
        gap = unit(4, "points"),
        annotation_name_side = "left",
        annotation_name_gp = gpar(fontface = "bold",fontsize = 12),
        gp = gpar(col = "white"))




panel_fun = function(index, nm) {
    pushViewport(viewport(yscale = range(df$total.hit_clones_G), xscale = c(0, 2)))
    grid.rect()
    grid.xaxis(main = FALSE, gp = gpar(fontsize = 6))
    grid.boxplot(df[index,"total.hit_clones_G"], 
      pos = 1, 
      gp = gpar(fill = "#D3D3D3"))
    popViewport()
}

bar = HeatmapAnnotation(
  foo = anno_link(
    align_to = df$Type, which = "row", panel_fun = panel_fun,
    side = "right",
    size = unit(1, "cm"), 
    link_gp = gpar(fill = c("#F06869","#7775CF","#3C97CD","#4AB79E","#B6B6B9"), alpha=0.2)
  ))


lgd_list = list(
    Legend(labels = c("Female","Male"), legend_gp = gpar(fill = c("#ee3017c9","#2068bbc3")), title = "Sex", nr = 1,grid_width = unit(3, "mm"),grid_height = unit(5, "mm")),
    Legend(labels = c("Death","Alive","FU loss"), legend_gp = gpar(fill = c("#256613","#46AF58","#AFD093")), title = "Survival", nr = 1,grid_width = unit(3, "mm"),grid_height = unit(5, "mm")),
    Legend(labels = c("Present","Absent","N/A"), legend_gp = gpar(fill = c("#46ADEC","#A6DCE5","#e9dede")), title = "Distant metastasis", nr = 1,grid_width = unit(3, "mm"),grid_height = unit(5, "mm")),
    # Legend(col_fun = col_fun2, title = "Infection", at = c(0, 2000, 4000, 6000),direction = "horizontal",legend_width = unit(4, "cm")),
    Legend(col_fun = col_fun, title = "Shannon entropy", at = c(0, 2, 4, 6, 8),direction = "horizontal",legend_width = unit(4, "cm")),
    Legend(col_fun = col_fun1, title = "Clonality", at = c(0, 0.2, 0.4, 0.6, 0.8),direction = "horizontal",legend_width = unit(4, "cm")),
    Legend(labels = c("[1,10]", "[11,100]","[101,1000]","[1001,)"), title = "Fraction",
        # legend_gp = gpar(fill = c("#FEDC7B","#49B498","#4C5997","#73407a"))
        legend_gp = gpar(fill = c("#6f99ad","#3C5488FF","#00A087FF","#4DBBD5FF"))
        )
)

pdf('test_heatmap.pdf', height=20*0.4, width=30*0.2 )

Heatmap(hp, 
        column_split = df$Type,
        top_annotation = column_ha1,
        bottom_annotation = bar,
        col = col_fun,
        height = unit(0, "cm"),
        show_heatmap_legend = FALSE,
        row_names_gp = gpar(fontface = "bold",fontsize = 12),
        row_names_side  = "left",
        column_names_gp = gpar(fontsize = 10),
         heatmap_legend_param = list(
        title = " ", 
        direction =  "horizontal")
        )
        
draw(ht,annotation_legend_list = lgd_list)

dev.off()

