library(tidyverse)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(immunarch)
library(ggpubr)
library(ggsci)


############ 50% mean copy number cutoff ############

def_half <- function(immdata){
    immdata_cutoff <- list()
    for (i in  immdata$meta$Sample){
        print(i)
        cutoff_ <- immdata$data[[i]]$Clones %>% sum / nrow( immdata$data[[i]] )
        immdata_cutoff$data[[i]] <- immdata$data[[i]] %>% filter(Clones >= cutoff_/2)
    }
    immdata_cutoff$meta <- immdata$meta
    return(immdata_cutoff)
}

immdata_TRG_pos_cutoff <- def_half(immdata_TRG_pos)
immdata_TRG_neg_cutoff <- def_half(immdata_TRG_neg)

immdata_TRD_pos_cutoff <- def_half(immdata_TRD_pos)
immdata_TRD_neg_cutoff <- def_half(immdata_TRD_neg)




cosine <- data.frame(type = 'ttt')
for (i in immdata_TRG_pos_cutoff$meta$Sample){
    for (j in immdata_TRG_pos_cutoff$meta$Sample){
        print(i)
        print(j)
        tmp <- immdata_TRG_pos_cutoff$data[[i]] %>% select(1,4) %>% 
            group_by(CDR3.aa) %>% summarise(clones = sum(Clones)) %>%
            inner_join(
            immdata_TRG_pos_cutoff$data[[j]]%>% select(1,4)%>% 
            group_by(CDR3.aa) %>% summarise(clones = sum(Clones)), by ="CDR3.aa")
        
            a <- as.vector(tmp$clones.x)
            b <- as.vector(tmp$clones.y)
            #a * b / (||a|| * ||b||)
            s <- a %*% b/ (norm(as.matrix(a),"2") * norm(as.matrix(b),"2"))%>% 
                as.data.frame %>% setnames(paste0(i,"_",j))
        
        cosine <- cbind(cosine, s)

}
}


######### heatmap ######### 
library(circlize)
library(ComplexHeatmap)

d3 <- cosine %>% as.matrix()
rownames(d3)<-colnames(d3)



col_fun = colorRamp2(c(0,0.1,0.3,0.5,0.7,0.9,1), 
c("#04427dce","#3665a7e7", "#96BBD7","#F4F5F6", "#EFB69A","#D3785C","#A42A31"))


pdf('cosine_heatmap_HC.pdf', height=10*0.5, width=11*0.45)
Heatmap(d3[1:7,1:7], 
    col=col_fun,
    show_heatmap_legend =FALSE,
    cluster_rows = FALSE,
    cluster_columns =FALSE,
    # border = TRUE,
    rect_gp = gpar(type = "none"),
    cell_fun = function(j, i, x, y, width, height, fill) {
        if(i >= j) {
			grid.rect(x, y, width, height, gp = gpar(fill = fill, col = "white"))
            grid.text(sprintf("%.1f",d3[1:7,1:7][i, j]), x, y, gp = gpar(fontsize = 14))
		}
        },
    row_names_side = "left",
    # top_annotation = column_ha,
    # column_split = df_$Type,
    # row_split = df_$Type,
    row_names_gp = gpar(fontsize = 10, fontface = "plain"),
    column_names_gp = gpar(fontsize = 10, fontface = "plain")
)
dev.off()



pdf('cosine_heatmap_miPTC.pdf', height=10*0.5, width=11*0.45)
Heatmap(d3[8:15,8:15], 
    col=col_fun,
    show_heatmap_legend =FALSE,
    cluster_rows = FALSE,
    cluster_columns =FALSE,
    # border = TRUE,
    rect_gp = gpar(type = "none"),
    cell_fun = function(j, i, x, y, width, height, fill) {
        if(j >= i) {
			grid.rect(x, y, width, height, gp = gpar(fill = fill, col = "white"))
            grid.text(sprintf("%.1f",d3[8:15,8:15][i, j]), x, y, gp = gpar(fontsize = 14))
		}
        },
    row_names_side = "right",
    column_names_side = "top",
    # top_annotation = column_ha,
    # column_split = df_$Type,
    # row_split = df_$Type,
    row_names_gp = gpar(fontsize = 10, fontface = "plain"),
    column_names_gp = gpar(fontsize = 10, fontface = "plain")
)
dev.off()



pdf('cosine_heatmap_wiPTC.pdf', height=10*0.5, width=11*0.45)
Heatmap(d3[16:27,16:27], 
    col=col_fun,
    show_heatmap_legend =FALSE,
    cluster_rows = FALSE,
    cluster_columns =FALSE,
    # border = TRUE,
    rect_gp = gpar(type = "none"),
    cell_fun = function(j, i, x, y, width, height, fill) {
        if(i >= j) {
			grid.rect(x, y, width, height, gp = gpar(fill = fill, col = "white"))
            grid.text(sprintf("%.1f",d3[16:27,16:27][i, j]), x, y, gp = gpar(fontsize = 14))
		}
    },
    row_names_side = "left",
    # top_annotation = column_ha,
    # column_split = df_$Type,
    # row_split = df_$Type,
    row_names_gp = gpar(fontsize = 10, fontface = "plain"),
    column_names_gp = gpar(fontsize = 10, fontface = "plain")
)
dev.off()



pdf('cosine_heatmap_PDTC.pdf', height=10*0.5, width=11*0.45)
Heatmap(d3[28:34,28:34], 
    col=col_fun,
    show_heatmap_legend =FALSE,
    cluster_rows = FALSE,
    cluster_columns =FALSE,
    # border = TRUE,
    rect_gp = gpar(type = "none"),
    cell_fun = function(j, i, x, y, width, height, fill) {
        if(i >= j) {
			grid.rect(x, y, width, height, gp = gpar(fill = fill, col = "white"))
            grid.text(sprintf("%.1f",d3[28:34,28:34][i, j]), x, y, gp = gpar(fontsize = 14))
		}
        },
    row_names_side = "left",
    # column_names_side = "top",
    # top_annotation = column_ha,
    # column_split = df_$Type,
    # row_split = df_$Type,
    row_names_gp = gpar(fontsize = 10, fontface = "plain"),
    column_names_gp = gpar(fontsize = 10, fontface = "plain")
)
dev.off()



pdf('cosine_heatmap_ATC.pdf', height=10*0.5, width=11*0.45)
Heatmap(d3[35:45,35:45], 
    col=col_fun,
    show_heatmap_legend =FALSE,
    cluster_rows = FALSE,
    cluster_columns =FALSE,
    # border = TRUE,
    rect_gp = gpar(type = "none"),
    cell_fun = function(j, i, x, y, width, height, fill) {
        if(j >= i) {
			grid.rect(x, y, width, height, gp = gpar(fill = fill, col = "white"))
            grid.text(sprintf("%.1f",d3[35:45,35:45][i, j]), x, y, gp = gpar(fontsize = 14))
		}
        },
    row_names_side = "right",
    column_names_side = "top",
    # top_annotation = column_ha,
    # column_split = df_$Type,
    # row_split = df_$Type,
    row_names_gp = gpar(fontsize = 10, fontface = "plain"),
    column_names_gp = gpar(fontsize = 10, fontface = "plain"),
    heatmap_legend_param = list(
        title = "Cosine (%)", 
        # title = "Morisita (%%)", 
        # title_gp = gpar(fontsize = 10, fontface = "bold"),
        # legend_height = unit(6, "cm"),
        # grid_width = unit(0.5, "cm"),
        # legend_height = unit(4, "cm"),
        # grid_height = unit(1, "cm"),
        # # direction =  "horizontal",
        title_position = "topcenter")
)
dev.off()


############ statistic ############

df <- rbind(
cosine[1:7,1:7]    %>%as.matrix() %>% `[`(lower.tri(.)) %>% as.data.frame() %>% mutate(type = "HC"),
cosine[8:15,8:15]  %>%as.matrix() %>% `[`(lower.tri(.))%>% as.data.frame() %>% mutate(type = "miPTC"),
cosine[16:27,16:27]%>%as.matrix() %>% `[`(lower.tri(.))%>% as.data.frame() %>% mutate(type = "wiPTC"),
cosine[28:34,28:34]%>%as.matrix() %>% `[`(lower.tri(.))%>% as.data.frame() %>% mutate(type = "PDTC"),
cosine[35:45,35:45]%>%as.matrix() %>% `[`(lower.tri(.))%>% as.data.frame() %>% mutate(type = "ATC")
# TRG_pos_p[46:54,46:54]%>%as.matrix() %>% `[`(lower.tri(.))%>% as.data.frame() %>% mutate(type = "HC")
) %>% rename("value"=1)


my_comparisons <- list( 
                        c("ATC", "HC"), 
                        c("ATC", "miPTC"),
                        c("ATC", "wiPTC") ,
                        c("ATC", "PDTC")
                        )

df$type <- factor(df$type,levels=c("HC","miPTC","wiPTC","PDTC","ATC"))



p<-ggplot(df, aes(x=type, y=value, fill=type))+
  geom_jitter(
        # position = position_dodge(0.3),
        position = position_jitterdodge(
            dodge.width = 0.3,
            jitter.width = 1),
       shape = 21,size = 3,color="black")+  
  geom_pointrange(stat="summary", fun.data="mean_sdl",fun.args = list(mult=1),
                  color = "black",size = 1)+
  geom_point(stat="summary", fun.y="mean",fun.args = list(mult=1),
                color = "white",size = 3)+
  stat_compare_means(method = "kruskal.test", label.x = 1,
    aes(label = paste0("p = ", ..p.format..)))+
    scale_color_npg(alpha = 1) +
    scale_fill_npg(alpha = 0.8) +
    # ylab("Overlap of Vg9-")+
    ylab("Cosine similairty of Vd2+")+
    xlab("")+
    theme_bw()+
    theme(
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 1),
        plot.margin = unit(c(1, 0.5, 0.5, 1), "cm"),
        axis.ticks = element_line(linewidth=0.75, color = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        # axis.ticks.x = element_blank(),
        axis.text.x=element_text(face="plain",size=12,angle=0,color="Black",vjust=1),
        axis.text.y=element_text(face="plain",size=12,angle=0,color="Black"),
        axis.title.y = element_text(size = 15, face="plain", color = "Black"),
        axis.title.x = element_text(size = 15, face="plain", color = "Black"),
        panel.grid= element_blank(),
        legend.position = 'none'
    )+
    expand_limits(y = c(0,1))+ 
    scale_y_continuous(expand = c(0, 0))




