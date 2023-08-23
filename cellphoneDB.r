library(tidyverse)
library(data.table) 
library(ggsci)
library(Seurat)
library(ggpubr)

scRNA_harmony <- readRDS("scRNA_harmony.vd2.celltype.rds")

d_pre <- subset(scRNA_harmony, subset =  group == "Pre")
d_post <- subset(scRNA_harmony, subset = group == "Post")

write.table(as.matrix(d_pre@assays$RNA@data), 'cellphonedb/pre/cellphonedb_count.txt', sep='\t', quote=F)
write.table(as.matrix(d_post@assays$RNA@data), 'cellphonedb/post/cellphonedb_count.txt', sep='\t', quote=F)

meta_data <- cbind(rownames(d_pre@meta.data), d_pre@meta.data[,'celltype_a', drop=F]) 
meta_data <- as.matrix(meta_data)
write.table(meta_data, 'cellphonedb/pre/cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)

meta_data <- cbind(rownames(d_post@meta.data), d_post@meta.data[,'celltype_a', drop=F]) 
meta_data <- as.matrix(meta_data)
write.table(meta_data, 'cellphonedb/post/cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)




net1 <- fread("cellphonedb/pre/out/count_network.txt") %>% 
    filter(SOURCE != "ATC"  & SOURCE != "Thyrocyte" &
    SOURCE != "Fibroblast" &  SOURCE != "Endothelial"&  SOURCE != "B cell"
     & SOURCE != "Erythrocyte" & SOURCE != "Mast cells") %>%
    filter(TARGET != "ATC" & TARGET != "Thyrocyte" &
    TARGET != "Fibroblast" & TARGET != "Endothelial"&  TARGET != "B cell"
     & TARGET != "Erythrocyte"& TARGET != "Mast cells")


hp <- net1 %>% select(count)%>% log() %>% as.matrix()
dim(hp) <- c(14,14)
# dim(hp) <- c(11,11)
colnames(hp) =  net1 %>% select(SOURCE) %>% unique %>% unlist()
rownames(hp) =  net1 %>% select(SOURCE) %>% unique %>% unlist()

net2 <- fread("cellphonedb/post/out/count_network.txt") %>% 
   filter(SOURCE != "ATC"  & SOURCE != "Thyrocyte" &
    SOURCE != "Fibroblast" &  SOURCE != "Endothelial"&  SOURCE != "B cell"
     & SOURCE != "Erythrocyte" & SOURCE != "Mast cells") %>%
    filter(TARGET != "ATC" & TARGET != "Thyrocyte" &
    TARGET != "Fibroblast" & TARGET != "Endothelial"&  TARGET != "B cell"
     & TARGET != "Erythrocyte"& TARGET != "Mast cells")


hp <- net2 %>% select(count) %>% log() %>% as.matrix()
 hp[hp<0] <- 0
dim(hp) <- c(13,13)
# dim(hp) <- c(11,11)
colnames(hp) =  net2 %>% select(SOURCE) %>% unique %>% unlist()
rownames(hp) =  net2 %>% select(SOURCE) %>% unique %>% unlist()


############## Fig 4E ##############

library(circlize)
library(ComplexHeatmap)

col_fun = colorRamp2(c(0,20,40,60,90),c("#5A79A6","#B1B6BC","#E9DACA","#C892A0","#A0567C"))
col_fun = colorRamp2(c(0.5,1.5,2.5,3.5,4.5),c("#5D7CA6",
# "#95A4B5","#DCD0C7","#E5C3BB","#C28A9C","#A15A7F"
"#B1B6BC","#E9DACA","#C892A0","#A0567C"
))

col_fun = colorRamp2(c(0.5,1.5,2.5,3.5,4.5),c("#3f0751e2","#432671b5","#468b8cae","#7ec96ba7","#f9e855b3"))

pdf('test_heatmap.pdf', height=15*0.5, width=16*0.5)
Heatmap(hp,
        row_names_gp = gpar(fontsize = 18),
        row_names_side  = "right",
        col = col_fun,
        gap = unit(40, "cm"),
        rect_gp = gpar(col = "white"),
        column_names_gp = gpar(fontsize = 18),
         heatmap_legend_param = list(
        title = "log(Count)", 
        title_gp = gpar(fontsize = 10, fontface = "bold"),
        # legend_height = unit(6, "cm"),
        # grid_width = unit(0.5, "cm"),
        legend_height = unit(5, "cm"),
        grid_width = unit(0.5, "cm"),
        # direction =  "horizontal")
        title_position = "leftcenter-rot"
        # title_position = "topcenter"
        ))

dev.off()


############## Fig 4G ##############

meansdf_post <- fread("cellphonedb/post/out/means.txt") %>% 
  dplyr::select(interacting_pair,"Vδ2-|ISG+ T","ISG+ T|Vδ2-") %>%
  reshape2::melt()
colnames(meansdf_post)<- c("interacting_pair","CC","means")

pvalsdf_post <- fread("cellphonedb/post/out/pvalues.txt") %>% 
  dplyr::select(interacting_pair,"Vδ2-|ISG+ T","ISG+ T|Vδ2-") %>%
  reshape2::melt()
colnames(pvalsdf_post)<- c("interacting_pair","CC","pvals")

meansdf_post$joinlab<- paste0(meansdf_post$interacting_pair,"_",meansdf_post$CC)
pvalsdf_post$joinlab<- paste0(pvalsdf_post$interacting_pair,"_",pvalsdf_post$CC)
pldf_post <- merge(pvalsdf_post,meansdf_post,by = c("joinlab","interacting_pair","CC"))

meansdf_pre <- fread("cellphonedb/pre/out/means.txt") %>% 
  dplyr::select(interacting_pair,"Vδ2-|ISG+ T","ISG+ T|Vδ2-") %>%
  reshape2::melt()
colnames(meansdf_pre)<- c("interacting_pair","CC","means")

pvalsdf_pre <- fread("cellphonedb/pre/out/pvalues.txt") %>% 
  dplyr::select(interacting_pair,"Vδ2-|ISG+ T","ISG+ T|Vδ2-") %>%
    reshape2::melt()
colnames(pvalsdf_pre)<- c("interacting_pair","CC","pvals")

pvalsdf_pre$joinlab<- paste0(pvalsdf_pre$interacting_pair,"_",pvalsdf_pre$CC)
meansdf_pre$joinlab<- paste0(meansdf_pre$interacting_pair,"_",meansdf_pre$CC)
pldf_pre <- merge(pvalsdf_pre,meansdf_pre,by = c("joinlab","interacting_pair","CC"))

df2 <- rbind(pldf_pre %>% mutate(time="Pre"),pldf_post%>% mutate(time="Post"))
df1 <- pldf_pre %>% mutate(time="Pre") %>% 
  right_join(pldf_post %>% mutate(time="Post"), by = c("joinlab","interacting_pair","CC")) %>%
  filter(means.y > means.x | pvals.y < pvals.x) %>%
  filter(means.y >= 0.5 & pvals.y < 0.05) %>%
  drop_na()

new_df <- rbind(df1 %>% select(1:6) %>% rename("pvals"="pvals.x","means"="means.x","time"="time.x"),
    df1 %>% select(1:3,7:9) %>% rename("pvals"="pvals.y","means"="means.y",,"time"="time.y"))

new_df$time <- factor(new_df$time, levels = c("Pre","Post"))


p<-new_df %>% 
    filter(joinlab != "PTGER2_ProstaglandinE2_byPTGES3_Vδ2-|ISG+ T"&
    joinlab !="PTGER4_ProstaglandinE2_byPTGES3_Vδ2-|ISG+ T" &
    joinlab!= "ICAM3_integrin_aLb2_complex_ISG+ T|Vδ2-") %>%
    ggplot(aes(time,interacting_pair) )+
    xlab("")+ylab("")+
    geom_point(aes(color=means,size=-log10(pvals+0.0001)) ) +
    scale_size_continuous(range = c(1,5))+
    scale_color_gradient2(high="#FFCC33",mid = "#66CC66",low ="#330066",midpoint = 0.7  )+ theme_bw()+ 
    theme(axis.text.x = element_text(angle = 0
    # ,hjust = 0.1,vjust = 0.8
    ),
    legend.key.size = unit(25, "pt"),
            # legend.title =element_blank(),
            legend.text = element_text(size = 15,  color = "Black"),
        axis.text= element_text(size = 15,  color = "Black"))

