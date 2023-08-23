library(tidyverse)
library(data.table)
library(ggsci)
library(Seurat)
library(ggpubr)
library(harmony)


num_PC <- 25
scRNA_mouse <- NormalizeData(scRNA_harmony_mouse) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
scRNA_mouse <- RunHarmony(scRNA_mouse, group.by.vars = "orig.ident")
scRNA_mouse <- RunUMAP(scRNA_mouse, reduction = "harmony", dims = 1:num_PC) %>%
                 FindNeighbors(reduction = "harmony", dims = 1:num_PC) %>%FindClusters(resolution = 1.2)

# Idents(scRNA_mouse)

all_gene <- scRNA_mouse@assays$RNA@counts %>% rownames() %>% as.data.frame %>% rename("gene"=1)

########## Fig S5C & Fig 5J ##########

p<-Seurat::DimPlot(
    scRNA_mouse,
    # scRNA_harmony_mouse,
    reduction = "umap", 
    # group.by = "seurat_clusters",
    group.by = "group", order = rev(c("Control","Immunotherapy","Radiotherapy","Radio+Immunotherapy")),
    # label = TRUE,
    label.size = 3,
    pt.size = 0.3
    )+
    # scale_color_npg()+
    scale_color_manual(
        values = c(       
        "#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF"
        # celltype
        # "#B09C85","#7E6148","#ED0000B2","#0072B5","#E18727","#20854E", #T
        # "#7876B1",#B 
        # "#6F99AD",
        # "#BC3C29","#4DBBD5FF","#00A087FF","#3C5488FF",
        # "#ADB6B6", #DC
        # "#0099B4B2", #Pro
        # "#F39B7F","#8491B4","#91D1C2"
    ))+
    # theme_bw()+
     theme(
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid= element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(), 
        #   legend.key.size = unit(25, "pt"),
          legend.title =element_blank(),
          legend.text = element_text(size = 10,  color = "Black")
    )
  
scRNA_mouse <- RenameIdents(scRNA_mouse, 
        
        '0' = 'DNT','22' = 'Tn',
        '16' = 'γδ T',
        '5' = 'Tem',
        '23' = 'Isg+T',
        '15' = 'Treg',

        '13' ='B cell',
        
        '9' = 'Neutrophil','19' = 'Neutrophil',
        '7' = 'Macro_Rsad2',
        '6' = 'Macro_Chil3','12' = 'Macro_Chil3',
        '2' = 'Macro_Retnla',
        '1' = 'Monocyte',
        '17' = 'DC',
        '10' = 'Proliferative',
    
        '11' = 'Endothelial','21' = 'Endothelial',

        '3' = 'Fibroblast', '18' = 'Fibroblast',

        '14' = 'Epithelial', '4' = 'Epithelial', '8' = 'Epithelial','24' = 'Epithelial',  '20' = 'Epithelial'
)


celltypes_col <- Idents(scRNA_mouse) %>% as.data.frame %>% rename("celltype"=1)
scRNA_mouse <- AddMetaData(object = scRNA_mouse, metadata = celltypes_col, col.name = 'celltype')


ident2group <-c(
    "TBP_mice_1"="Control",
    "TBP_mice_3"="Immunotherapy",
    "TBP_mice_4-2"="Radiotherapy",
    "TBP_mice_2"="Radio+Immunotherapy"
    )

scRNA_mouse[['group']] = unname(ident2group[scRNA_mouse@meta.data$orig.ident])

scRNA_mouse$group <- factor(scRNA_mouse$group, levels = c("Control","Immunotherapy","Radiotherapy","Radio+Immunotherapy"))




############# Fig 4K ############# 
library(ggalluvial)

T_cell_ <-  subset(scRNA_mouse, subset = celltype %in% c("DNT","Tn","Tem","Isg+T","γδ T",'Treg'))

T_cell_@meta.data %>% head()

pp<-table(scRNA_mouse@meta.data[, c("group", "celltype")]) / as.vector(table(scRNA_mouse$group))

pp <- pp %>% as.data.frame() %>% melt()

pp$group <- factor(pp$group, c("Control",
            "Radiotherapy",
            "Immunotherapy",
            "Radio+Immunotherapy"))

# prop$Var1 <- factor(prop$Var1, rev(c("T cell","Treg","B cell","Neutrophils","Mono/Macs","DC","Proliferative","Endothelial","Fibroblasts","Epithelial")))

pp$celltype <- factor(pp$celltype, c("DNT","Tn","Tem","Isg+T","Treg","γδ T",
"B cell","Neutrophil","Macro_Rsad2","Macro_Chil3","Macro_Retnla","Monocyte","DC","Proliferative",
"Endothelial","Fibroblast" ,"Epithelial"))

p<-ggplot(pp, aes(x = group, y = value, fill = celltype, 
        label = paste0(round(value*100,digits=2)),
        stratum= celltype, alluvium = celltype)
  ) +
  geom_col(width = 0.5,color=NA)+
geom_flow(width = 0.5, alpha = 0.1)+
  geom_text(position=position_stack(vjust=0.5), color="white", size = 4.5, angle=0) +
 scale_fill_manual(
        values = 
        c(
        "#ad9d88cd","#79624bd1","#3070b0d4", "#d48b3ed7","#418353cd", "#d92f20e0",
        # "#ADCAD7", "#5289A8","#91D1C2","#ADB6B6",
        "#EBEBEB","#EBEBEB","#EBEBEB","#EBEBEB",
        "#EBEBEB","#EBEBEB","#EBEBEB","#EBEBEB","#EBEBEB","#EBEBEB","#EBEBEB","#EBEBEB","#EBEBEB","#EBEBEB","#EBEBEB","#EBEBEB")
        # c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF")
        )+
ylab("")+xlab("")+
  theme_bw()+
  theme(
        panel.grid= element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key.size = unit(25, "pt"),
        legend.title =element_blank(),
        legend.text = element_text(size = 20,  color = "Black"),
        axis.text.y = element_text(size = 20,  angle = 0, color = "Black"),
       axis.text.x= element_text(size = 20,  angle = 90, color = "Black",vjust=0.5,hjust = 1),
      plot.title = element_text(face="bold",size=20,hjust = 0.5)
  )
