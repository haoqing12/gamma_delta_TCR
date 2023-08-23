library(tidyverse)
library(data.table)
library(ggsci)
library(Seurat)
library(ggpubr)
library(SingleR)
# library(celldex)
library(harmony)





#=============Harmony流程=============
num_PC <- 30
scRNA_harmony <- NormalizeData(scRNA_) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident")
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:num_PC) %>%
                 FindNeighbors(reduction = "harmony", dims = 1:num_PC) %>%FindClusters(resolution = 1)


scRNA_harmony.markers <- FindAllMarkers(scRNA_harmony, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25)


Idents(scRNA_harmony) <- scRNA_harmony@meta.data$seurat_clusters


scRNA_harmony <- RenameIdents(scRNA_harmony, 
        '2' = 'γδ T',
        '0' = 'Tem',
        '1' = 'Tex','4' = 'Tex',
        '6' = 'Tn',
        '12' = 'Treg',
        '14' = 'ISG+ T',
        '9' = 'Proliferative T',

        '19' ='B cell',
        '17' = 'pDC',
        '20' = 'cDC',
        '8' = 'Mono_CD14',
        '7' = 'Mono_CD16', 
        '15' = 'Macro_IL1RN', 
        '13' = 'Macro_C1QC','5' = 'Macro_C1QC',
        '22' = 'Mast cells',

        '18' = 'Endothelial',
        '11' = 'Fibroblast', '16' = 'Fibroblast',
        '23' = 'Thyrocyte',
        '10' = 'ATC', '3' = 'ATC',
        '21' = 'Erythrocyte'
)


scRNA_harmony$celltype <- NULL
celltypes_col <- Idents(scRNA_harmony) %>% as.data.frame %>% rename("celltype"=1)
scRNA_harmony <- AddMetaData(object = scRNA_harmony, metadata = celltypes_col, col.name = 'celltype')
Idents(scRNA_harmony) <- scRNA_harmony$celltype


######### Fig 5A #########

p<-Seurat::DimPlot(scRNA_harmony,
    reduction = "umap", 
    # group.by = "seurat_clusters",
    group.by = "celltype",
    label = TRUE,
    label.size = 3,
    pt.size = 0.25,
    )+
    scale_color_manual(values = c(
        "#DD704D",
         "#D7E9C8","#A5C899", "#FFDC91B2",
        "#5289A8","#7ECADA","#ADCAD7",
         "#539C4C", #B
        "#0099B4B2", "#7E6148FF", #DC
        "#E1A83A","#E2DBA1", # Mono
        "#D977A6","#91D1C2FF", #Macro
        "#88679C","#C7AEC3", 
        "#B29671", "#FDAF91",
    "#E64B35FF",
    "#ADB6B6B2" , "#00468BB2",#useless
    "#1B1919B2"
    ))+
    theme_bw()+
    theme(
    axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid= element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        # legend.key.size = unit(25, "pt"),
        legend.title =element_blank(),
        legend.text = element_text(size = 10,  color = "Black"),
        # axis.text.y = element_text(size = 10,  angle = 0, color = "Black"),
    #    axis.text.x= element_text(size = 10,  angle = 0, color = "Black"),
    #   plot.title = element_text(face="bold",size=10,hjust = 0.5)
  )


######### Fig S5A #########

scRNA_harmony_T <- subset(scRNA_harmony, subset =  celltype %in% c("γδ T","Tem","Tex","Tn","Treg","ISG+ T","Proliferative T"))

p<-VlnPlot(scRNA_harmony_T,
        features = c(
          "CD3G","CD4","CD8A","TRDC"),
        ncol=1, 
        pt.size = 0,
        split.by = 'group',
        group.by = 'celltype',
        cols = c("#BC3C29","#0072B5")
        )



######### Fig 5B #########
p<-DotPlot(
    scRNA_harmony,
    features = c(
      "GNLY","TRDC","KLRF1",
      "GZMA","CCL5","GZMK",
      "IFNG","CXCL13","XCL2",
      "IL7R","LEF1","CCR7",
      "CTLA4","IL2RA","FOXP3",
      "ISG20","IFI44L","MX1",
      "STMN1","PCLAF","MKI67",
      "CD79A","CD79B","IGHM",
    #myeloid
      "JCHAIN","IRF7","PTGDS",
      "HLA-DQA2","HLA-DQA1","CD1C",

      "LYZ","S100A8","S100A9",
      "FCGR3A","LILRB2","PILRA",

      "CCL8","APOBEC3A","IL1RN",
      "C1QB","C1QC",'C1QA',

      "TPSB2","TPSAB1","HPGDS",

      #endo
      "PECAM1","VWF","CDH5",
      "ACTA2","SPARCL1","MGP",
      "SLPI","TG","LCN2",
      "AXL","PBK","UBE2C",
      # "FLNA","IGFBP5","PKN1",
      # "CCND1","GNG11",
      "HBB","HBA2","HBA1"
    )
    )+
    coord_flip()+
    theme_bw()+
    theme(
        panel.grid = element_blank(), 
        axis.text.y=element_text(size=12,angle=0,color="Black",face="italic"
        ,hjust=1,vjust=0.5
        ),
        axis.text.x=element_text(size=12,angle=90,color="Black",hjust=1,vjust=0.5)
    )+
    labs(x=NULL,y=NULL)+
    guides(size=guide_legend(order=3))+
    scale_color_gradientn(
        values = seq(-0.5,1.5,0.5),
        colours = c('#F6FBFD','#C1D4E6','#8F99C6','#6E3075','#834B9A')
    )


######### Fig 5C #########

temp <- subset(scRNA_harmony, subset = celltype == "γδ T")

num_PC <- 5
scRNA_T <- NormalizeData(temp) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
scRNA_T <- RunHarmony(scRNA_T, group.by.vars = "orig.ident")
scRNA_T <- RunUMAP(scRNA_T, reduction = "harmony", dims = 1:num_PC) %>%
                 FindNeighbors(reduction = "harmony", dims = 1:num_PC) %>%FindClusters(resolution = 0.2)

scRNA_T <- RenameIdents(scRNA_T, 
        '2' = 'Vδ2+',
        '0' = 'Vδ2-', '1' = 'Vδ2-')

p<-Seurat::DimPlot(scRNA_T,
      #  reduction = "umap", 
      # group.by = "seurat_clusters",
      pt.size = 1,label.size = 3,
      label = FALSE
      )+
      scale_color_nejm()+
      theme(
        legend.title =element_blank(),
        axis.text = element_text(size = 7,  color = "Black"),
        axis.title = element_text(size = 7,  color = "Black"),
  )

######### Fig S5B #########
p<-VlnPlot(scRNA_T, features = c("KLRC1","FCER1G","TRGC2","IL32"),
  ncol=1,
  pt.size = 0,
  cols = c("#BC3C29B2","#0072B5B2"))

p<-FeaturePlot(scRNA_T, features = c("KLRC1","GZMK","IL32","FGFBP2"),ncol=2,pt.size=1.2)




barcode <- Idents(scRNA_T) %>% as.data.frame %>% rename("celltype_"=1)%>% rownames_to_column()
temp_ <- scRNA_harmony@meta.data  %>% rownames_to_column() %>% 
    left_join(barcode, by = "rowname")

t1 <- temp_ %>% filter(celltype == "γδ T") %>% mutate(celltype_a = celltype_)
t2 <- temp_ %>% filter(celltype != "γδ T") %>% mutate(celltype_a = celltype)

temp_ <- rbind(t1,t2) %>%select(rowname,celltype_a)
rownames(temp_) <- temp_$rowname
temp_ <- temp_ %>% select(-rowname)

scRNA_harmony <- AddMetaData(object = scRNA_harmony, metadata = temp_, col.name = 'celltype_a')
saveRDS(scRNA_harmony,'seurat/scRNA_harmony.vd2.celltype.rds')

 ############# Fig 5C KIR2DL1 and KLRC2 expression #############

scRNA_gdT <- subset(scRNA_harmony, subset =  celltype %in% c("γδ T"))

expression_matrix <- GetAssayData(scRNA_gdT, assay = "RNA",slot = "counts")

expression_matrix <- expression_matrix  %>% as.data.frame %>% 
            rownames_to_column() %>% filter(rowname %in% cyto) %>%
             t() %>% as.data.frame()

colnames(expression_matrix) <- expression_matrix[1,]
expression_matrix <- expression_matrix[-1,]
     
df <- expression_matrix %>% rownames_to_column() %>% left_join(
    scRNA_gdT@meta.data %>% select(group,celltype_a) %>% rownames_to_column(),
    by = "rowname")

i <- 'KIR2DL1'
i <- 'KLRC2'

p<-ggboxplot(df, x = "celltype_a", 
    y = i,color = "celltype_a",
    width = 0.4,
    add = "jitter",
    add.params = list(alpha = 0.6, size = 3),
    bxp.errorbar = T,
    bxp.errorbar.width = 0.2)+
    scale_color_manual(values = c("#f8e1a4","#84B8D5"))+
    # ylim(0,0.5)+
    stat_compare_means(method = "wilcox.test", label.x = 1.5,
        aes(label = ..p.signif..)
    )+
  ylab("Score")+xlab("")+
  facet_wrap(~group)+
  labs(title = paste0(i," expression in γδ subset"))+
#   labs(title ="Vδ2- Ligand in ISG+ T cells")+
  guides(color = "none") +
  theme_bw()+
  theme(
    panel.border = element_blank(),
    plot.title = element_text(face="bold",size=13,hjust = 0.5),
    axis.line = element_line(color = "black", size = 1),
    plot.margin = unit(c(1, 0.5, 0.5, 1), "cm"),
    axis.ticks = element_line(linewidth=1, color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    # axis.ticks.x = element_blank(),
    axis.text.x=element_text(face="plain",size=12,angle=0,color="Black",vjust=1),
    axis.text.y=element_text(face="plain",size=12,angle=0,color="Black"),
    axis.title.y = element_text(size = 15, face="plain", color = "Black"),
    axis.title.x = element_text(size = 15, face="plain", color = "Black"),
    panel.grid= element_blank()
    )


##### STEP2 Pre vs Post #####

##### Fig 5D #####

scRNA_harmony_T <- subset(scRNA_harmony, subset =  celltype %in% c("γδ T","Tem","Tex","Tn","Treg","ISG+ T"
# ,"Proliferative T"
))

library(ggalluvial)

prop$Var1 <- factor(prop$Var1, levels=rev(c("Tem","Tex","Tn","Treg","ISG+ T","Vδ2+","Vδ2-")))
prop$Var2 <- factor(prop$Var2, levels=rev(c("Pre","Post")))

p<-ggplot(prop, aes(x = Var2, y = Freq, fill = Var1, 
        label = paste0(round(Freq*100,digits=2)),
        stratum= Var1, alluvium = Var1)
  ) +
  geom_col(width = 0.45,color=NA)+
  geom_stratum(width=0.45, size = 0.5, alpha = 0.7, color="white") +
geom_flow(width = 0.45, alpha = 0.2)+
#   geom_alluvium(width = 0.4,alpha = 0.2,knot.pos = 0) + #与flow效果相同
  geom_text(position=position_stack(vjust=0.5), color="white", size = 4.5, angle=0) +
   scale_fill_manual(values = rev(c("#82C1C1","#497F7C","#53798E","#D5BBA1","#B09C85FF","#B36F6D","#FF5940")))+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid= element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key.size = unit(25, "pt"),
        legend.title =element_blank(),
        legend.text = element_text(size = 15,  color = "Black"),
       axis.text= element_text(size = 15,  color = "Black"),
      legend.position="bottom",
      plot.title = element_text(face="bold",size=20,hjust = 0.5)
  )+coord_flip()
  guides(fill=guide_legend(nrow = 1))   
p
# ggsave(p, filename = "plot/t1.pdf", width = 5*1, height = 7*1.2)
ggsave(p, filename = "plot/t1.pdf", width = 8*1, height = 3*1.2)



##### Fig 5F #####

ligand_Vd2_neg <- c(
  "CD1D","CD1C",
  "MR1",
  "PROCR",
  "MICA","MICB",
  # "ULBP1","ULBP2",
  "ANXA2","EPHA2",
  "BTNL3"
)

ligand_Vd2_pos <- c(
  "MSH2",
  "BTN2A1","BTN3A1"
)

scRNA_harmony <- AddModuleScore(scRNA_harmony, features = list(ligand_Vd2_neg), name="ligand_Vd2_neg")
scRNA_harmony <- AddModuleScore(scRNA_harmony, features = list(ligand_Vd2_pos), name="ligand_Vd2_pos")


d1 <- subset(scRNA_harmony, subset =  group == "Pre")
d2 <- subset(scRNA_harmony, subset = group == "Post")


new_f <- rbind(
    data.frame(d1$ISG1) %>% setnames("value") %>% mutate(variable = "Pre"),
    data.frame(d2$ISG1) %>% setnames("value")  %>% mutate(variable = "Post")
    )


p<-ggboxplot(new_f, x = "variable", 
    y = "value",color = "variable",
    width = 0.4,
    # palette = c("#00AFBB", "#E7B800"),
    add = "jitter",
    add.params = list(alpha = 0.6, size = 3),
    bxp.errorbar = T,
    bxp.errorbar.width = 0.2)+
scale_color_manual(values = c("#f8e1a4","#84B8D5"))+
# ylim(0,0.5)+
stat_compare_means(method = "wilcox.test", label.x = 1.5,
#   ,label.y = 0.45
    # label.y = 0.85,
    aes(label = ..p.signif..)
  )+
  ylab("Score")+xlab("")+
  # labs(title = paste0("γδ T cell features score in ", i))+
  labs(title ="Vδ2- Ligand in ISG+ T cells")+
  guides(color = "none") +
  theme_bw()+
  theme(
    plot.title = element_text(face="bold",size=14,hjust = 0.5),
    # plot.subtitle = element_text(face="bold",size=15,hjust = 0.5),
    legend.title=element_blank(),
    panel.border = element_rect(size=1.5,color = "black"),
    # legend.position = c(0.8,0.8),
    panel.grid= element_blank(),
    legend.text=element_text(colour="black", size = 12, face = "bold"),
    # legend.key.size = unit(35, "pt"),
    axis.text=element_text(size=12,angle=0,color="Black", vjust=0.5),
    axis.title=element_text(face="bold",size=15,angle=0,color="Black")
    )
