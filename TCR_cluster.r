library(tidyverse)
library(data.table)
library(kernlab)
library(igraph)
library(ggraph)
library(stringr)
# library(colormap)

library(wesanderson)
library(reshape2)
# library(stringdist)

######   cluster   ######

getplot <- function(immdata,type,cutoff){

for (sample1 in immdata$meta$Sample){

      t1 <- immdata$data[[sample1]] %>% select(1,2,4) %>%
      distinct(CDR3.aa)
# slice_sample(n = 200)

stringkern <- stringdot(type = 'spectrum', length = 3, normalized = TRUE)


q <- kernelMatrix(stringkern,t1$CDR3.aa)
colnames(q) <- t1$CDR3.aa
rownames(q) <- t1$CDR3.aa



connect <- q %>% as.data.frame() %>% rownames_to_column() %>% rename("from" = 1) %>%
  tidyr::gather(key="to", value="value", -1) %>%
  # arrange(value) %>% distinct(value, .keep_all = TRUE) %>% 
  # dplyr::filter(value >= 0.75)
  dplyr::filter(value >= cutoff)
  # dplyr::filter(value < 1)


c(as.character(connect$from), as.character(connect$to)) %>%
  as_tibble() %>%
  group_by(value) %>%
  summarize(n=n()) -> vertices
colnames(vertices) <- c("name", "n")

mygraph <- graph_from_data_frame(connect, vertices = vertices, directed = FALSE )
# plot(mygraph)

com <- walktrap.community(mygraph)

vertices_nofilter <- vertices %>% 
  mutate(group = com$membership) %>%
  arrange(group,desc(n)) %>%
  mutate(name=factor(name, name))

vertices_filter <- vertices %>% 
  mutate(group = com$membership) %>%
  filter(group<=100) %>%
  arrange(group,desc(n)) %>%
  mutate(name=factor(name, name))

# vertices <- vertices %>% 
#   mutate(group = com$membership) %>%
#   mutate(group=as.numeric(factor(group,
#                                  levels=sort(summary (as.factor(group)),index.return=TRUE,decreasing = T)$ix,
#                                  order=TRUE)))%>%
#   filter(group<100) %>%
#   arrange(group,desc(n)) %>%
#   mutate(name=factor(name, name))

connect <- connect %>%
  filter(from %in% vertices_filter$name) %>%
  filter(to %in% vertices_filter$name)%>%
 left_join(vertices_filter,by=c('from'='name'))

# Create a graph object with igraph
mygraph <- graph_from_data_frame(connect, vertices = vertices_filter, directed = FALSE )
mycolor <- wes_palette("Darjeeling1", max(vertices_filter$group), type = "continuous")
mycolor <- sample(mycolor, length(mycolor))

# p<-ggraph(mygraph,layout='fr') + 
p<-ggraph(mygraph,layout='nicely') + 
  geom_edge_link(edge_colour="black", edge_alpha=0.7, edge_width=0.7) +
  geom_node_point(aes(size=n, fill=as.factor(group)), shape=21,color='black',alpha=0.9) +
  scale_size_continuous(range=c(3,9)) +
  scale_fill_manual(values=mycolor) +
  geom_node_text(aes(label=ifelse(n>2, as.character(name), "")), size=0, color="black") +
  # geom_node_label(aes(label = V(mygraph)$name), size = 2)+
  # expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2))+
  theme_minimal() +
  theme(
    legend.position="none",
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks =element_blank(),
    axis.text =element_blank(),
    axis.title = element_blank()
    #plot.margin=unit(c(0,0,0,0), "null"),
    #panel.spacing=unit(c(0,0,0,0), "null")
  )

# ggsave(filename = "plot/t1.TIFF",p, width = 5,height = 5)

# cluster <- connect %>% filter(n > 2) %>% pull(group) %>% sort %>% unique %>% max
# n_cluster <- connect$group %>% sort %>% unique %>% max
# prop <- cluster/n_cluster

cluster <- vertices_nofilter %>% filter(n > 2) %>% pull(group) %>% sort %>% unique %>% max
# n_cluster <- connect$group %>% sort %>% unique %>% max
n_cluster <- vertices_nofilter$group %>% max()
TCR_cluster <- vertices_nofilter %>% filter(n > 2) %>% nrow()

data.frame(sample = sample1,
  TCRs = nrow(t1),
  multi_cluster = cluster,
  total_cluster = n_cluster,
  TCR_cluster = TCR_cluster/nrow(t1),
  prop1 = cluster/nrow(t1), 
  prop2 =  cluster/n_cluster 
  ) %>% fwrite(paste0("cluster/500/",type,"/",sample1,"_",cutoff,".tsv"), sep = "\t")

ggsave(filename = paste0("cluster/500/",type,"/",sample1,"_",cutoff,".TIFF"),p, width = 5,height = 5)
}
}

#expanded
getplot(immdata_TRD_cutoff,"TRD_exp",0.8)
getplot(immdata_TRG_cutoff,"TRG_exp",0.8)

#Vd2+/- 
getplot(immdata_TRD_pos_cutoff,"TRD_pos",0.8)
getplot(immdata_TRD_neg_cutoff,"TRD_neg",0.8)


#Vg9+/- 
getplot(immdata_TRG_pos_cutoff,"TRG_pos",0.8)
getplot(immdata_TRG_neg_cutoff,"TRG_neg",0.8)


# #random add "slice_sample(n = 100)"
# getplot(immdata_TRD,"TRD_random",0.8)
# getplot(immdata_TRG,"TRG_random",0.8)


######### Fig 3E & FigS3D ######### 

Pos <- 
  fread("cluster/500/Merge_TRG_pos_0.8.tsv") %>%
  # fread("cluster/TRD_pos/Merge_TRD_pos_0.75.tsv") %>%
  mutate(multi_cluster=ifelse(multi_cluster == "-Inf", 0, multi_cluster)) %>%
  # mutate(prop=ifelse(prop == "-Inf", 0, prop)) %>% 
  rename_with(.fn=~ gsub("$","_Pos",.x),.cols=c(multi_cluster:prop2)) %>%
  mutate(Norm_nub_clusters_Exp = scale(multi_cluster_Pos,center = FALSE))


Neg <- 
  fread("cluster/500/Merge_TRG_neg_0.8.tsv") %>%
  # fread("cluster/TRD_neg/Merge_TRD_neg_0.75.tsv") %>%
  mutate(multi_cluster=ifelse(multi_cluster == "-Inf", 0, multi_cluster)) %>%
  # mutate(prop=ifelse(prop == "-Inf", 0, prop)) %>% 
  rename_with(.fn=~ gsub("$","_Neg",.x),.cols=c(multi_cluster:prop2)) %>%
  mutate(Norm_nub_clusters_RNM = scale(multi_cluster_Neg,center = FALSE))


library(stringr)

df1 <- Pos %>% left_join(Neg, by = c("sample")) %>%
      mutate(type = str_extract(sample, "([a-z]*[A-Z]+)")) %>%
      mutate(sample = gsub("_G.TRG", "", sample))


df2 <- reshape2::melt(df1 %>% select(sample,type,TCR_cluster_Pos,TCR_cluster_Neg),  id.vars=c("sample","type"))
# df3 <- df2 %>% filter(type == "ftc")

mycolor <- wes_palette("Darjeeling1", df2$sample %>% n_distinct, type = "continuous")
mycolor <- sample(mycolor, length(mycolor))


p<-ggplot(df2, aes(x=variable, y= 100*value))+ 
  geom_boxplot(fill = NA,width = 0.3, colour = "black")+
  geom_line(aes(group = sample),color = 'grey',size=0.5) +
  geom_point(aes(group = sample, color = sample), size=3) +
  guides(color="none")+
  scale_color_manual(values = mycolor)+
  theme_bw()+
  ylab("Normalized proportion of clustered TCRs (%)") +
  # xlab("γ-chain")+
  xlab("δ-chain")+
  stat_compare_means(method = "wilcox.test", paired=TRUE, label.x = 1.45,
  aes(label = paste0("p = ", ..p.format..))
  )+
  # annotate("text",x=1.5,y=0.42,size = 6,label ="paste(~italic(p), \" < 0.0001\")",parse = TRUE)+
  scale_x_discrete(labels = c("Vγ9+", "Vγ9-")) +
  # scale_x_discrete(labels = c("Vδ2+", "Vδ2-")) +
  # scale_x_discrete(labels = c("Expanded", "Random")) +
  theme(strip.background = element_rect(fill="grey90"),
        strip.text = element_text(size=15,face="plain",color="black"),
        axis.title=element_text(size=15,face="plain",color="black"),
        axis.text = element_text(size=15,face="plain",color="black"),
        panel.background=element_rect(colour="black",fill=NA,size = 2),
        # panel.grid=element_blank(),
        legend.background=element_rect(colour=NA,fill=NA),
        axis.ticks=element_line(colour="black"),
        legend.text=element_text(colour="black", size = 15, face = "bold"),
)
p
ggsave(filename = 'plot/t1.pdf',p,width = 5,height = 5)


df3 <- fread("cluster/500/Merge_TRG_neg_0.8.tsv") %>% 
      select(1,2,TCR_cluster) %>% 
      mutate(type = str_extract(sample, "([a-z]*[A-Z]+)")) %>%
mutate(sample = gsub("_G.TRG", "", sample))


df3$type<-factor(df3$type,levels=c("miPTC","wiPTC","PDTC","ATC"))


p<-ggboxplot(df3 %>% filter(type != "MTC"), 
          x = "type", y = "TCR_cluster",
          # fill = "type",
          alpha=0.5,
          # palette = c("#00AFBB", "#E7B800"),
          width = 0.25,
          palette = 'npg',
          add = "jitter",
          # add.params = list(size = 2, color=Type),
          bxp.errorbar = T,
          bxp.errorbar.width = 0.3)+
  stat_compare_means(method = "kruskal",label.x = 1,aes(label = paste0("p = ", ..p.format..)))+
  # stat_compare_means(
  #       comparisons = list(c("miPTC","wiPTC"),c("wiPTC","PDTC"),c("PDTC","ATC")), 
  #       aes(label = ..p.signif..),
  #       size=5,
  #         # label.y=c(0.64,0.68,0.72,0.76),
  #       method = "wilcox.test")+
  geom_violin(aes(color = type, fill = type), alpha = 0.3, width = 0.8, trim = TRUE)+
  theme_light()+
    # labs(title = "Random TCRs in δ-chain")+
  xlab(" ")+
  ylab("Proportion of clustered Vγ9- TCRs")+
  # ylab("Proportion of clustered Vδ2- TCRs")+
  # ylab("Normalized proportion of\n clustered TRD random")+
  theme(
    legend.position="none",
    plot.title = element_text(size=15,hjust = 0.5),
    panel.border = element_rect(size=2,color = "black"),
    axis.text.x=element_text(size=12,angle=0,color="Black"),
    axis.text.y=element_text(size=12,angle=0,color="Black"),
    axis.title.y = element_text(size = 15,  color = "Black"),
    axis.title.x = element_text(size = 15,  color = "Black")
        )


########### Fig 3F-3G & Fig S3E-S3F###########

expanded <- fread("Pipeline/cluster/Merge_G_expanded_0.75.tsv")

expanded$type<-factor(expanded$type,levels=c("hc","mtc","miptc","wiptc","pdtc","atc"))

my_comparisons <- list(
    c('miptc',"atc"),
    c('wiptc',"atc"),
    c('pdtc',"atc")
)


p<-ggboxplot(expanded %>% filter(type!="mtc"), 
          x = "type", y = "TCRs",
          # fill = "type",
          alpha=0.5,
          # palette = c("#00AFBB", "#E7B800"),
          width = 0.25,
          palette = 'npg',
          add = "jitter",
          # add.params = list(size = 2, color=Type),
          bxp.errorbar = T,
          bxp.errorbar.width = 0.3)+
  scale_x_discrete(labels = c("HC","miPTC","wiPTC","PDTC","ATC"))+
  #  stat_compare_means(method = "kruskal",label.x = 1,aes(label = paste0("p = ", ..p.format..)))+
  annotate("text", x = 1.5, y = 470, size = 5,label ="paste(~italic(p), \" = 0.017\")",parse = TRUE)+
  stat_compare_means(comparisons = my_comparisons, 
  aes(label = ..p.signif..),
  size=4,
    # label.y=c(580,620,664),
    label.y = c(390,420,450),
  method = "wilcox.test")+
  geom_violin(aes(color = type, fill = type), alpha = 0.3, width = 0.8, trim = TRUE)+
  theme_light()+
  labs(title = "Number of expanded clonotypes in in γ-chain")+
    # labs(title = "Number of expanded clonotypes in δ-chain")+
  xlab(" ")+
  ylab("Normalized number of clustered TCRs")+
  theme(
    legend.position="none",
    plot.title = element_text(size=15,hjust = 0.5),
    panel.border = element_rect(size=2,color = "black"),
    axis.text =element_text(size=12,angle=0,color="Black"),
    axis.title = element_text(size = 15,  color = "Black")
        )



########### FigS 3G ###########
library(ggbeeswarm)
library(ggsci)


df <- rbind(
fread("cluster/500/Merge_TRD_neg_0.8.tsv") %>% 
      select(1,2,TCR_cluster) %>% 
      mutate(type = str_extract(sample, "([a-z]*[A-Z]+)"),
          type = "Vd2-"),
fread("cluster/500/Merge_TRD_pos_0.8.tsv") %>% 
      select(1,2,TCR_cluster) %>% 
      mutate(type = str_extract(sample, "([a-z]*[A-Z]+)"),
          type = "Vd2+"),
fread("cluster/500/Merge_TRG_neg_0.8.tsv") %>% 
      select(1,2,TCR_cluster) %>% 
      mutate(type = str_extract(sample, "([a-z]*[A-Z]+)"),
          type = "Vg9-"),
fread("cluster/500/Merge_TRG_pos_0.8.tsv") %>% 
      select(1,2,TCR_cluster) %>% 
      mutate(type = str_extract(sample, "([a-z]*[A-Z]+)"),
          type = "Vg9+"))



p<-ggplot(df, aes(x=type, y=100*TCR_cluster, fill=type))+
    geom_beeswarm(aes(fill = type),shape=21,colour="black",size=4,cex=1.5)+
  geom_pointrange(stat="summary", fun.data="mean_sdl",fun.args = list(mult=1),
                  color = "black",size = 1)+
  geom_point(stat="summary", fun.y="mean",fun.args = list(mult=1),
                color = "white",size = 3)+
  stat_compare_means(method = "kruskal.test", label.x = 1,
    # label.y = 5,
    aes(label = paste0("p = ", ..p.format..)))+
  stat_compare_means(comparisons = 
    list(c("Vd2-","Vd2+"),c("Vd2-","Vg9-"),c("Vd2-","Vg9+")), 
    method = "wilcox.test",
        aes(label = ..p.signif..),
        size=6
        )+
    scale_color_npg(alpha = 1) +
    scale_fill_npg(alpha = 1) +
    ylab("Proportion of clustered TCRs(%)")+
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
    scale_y_continuous(expand = c(0, 0))


