library(immunarch)
library(ggpubr)
library(ggsci)
library(ggpubr)

file_path_TRG = "clonotype/TRG"
file_path_TRD = "clonotype/TRD"

immdata_TRG <- repLoad(file_path_TRG)
immdata_TRD <- repLoad(file_path_TRD)

##### input data #####

immdata_TRG <-repFilter(immdata_TRG, .method = "by.clonotype", .query = list(V.name = exclude("TRGV1","TRGV10", "TRGV11","TRGV5P","TRGV6","TRGV7")), .match="substring")
immdata_TRG <- repFilter(immdata_TRG, .method = "by.meta", .query = list(Sample = exclude(c("_2","_3","_4","_5","_6"))), .match="substring")

immdata_TRD <- repFilter(immdata_TRD, .method = "by.meta", .query = list(Sample = exclude(c("_2","_3","_4","_5","_6"))), .match="substring")
immdata_TRD <-repFilter(immdata_TRD, .method = "by.clonotype", .query = list(V.name = include("DV4","DV5","DV6","DV7","DV8","TRDV1","TRDV2","TRDV3")), .match="substring")
immdata_TRD <-repFilter(immdata_TRD, .method = "by.clonotype",.query =  list(J.name = exclude("AJ")), .match = "substring")



#####Fig 1B #####
TRD_len <- repExplore(immdata_TRD$data, .method = "len", .col = "aa")
TRG_len <- repExplore(immdata_TRG$data, .method = "len", .col = "aa")

vis(TRG_len, .by = "type", .meta = info)

seq(3,45, by = 3)

TRG_len <- TRG_len %>% mutate(group = case_when(
    Length >= 3 & Length < 6  ~ "4",
    Length >= 6 & Length < 9  ~ "7",
    Length >= 9 & Length < 12  ~ "10",
    Length >= 12 & Length < 15  ~ "13",
    Length >= 15 & Length < 18  ~ "16",
    Length >= 18 & Length < 21  ~ "19",
    Length >= 21 & Length < 24  ~ "22",
    Length >= 24 & Length < 27  ~ "25",
    Length >= 27 & Length < 30  ~ "28",
    Length >= 30 & Length < 33  ~ "31",
    Length >= 33 & Length < 36  ~ "34",
    Length >= 36 & Length < 39  ~ "37",
    Length >= 39 & Length < 42  ~ "40",
    Length >= 42 & Length < 45  ~ "43",
    Length >= 45 & Length < 48  ~ "46"
))

TRD_len <- TRD_len %>% mutate(group = case_when(
    Length >= 3 & Length < 6  ~ "4",
    Length >= 6 & Length < 9  ~ "7",
    Length >= 9 & Length < 12  ~ "10",
    Length >= 12 & Length < 15  ~ "13",
    Length >= 15 & Length < 18  ~ "16",
    Length >= 18 & Length < 21  ~ "19",
    Length >= 21 & Length < 24  ~ "22",
    Length >= 24 & Length < 27  ~ "25",
    Length >= 27 & Length < 30  ~ "28",
    Length >= 30 & Length < 33  ~ "31",
    Length >= 33 & Length < 36  ~ "34",
    Length >= 36 & Length < 39  ~ "37",
    Length >= 39 & Length < 42  ~ "40",
    Length >= 42 & Length < 45  ~ "43",
    Length >= 45 & Length < 48  ~ "46"
))



df <- rbind(
TRD_len %>% mutate(chain = "Delta"),
TRG_len %>% mutate(chain = "Gamma")) %>%
  mutate(group = as.numeric(group))


df_summary <- df %>%
  group_by(group,chain) %>%
  summarize(mean=mean(Count),sd = sd(Count))


df$chain <- factor(df$chain, levels = c("Gamma","Delta"))

p<-ggplot(df, aes(x=group, y=Count,fill = chain)) +
  geom_col(position = 'dodge',width = 2.5) +
  scale_fill_manual(values = c('#5D6A8C','#6E9EAC'))+
  scale_color_manual(values = c('#5D6A8C','#6E9EAC'))+
  labs(x="CDR3 length (amino acid)",y="Number of clonotyoes") +
  theme_bw()+
  theme(
    legend.title = element_blank(),
    panel.border = element_rect(size=2,color = "black"),
    # panel.grid= element_blank(),
    # plot.title = element_text(face="bold",size=25,hjust = 0.5),
    legend.text=element_text(colour="black", size = 15),
    axis.text=element_text(size=12,angle=0,color="Black", vjust=0.5),
    axis.title=element_text(size=15,angle=0,color="Black"),
    legend.position = c(0.8, 0.8)
    )+
  scale_x_continuous(breaks=c(4,7,10,13,16,19,22,25,28,31,34,37,40,43,46))
p
ggsave(filename = 'plot/t1.pdf',p, height = 4, width = 8)



new <- df_summary %>%as.data.frame %>% uncount(sum) 


p<-ggplot(new %>% filter(Length <= 46), aes(x=Length,color=chain,fill=chain))+
  geom_density(alpha=0.3,bw=1,size=1)+
  scale_color_manual(values = c('#5D6A8C','#6E9EAC'))+
  scale_fill_manual(values = c('#5D6A8C','#6E9EAC'))+
  # xlim(0,45)+
  theme_bw()+
  labs(x="CDR3 length (amino acid)",y="Density") +
  theme_bw()+
  theme(
    legend.title = element_blank(),
    panel.border = element_rect(size=2,color = "black"),
    # panel.grid= element_blank(),
    # plot.title = element_text(face="bold",size=25,hjust = 0.5),
    legend.text=element_text(colour="black", size = 15, face = "plain"),
    axis.text=element_text(face="plain",size=12,angle=0,color="Black", vjust=0.5),
    axis.title=element_text(face="plain",size=15,angle=0,color="Black"),
    legend.position = c(0.8, 0.8)
    )+
  scale_x_continuous(breaks=c(4,7,10,13,16,19,22,25,28,31,34,37,40,43,46))

p
ggsave(filename = 'plot/t2.pdf',p, height = 4, width = 8)


#####Fig 1C #####
imm_gu <- geneUsage(immdata_TRG$data, "HomoSapiens.TRGV", .norm = T)
# imm_gu <- geneUsage(immdata_TRD$data, "HomoSapiens.TRDV", .norm = T)


tmp <- imm_gu %>% t %>% as.data.frame
colnames(tmp) <- tmp[1,]
tmp <- tmp[-1,]
tmp <- tmp %>% mutate_at(vars(1:6),as.numeric) %>% 
# tmp <- tmp %>% mutate_at(vars(1:8),as.numeric) %>% 
    tibble::rownames_to_column(var = "Sample") %>% 
    reshape2::melt(id.vars="Sample") %>%
    left_join(
      fread("metadata.txt") %>% select(Type2,2), 
      by = "Sample")

tmp$Status <- factor(tmp$Status, levels = c("HC","miPTC","wiPTC","PDTC","ATC"))
tmp$Type2 <- factor(tmp$Type2, levels = c("HC","Patients"))

p <- ggplot(tmp, aes(x = variable, y = value))+
  geom_boxplot(aes(fill = Type2), outlier.size = 0,outlier.shape = NA,
      position = position_dodge(0.85), lwd = 0.9, width = 0.7)+
  geom_jitter(
    aes(fill = Type2),
        # position = position_dodge(0.85),
        position = position_jitterdodge(
            dodge.width = 0.85,
            jitter.width = 0.2),
       shape = 21,size = 4,color="black")+  
  scale_fill_npg(alpha=0.9)+
  theme_bw()+
  ylab("V chain usage")+
  xlab(" ")+
theme(  panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 1),
        axis.line.x = element_line(color = "black", size = 1),
        axis.line.y = element_line(color = "black", size = 1),
        plot.margin = unit(c(1, 0.5, 0.5, 1), "cm"),
        axis.ticks = element_line(linewidth=1, color = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        # axis.ticks.x = element_blank(),
        legend.position =  c(0.2,0.8),
        axis.text.x=element_text(face="plain",size=18,angle=0,color="Black",vjust=1),
        axis.text.y=element_text(face="plain",size=18,angle=0,color="Black"),
        axis.title.y = element_text(size = 20, face="plain", color = "Black"),
        axis.title.x = element_text(size = 20, face="plain", color = "Black"),
        panel.grid= element_blank()
)+
   expand_limits(y = c(0,0.4002))+ 
   scale_y_continuous(expand = c(0, 0))




####### public TRG (Fig S1B)#######

tc <- trackClonotypes(immdata_TRG$data, "CALWEVQELGKKIKVF", .col = "aa")


JP <- t(tc) %>% as.data.frame() %>%  
rownames_to_column() %>%
  left_join(
      fread("metadata.txt") %>% select(Type2,2), by = c("rowname"="Sample")) %>%
  mutate(V1 = as.numeric(V1))



JP$Status <- factor(JP$Status, levels = c("HC","MTC","miPTC","wiPTC","PDTC","ATC"))
TT$chain <- factor(TT$chain, levels = c("Gamma","Delta"))

p <- ggplot(JP%>% filter(V1<0.05), aes(x = Type2, y = V1*100, fill = factor(Type2)))+
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 1,
      position = position_dodge(1))+
  geom_pointrange(stat="summary", fun.data="mean_sdl",fun.args = list(mult=1),
                  color = "black",size = 1.2, 
                  position = position_dodge(1))+
  geom_point(stat="summary", fun.y="mean",fun.args = list(mult=1),
                color = "white",size = 3.5, 
                position = position_dodge(1))+
  scale_fill_nejm(alpha = 0.86)+
    theme_bw()+
  ylab("public Vγ9-JγP clonotype(%)")+
  xlab("")+
theme(  panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 1),
        axis.line.x = element_line(color = "black", size = 1),
        axis.line.y = element_line(color = "black", size = 1),
        plot.margin = unit(c(1, 0.5, 0.5, 1), "cm"),
        axis.ticks = element_line(linewidth=1, color = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        # axis.ticks.x = element_blank(),
        # legend.position =  c(0.2,0.8),
        axis.text.x=element_text(face="plain",size=12,angle=0,color="Black",vjust=1),
        axis.text.y=element_text(face="plain",size=12,angle=0,color="Black"),
        axis.title.y = element_text(size = 12, face="plain", color = "Black"),
        axis.title.x = element_text(size = 12, face="plain", color = "Black"),
        panel.grid= element_blank()
)+
  stat_compare_means(method = "wilcox.test",aes(label = paste0("p = ", ..p.format..)))+
  scale_y_continuous(expand = c(0, 0)) 




############ clonal space homeostasis (Fig 1E and Fig S1C)############

immdata_TRG_pos <- repFilter(immdata_TRG, .method = "by.clonotype", .query = list(V.name = include("TRGV9")), .match="substring")
immdata_TRG_neg <- repFilter(immdata_TRG, .method = "by.clonotype", .query = list(V.name = exclude("TRGV9")), .match="substring")

immdata_TRD_pos <- repFilter(immdata_TRD, .method = "by.clonotype", .query = list(V.name = include("TRDV2")), .match="substring")
immdata_TRD_neg <- repFilter(immdata_TRD, .method = "by.clonotype", .query = list(V.name = exclude("TRDV2")), .match="substring")



imm_hom <- repClonality(immdata_TRG_pos$data,
  .method = "homeo",
  .clone.types = c(Rare = .00001, Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1)
)


d_hom <- imm_hom %>% melt %>%
  left_join(
      fread("metadata.txt") %>% select(Type2,2), 
      by = c("Var1"="Sample")) %>%
  # mutate(Var1 = gsub("_G\\.TRG", "", Var1))
  mutate(Var1 = gsub("_D\\.TRD", "", Var1))


p2 <- ggplot(d_hom 
# %>% filter(variable == "Hyperexpanded")
, aes(x=Var1,y=value,fill=Var2))+
  geom_bar(stat="identity",width=0.9,
  # position='fill'
  position = "stack"
  )+
  theme_bw()+
  scale_fill_manual(values = c("#B5D8E0","#65C5DC","#0095c7","#0877b7","#1a548b"))+
  ylab("Relative abundance of V9+")+
  xlab("")+
  theme(
    legend.title =element_blank(),
    legend.text = element_text(colour="black", size = 8, face = "plain"),
    axis.title.y=element_text(face="plain",size=10,angle=90,color="Black", vjust=0.5),
    axis.text.y=element_text(face="plain",size=8,angle=0,color="Black"),
    axis.text.x=element_text(size=7,angle=90,color="Black",vjust=0.5,hjust = 1),
    plot.title = element_text(face="bold",size=13,hjust = 0.5),
    panel.grid= element_blank(),
    panel.background = element_rect(colour = "black", size=0.2, fill=NA),
    strip.text.x = element_text(size = 10, face = "bold"),
    strip.background.x = element_rect(fill = "#DFEDFA"),
    legend.position = "bottom"
  )+
  # scale_y_continuous(labels=c("0","25","50","75","100"))+
  facet_grid(. ~ Type2, scales="free_x", space = "free")



d3 <- imm_hom %>% melt %>%
  left_join(
      fread("metadata.txt") %>% select(Type2,2), 
      by = c("Var1"="Sample")) %>%
  mutate(Var1 = gsub("_D\\.TRD", "", Var1))

d3$Type2 <- factor(d3$Type2, levels = c("HC","Patients"))
d3$Status <- factor(d3$Status, levels = c("HC","MTC","miPTC","wiPTC","PDTC","ATC"))

p<-ggplot(d3, aes(x=Var2, y =value))+
  geom_boxplot(aes(fill = Type2), outlier.size = 0,outlier.shape = NA,
        position = position_dodge(0.8), lwd = 0.9, width = 0.6)+
  geom_jitter(
        aes(fill = Type2),
        position = position_jitterdodge(
        dodge.width = 0.8,
        jitter.width = 0.3),
        shape = 21,size = 3,color="black")+
    scale_fill_nejm(alpha = 0.65)+
    ylab("") +
    xlab(" ")+
    theme_bw()+
    theme(
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 1),
        axis.line.x = element_line(color = "black", size = 1),
        axis.line.y = element_line(color = "black", size = 1),
        plot.margin = unit(c(1, 0.5, 0.5, 1), "cm"),
        axis.ticks = element_line(linewidth=1, color = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        # axis.ticks.x = element_blank(),
        axis.text.x=element_text(face="plain",size=12,angle=0,color="Black",vjust=1),
        axis.text.y=element_text(face="plain",size=12,angle=0,color="Black"),
        axis.title.y = element_text(size = 10, face="plain", color = "Black"),
        axis.title.x = element_text(size = 10, face="plain", color = "Black"),
        panel.grid= element_blank(),
        legend.position = c(0.2,0.7)
    )+
    stat_compare_means(
        aes(group=Type2, label = ..p.signif..), 
        label.y = 0.83,
        size=5,
        # aes(group=CLASS,label = paste0("p = ", ..p.format..)),
        method = "wilcox.test")+
    expand_limits(y = c(0,1))+ 
# scale_x_continuous(expand = c(0, 0)) + 
scale_y_continuous(expand = c(0, 0)) 

