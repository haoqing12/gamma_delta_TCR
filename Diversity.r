

library(ggpubr)
library(DescTools)
library(vegan)

getdiv <- function(immdata,v) {
TRG_d <- data.frame()
for (n in immdata$meta$Sample){
      shannon_entropy <- vegan::diversity(immdata$data[[n]]$Clones , "shannon")
      Pielou_index <- shannon_entropy / log(vegan::specnumber(immdata$data[[n]]$Clones))
      Clonality = 1 - Pielou_index
      Gini.Simpson <- vegan::diversity(immdata$data[[n]]$Clones, "simpson")
      inverse.Simpson <- vegan::diversity(immdata$data[[n]]$Clones, "invsimpson")
      Gini_coefficient <- Gini(immdata$data[[n]]$Clones)

      diversity <- data.frame(sample = n,
                              type = v,
                              unique.clonotypes = nrow(immdata$data[[n]]),
                              total.clonotypes = sum(immdata$data[[n]]$Clones),
                              shannon_entropy = shannon_entropy,
                              Gini_Simpson = Gini.Simpson,
                              inv_Simpson = inverse.Simpson,
                              Clonality = Clonality,
                              Gini_coefficient = Gini_coefficient
      )
     TRG_d <- rbind(TRG_d, diversity)
}
  return(TRG_d)
}

div_trd <- getdiv(immdata_TRD,"Delta")%>% mutate(Type2 = str_extract(sample, "([a-z]*[A-Z]+)"))
div_trg <- getdiv(immdata_TRG,"Gamma")%>% mutate(Type2 = str_extract(sample, "([a-z]*[A-Z]+)"))


### Fig S1A ###
df2 <- cbind(div_trg %>% dplyr::rename_with(.fn=~ gsub("$","_G",.x),.cols=c(sample:Gini_coefficient)),
            div_trd %>% dplyr::rename_with(.fn=~ gsub("$","_D",.x),.cols=c(sample:Gini_coefficient)))


newdf <- melt(df2%>% select(total.clonotypes_G, total.clonotypes_D) %>% rename("Gamma"=1, "Delta"=2))


p<-ggplot(newdf, aes(x=variable, y=value, fill=variable))+
  geom_beeswarm(aes(fill = variable),shape=21,colour="black",cex=4,size=6.5)+
  geom_pointrange(stat="summary", fun.data="mean_sdl",fun.args = list(mult=1),
        color = "black",size = 1.4)+
  geom_point(stat="summary", fun.y="mean",fun.args = list(mult=1),color = "white",size = 5)+
  scale_fill_nejm(alpha = 0.9)+
  stat_compare_means(method = "wilcox.test", 
    aes(label = paste0("p = ", ..p.format..)))+
  theme_bw()+
  ylab("Total number of TCRs") +
  xlab(" ") +
  # scale_y_continuous(labels = c(0,"2x10^6","4x10^6","6x10^6","8x10^6"))+
  theme( panel.border = element_blank(),
    axis.line.x = element_line(color = "black", size = 1),
    axis.line.y = element_line(color = "black", size = 1),
    plot.margin = unit(c(1, 0.5, 0.5, 1), "cm"),
    axis.ticks = element_line(linewidth=1, color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    # axis.ticks.x = element_blank(),
    axis.text.x=element_text(face="plain",size=15,angle=0,color="Black",vjust=1),
    axis.text.y=element_text(face="plain",size=15,angle=0,color="Black"),
    axis.title.y = element_text(size = 18, face="plain", color = "Black"),
    axis.title.x = element_text(size = 18, face="plain", color = "Black"),
    panel.grid= element_blank(),
    legend.position = 'none'
  )



p<-ggscatter(df2,  x = "unique.clonotypes_G", y = "unique.clonotypes_D",
          add = "reg.line", 
          size=2, # Add regressin line
          add.params = list(color = "grey30", fill = "white") # Customize reg. line
          )+
  scale_color_npg(alpha = 0.7)+
  stat_cor(method = "spearman",cor.coef.name =  "rho")+
  ylab("Delta")+ xlab("Gamma")+
  theme_bw()+
  labs(title = "Correlation of γδTCRs clonotypes")+
  theme(
    # panel.grid= element_blank(),
    panel.border = element_rect(size=2),
    plot.title = element_text(face="bold",size=25,hjust = 0.5),
    axis.text =element_text(size=15,angle=0,color="Black"),
    axis.title =element_text(face="bold",size=20,angle=0,color="Black"),
  )




#### Fig 1F (HC vs TC)####

# div_trd_neg <- getdiv(immdata_TRD_neg,"neg")
# div_trd_pos <- getdiv(immdata_TRD_pos,"pos")

div_trg_neg <- getdiv(immdata_TRG_neg,"neg")
div_trg_pos <- getdiv(immdata_TRG_pos,"pos")


tmp <- div_trg_pos %>% left_join(div_trg_neg, by = "sample") %>% 
        left_join(fread("meta.data") %>% mutate(sample = paste0(Patient,"_G.TRG")), by = "sample")

tmp_ <- tmp %>% 
        select(1,CLASS,shannon_entropy.x,shannon_entropy.y) %>% 
        rename("Vδ2+"=3,"Vδ2-"=4) %>% melt()
        # rename("Vγ9+"=3,"Vγ9-"=4) %>% melt() 


p<-ggplot(tmp_, aes(x=variable, y =value))+
    geom_violin(aes(fill=factor(CLASS)),color = "white", position = position_dodge(0.8),alpha = 0.3, width = 1, trim = TRUE)+
    geom_boxplot(outlier.size = 1, width = 0.3, aes(color=factor(CLASS)),
               position = position_dodge(0.8),size=1)+
    scale_fill_manual(values = c('#5D6A8C','#6E9EAC'))+
    scale_color_manual(values = c('#5D6A8C','#6E9EAC'))+
    # scale_color_manual(values = c('#333462d3','#97423acd'))+
    # ylab("Shannon entropy of γδTCR repertoires") +
    ylab("Gini coeff of δTCR repertoires") +
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
axis.title.y = element_text(size = 15, face="plain", color = "Black"),
axis.title.x = element_text(size = 15, face="plain", color = "Black"),
panel.grid= element_blank()
    )+
    # ylim(2.5,7.8)+
    stat_compare_means(
        aes(group=CLASS
        , label = ..p.signif..
        ), 
        # label.y = 37,
        label.y = 0.95,
        size=5,
        # aes(group=CLASS,label = paste0("p = ", ..p.format..)),
        method = "wilcox.test")+
  # expand_limits(y = c(0,40))+ 
# scale_x_continuous(expand = c(0, 0)) + 
scale_y_continuous(expand = c(0, 0)) 



### Fig2C and FigS2

p<-ggboxplot(div_trg_pos %>% filter(Type != "MTC"),
          x = "Type", y = "Clonality",fill = "Type",
          alpha=0.5,
          width = 0.25,
          palette = 'npg',
          add = "jitter",
          bxp.errorbar = T,
          bxp.errorbar.width = 0.3)+
  stat_compare_means(method = "kruskal.test",
    aes(label = paste0("p = ", ..p.format..)))+
  # stat_compare_means(method = "wilcox.test",aes(label = ..p.signif..),
  #   comparisons = list(
  #     c("miPTC","HC"),c("wiPTC","PDTC"),
  #     c("HC","ATC"),
  #     c("miPTC","ATC"),
  #     c("wiPTC","ATC"),
  #   c("PDTC","ATC")
  #   )
  #   )+
  geom_violin(aes(color = Type, fill = Type), alpha = 0.3, width = 1, trim = TRUE)+
    theme_light()+
  xlab("")+
  # ylab("Shannon entropy of Vδ2+ repertoires")+
  # ylab("Shannon entropy of Vγ9+ repertoires")+
  # ylab("Clonality of Vδ2+ repertoires")+
  ylab("Clonality of Vγ9+ repertoires")+
  # ylab("Gini coefficient of Vδ2+ repertoires")+
  # ylab("Gini coefficient of Vγ9+ repertoires")+
  # ylab("D50 of Vδ2+ repertoires")+
  #  ylab("D50 of Vγ9+ repertoires")+
  theme(
    legend.position="none",
    panel.border = element_rect(size=1.5,color = "black"),
    axis.text.x=element_text(size=12,angle=0,color="Black"),
    axis.text.y=element_text(size=12,angle=0,color="Black"),
    axis.title.y = element_text(size = 12,  color = "Black"),
    axis.title.x = element_text(size = 12,  color = "Black")
        )
