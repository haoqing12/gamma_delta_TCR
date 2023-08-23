

########## Expanded TCR  ########## 
####### Fig 3A ####### 

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


getinfo <- function(imm,Chain,imm_all) {
  DF <- data.frame()
  for (i in imm$meta$Sample){
      df <- imm$data[[i]]
      temp <- data.frame(
        sample = i,
        chain = Chain,
        n_clonotype = nrow(df),
        prop_clonotype =  nrow(df)/nrow(imm_all$data[[i]]),
        prop_reads = (df$Clones %>% sum()) / (imm_all$data[[i]]$Clones%>% sum()),
        all_reads = df$Proportion %>% sum()
      )
      DF <- rbind(temp, DF)
  }
  return(DF)
}

df_50cutoff_info <- 
    rbind(
      rbind(getinfo(immdata_TRD_neg_cutoff,"TRD_neg",immdata_TRD_neg),
      getinfo(immdata_TRD_pos_cutoff,"TRD_pos",immdata_TRD_pos))%>% 
        left_join(
          fread("metadata.txt") %>% select(Status,2), 
          by = c("sample"="Sample")
          ),

      rbind(getinfo(immdata_TRG_neg_cutoff,"TRG_neg",immdata_TRG_neg),
      getinfo(immdata_TRG_pos_cutoff,"TRG_pos",immdata_TRG_pos))%>% 
          left_join(
            fread("metadata.txt") %>% select(Status,2), 
            by = c("sample"="Sample"))
    ) %>% mutate(chain_ = ifelse(grepl("TRG",chain), "TRG","TRD"))



p<-ggboxplot(df_50cutoff_info,
         x = "chain_",
          # y = "prop_clonotype",
          y =  "prop_reads",
          color = "chain_",
          # palette = c('#333462d3','#97423acd'),
          width = 0.5,
          add = "jitter",
          add.params = list(size = 4),
          bxp.errorbar = T,
          bxp.errorbar.width = 0.2)+
scale_color_nejm()+
  xlab("")+
  ylab('Percentage of total producitve reads represented by TCRexp')+
  # ylab('Percentage of total unique clonotypes\nrepresented by TCRexp')+
  theme_bw()+
  theme( panel.border = element_blank(),
      axis.line = element_line(color = "black", size = 1),
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
      legend.position = "none"
)+
# coord_cartesian(ylim = c(0, 0.8))+
# scale_y_continuous(breaks = seq(0, 0.8, by = 0.2))+
expand_limits(y = c(0.8,1))+ 
# scale_x_continuous(expand = c(0, 0)) + 
scale_y_continuous(expand = c(0, 0)) +
stat_compare_means(
    # aes(group=Type2, label = ..p.signif..), 
    # label.y = 0.63,
    size=5,
    aes(label = paste0("p = ", ..p.format..)),
    method = "wilcox.test")

