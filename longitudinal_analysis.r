library(immunarch)
library(ggalluvial)
library(tidyr)
library(ggsci)
library(ggpubr)


immdata_TRD <- repLoad(file_path_D)

immdata_TRD <-repFilter(immdata_TRD, .method = "by.clonotype", .query = list(V.name = include("DV4","DV5","DV6","DV7","DV8","TRDV1","TRDV2","TRDV3")), .match="substring")
immdata_TRD <-repFilter(immdata_TRD, .method = "by.clonotype",.query =  list(J.name = exclude("AJ")), .match = "substring")


immdata_TRG <- repLoad(file_path_G)
immdata_TRG <-repFilter(immdata_TRG, .method = "by.clonotype", .query = list(V.name = exclude("TRGV1","TRGV10", "TRGV11","TRGV5P","TRGV6","TRGV7")), .match="substring")




immdata_TRD_pos <- repFilter(immdata_TRD, .method = "by.clonotype", .query = list(V.name = include("TRDV2")), .match="substring")
immdata_TRD_neg <- repFilter(immdata_TRD, .method = "by.clonotype", .query = list(V.name = exclude("TRDV2")), .match="substring")


immdata_TRG_pos <- repFilter(immdata_TRG, .method = "by.clonotype", .query = list(V.name = include("TRGV9")), .match="substring")
immdata_TRG_neg <- repFilter(immdata_TRG, .method = "by.clonotype", .query = list(V.name = exclude("TRGV9")), .match="substring")


getcutoff <- function(immdata_TRD){
    immdata_TRD_cutoff <- list()
    for (i in immdata_TRD$meta$Sample){
      print(i)
      cutoff_ <- immdata_TRD$data[[i]]$Clones %>% sum/nrow(immdata_TRD$data[[i]])
      immdata_TRD_cutoff$data[[i]] <- immdata_TRD$data[[i]] %>% filter(Clones >= cutoff_/2) %>% arrange(desc(Clones))
    }
    immdata_TRD_cutoff$meta <- immdata_TRD$meta

  return(immdata_TRD_cutoff)
}


immdata_TRD_neg_cutoff <- getcutoff(immdata_TRD_neg)
immdata_TRD_pos_cutoff <- getcutoff(immdata_TRD_pos)

immdata_TRG_neg_cutoff <- getcutoff(immdata_TRG_neg)
immdata_TRG_pos_cutoff <- getcutoff(immdata_TRG_pos)




get_all <- function(immdata){
    
    t1 =  immdata$data[[1]] %>% mutate(Proportion = Proportion*100) %>% rename(Clones_t1 = 1, Proportion_t1 = 2) %>% select(1,2,4) %>% group_by(CDR3.aa) %>% summarise(Clones_t1 = sum(Clones_t1),Proportion_t1= sum(Proportion_t1))
    t2 =  immdata$data[[2]] %>% mutate(Proportion = Proportion*100) %>% rename(Clones_t2 = 1, Proportion_t2 = 2) %>% select(1,2,4) %>% group_by(CDR3.aa) %>% summarise(Clones_t2 = sum(Clones_t2),Proportion_t2= sum(Proportion_t2))
    t3 =  immdata$data[[3]] %>% mutate(Proportion = Proportion*100) %>% rename(Clones_t3 = 1, Proportion_t3 = 2) %>% select(1,2,4) %>% group_by(CDR3.aa) %>% summarise(Clones_t3 = sum(Clones_t3),Proportion_t3= sum(Proportion_t3))
    t4 =  immdata$data[[4]] %>% mutate(Proportion = Proportion*100) %>% rename(Clones_t4 = 1, Proportion_t4 = 2) %>% select(1,2,4) %>% group_by(CDR3.aa) %>% summarise(Clones_t4 = sum(Clones_t4),Proportion_t4= sum(Proportion_t4))
    t6 =  immdata$data[[5]] %>% mutate(Proportion = Proportion*100) %>% rename(Clones_t6 = 1, Proportion_t6 = 2) %>% select(1,2,4) %>% group_by(CDR3.aa) %>% summarise(Clones_t6 = sum(Clones_t6),Proportion_t6= sum(Proportion_t6))
    t7 =  immdata$data[[6]] %>% mutate(Proportion = Proportion*100) %>% rename(Clones_t7 = 1, Proportion_t7 = 2) %>% select(1,2,4) %>% group_by(CDR3.aa) %>% summarise(Clones_t7 = sum(Clones_t7),Proportion_t7= sum(Proportion_t7))

PAT_all <- t1 %>% full_join(t2, by = "CDR3.aa")%>% full_join(t3, by = "CDR3.aa")%>% full_join(t4, by = "CDR3.aa")%>% 
full_join(t5, by = "CDR3.aa")%>% full_join(t6, by = "CDR3.aa") %>% 
mutate(
  type = case_when(
    Proportion_t1 >0 ~ "Pre-existed",
    is.na(Proportion_t1) & (Proportion_t2>0 & Proportion_t3>0 & Proportion_t4>0) ~ "TKI",
    is.na(Proportion_t1) & ((Proportion_t2>0 & Proportion_t3>0)|(Proportion_t2>0 & Proportion_t4>0 ) | (Proportion_t3>0 & Proportion_t4>0)) ~ "TKI",
    is.na(Proportion_t1) & ( Proportion_t2>0 | Proportion_t3>0| Proportion_t4>0 ) ~ "TKI_indi",
    is.na(Proportion_t1) & is.na(Proportion_t2) & is.na(Proportion_t3) & is.na(Proportion_t4) & (Proportion_t5>0)  ~ "Radio",
    is.na(Proportion_t1) & is.na(Proportion_t2) & is.na(Proportion_t3) & is.na(Proportion_t4)  & is.na(Proportion_t5) & Proportion_t6>0  ~ "PD-1"
    )
    )
    return(PAT_all)
}


########### Fig 4B ##############
df <- rbind(
    get_all(immdata_TRD_neg_cutoff) %>% filter(type %in% c("TKI","Radio","PD-1")) %>% mutate(class = "Vd2-"),
    get_all(immdata_TRD_pos_cutoff) %>% filter(type %in% c("TKI","Radio","PD-1")) %>% mutate(class = "Vd2+"),
    get_all(immdata_TRG_neg_cutoff) %>% filter(type %in% c("TKI","Radio","PD-1")) %>%mutate(class = "Vg9-"),
    get_all(immdata_TRG_pos_cutoff) %>% filter(type %in% c("TKI","Radio","PD-1")) %>%mutate(class = "Vg9+")
    )%>%
    select(CDR3.aa,Proportion_t1,Proportion_t2,Proportion_t3,Proportion_t4,Proportion_t5,Proportion_t6,class)

df[is.na(df)] <- 0

df2 <- df %>% melt() %>% mutate(value_ = log(value))
df2$class <- factor(df2$class, levels = c("Vd2-","Vd2+","Vg9-","Vg9+"))

p<-ggplot(df2
#  %>% filter(value_ != "-Inf") 
 ,aes(x=variable,y=value,fill=factor(class)))+

  stat_summary(aes(color=class, group = class, linetype = class),geom = 'line',fun='mean',cex=1,linewidth = 2)+
  stat_summary(aes(color=class, group = class, shape = class), fun.y = "mean", geom = "point",size=8) +        # adding data points
  stat_summary(aes(color=class),geom = 'errorbar',cex=1,width=0.2,linewidth = 1)+
  scale_shape_manual(values = c(15,16,17,18))+
  scale_linetype_manual(values = c( "solid", "dashed","twodash","dotted"))+
#   scale_color_manual(values = c('#0886CC','#DA4E33'))+
  theme_classic(base_size = 15)+
  scale_color_nejm()+
  ylab("Proportion of Neo TCRexp (%)")+
  xlab("")+
  scale_x_discrete(labels = c(
    "T1",
   "T2", "T3","T4","T5","T6"))+
  theme( panel.border = element_blank(),
  legend.position = c(0.2,0.8),
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
# expand_limits(y = c(0,0.035))+ 
scale_y_continuous(
  # limits = c(-6,-3),
  expand = c(0, 0))+
    annotate(geom = 'segment',x=6.5,xend = 6.5,y = 0.012,yend = 0.032,cex=1)+
  annotate(geom = 'text',label='****',x=6.7,y=0.02,size=8,angle=90)





########## Fig 4A ########## 

library(ggpubr)
library(DescTools)
library(vegan)

getdiv <-function(S1){
info <- data.frame()
for (i in c(1,2,3,4,6,7)){
    temp <- PAT_all  %>% filter(!is.na(paste0("Clones_t",i))) %>%
          mutate( 
            class = case_when(
              type == "Pre-existed" ~ "Pre-existed",
              # type == "TKI_indi" ~ "Neo",
              type == "TKI" ~ "Neo",
              type == "Radio" ~ "Neo",
              type == "PD-1" ~ "Neo")) %>% 
          filter(class == S1) %>%
          select(CDR3.aa,paste0("Clones_t",i),class,paste0("Proportion_t",i)) %>% drop_na() %>%
          rename("Clones"=2, "Proportion"=4)

        shannon.entropy <- vegan::diversity(temp$Clones, "shannon")
        Gini.Simpson <- vegan::diversity(temp$Clones, "simpson")
        inverse.Simpson <- vegan::diversity(temp$Clones, "invsimpson")
        Pielou_index <- shannon.entropy / log(vegan::specnumber(temp$Clones))
        Gini_coefficient <- Gini(temp$Clones)
        #Clonality
        Clonality = 1 - Pielou_index

    getcounts <- data.frame(
        sample = i,
        type_ = S1,
        number_of_clonotype = nrow(temp),
        Proportion_of_clones = temp %>% pull(Proportion) %>% sum(),
        shannon_entropy = shannon.entropy,
        Gini.Simpson = Gini.Simpson,
        inverse.Simpson = inverse.Simpson,
        Gini_coefficient = Gini_coefficient,
        Clonality = Clonality*7.5)

  info <- rbind(info,getcounts)
}
  return(info)
}

df_PAT <- rbind(
        getdiv("Pre-existed"),
        getdiv("Neo")) %>%
        mutate(TT = paste0(sample,"_",type_))

df_PAT[7,"shannon_entropy"] <- 0
df_PAT[7,"Clonality"] <- 0


df1 <- df_PAT %>% 
  mutate(sample = as.character(sample)) %>%
  select(sample,shannon_entropy,Clonality,type_)  %>% melt(id.vars=c("sample","type_")) %>%
  mutate(TT = paste0(type_,"_",variable))

p <- ggplot(df1)+
  geom_line(aes(x=sample, y = value,group = TT, color=TT),size=2)+
  geom_point(aes(x=sample, y = value,group = TT, color=TT),size=3)+
  # scale_color_npg()+
  scale_color_manual(values = rev(c("#BC3C29B2","#0072B5B2","#BC3C29","#0072B5")))+
  ylab("Shannon entropy")+
  scale_x_discrete(labels = c("T1","T2","T3","T4","T5","T6"))+
  scale_y_continuous(sec.axis = sec_axis(~./7.5, name = "Clonality"))+
  theme_minimal()+
  xlab("\nTime") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.3,0.8),
    # legend.position = "none",
    axis.ticks = element_line(linewidth=1, color = "black"),
    axis.text.x=element_text(size=15,angle=0,color="Black"),
    axis.text.y=element_text(size=15,angle=0,color="Black"),
    axis.title.y = element_text(size = 15,face="bold",color = "Black"),
    axis.title.x = element_text(size = 15, face="bold",color = "Black"),
    legend.text  = element_text(size=15,face="plain",color="black"),
    legend.title = element_blank(),
    # panel.grid= element_blank()
    )



######## FigS4A ########

PAT_all <- get_all(immdata_TRD_neg_cutoff)

PAT_ <- t1 %>% full_join(t2, by = "CDR3.aa")%>% full_join(t3, by = "CDR3.aa")%>% full_join(t4, by = "CDR3.aa")%>% 
      full_join(t5, by = "CDR3.aa")%>% full_join(t6, by = "CDR3.aa") %>%
      left_join(PAT_all %>% select(1,type))

PAT_[is.na(PAT_)] <- 0

PAT_diff <- PAT_ %>% arrange(desc(Proportion_t6)) %>% head(50) %>%
    mutate(diff = Proportion_t6 - Proportion_t1) %>% 
    mutate(
    class2 = case_when(type == "Pre-existed" ~ "Pre",
                      type %in% c("TKI","PD-1","Radio","TKI_indi") ~ "NEO"),
    class = case_when(
            diff == 0 & class2 != "NEO" ~ "stable",
            diff < 1 & class2 != "NEO" & diff > -1  ~ "stable",
            diff > 1  & class2 != "NEO"~ "Expanding",
            diff < 1  & class2 != "NEO"~ "Contracting",
            class2 == "NEO" ~ "NEO"
          )
    ) %>% 
    as.data.frame()

df_t <- PAT_diff %>%  select(1,Proportion_t1,Proportion_t6,diff,class)

df_m <- PAT_diff %>%select(1,Proportion_t1,Proportion_t2,Proportion_t3,Proportion_t4,Proportion_t5,Proportion_t6,class)%>% 
    rename("1"=2,"2"=3,"3"=4,"4"=5,"5"=6,"6"=7) %>%
    melt(vars=c(CDR3.aa,class))

df_m$CDR3.aa <- factor(df_m$CDR3.aa, levels =
      c(df_t %>% filter(class == "NEO") %>% pull(CDR3.aa),
      df_t %>% filter(class == "Expanding") %>% pull(CDR3.aa),
      df_t %>% filter(class %in% c("stable","Contracting")) %>% pull(CDR3.aa))
)

library(RColorBrewer)

display.brewer.pal(5,"Greys")

mycolor1 <- colorRampPalette(brewer.pal(5,"YlOrRd"))(df_t %>% filter(class=="NEO") %>% nrow())
mycolor2 <- colorRampPalette(brewer.pal(5,"Blues"))(df_t %>% filter(class=="Expanding") %>% nrow())
mycolor3 <- colorRampPalette(brewer.pal(5,"Greys"))(df_t %>% filter(class %in% c("stable","Contracting")) %>% nrow())



p<-ggplot(df_m, aes(x = variable, y = value,  group = CDR3.aa, fill = CDR3.aa))+
  geom_area(position = "stack")+
  geom_line(position = "stack", color = NA)+
  scale_fill_manual(values = c( mycolor1,mycolor2,rev(mycolor3)))+
  theme_bw()+
  expand_limits(y = c(0,75))+
  ylab("Productive reads of Vd2+(%)")+
  xlab(" ")+
  # scale_x_discrete(label = c("V"))
  theme( panel.border = element_blank(),
  axis.line.x = element_line(color = "black", size = 0.8),
  axis.line.y = element_line(color = "black", size = 0.8),
  plot.margin = unit(c(1, 0.5, 0.5, 1), "cm"),
  axis.ticks = element_line(linewidth=0.8, color = "black"),
  axis.ticks.length = unit(0.25, "cm"),
  legend.position = "none",
  # axis.ticks.x = element_blank(),
  axis.text.x=element_text(face="plain",size=15,angle=0,color="Black",vjust=1),
  axis.text.y=element_text(face="plain",size=15,angle=0,color="Black"),
  axis.title.y = element_text(size = 15, face="plain", color = "Black"),
  axis.title.x = element_text(size = 15, face="plain", color = "Black"),
  panel.grid= element_blank()
  )+
  scale_y_continuous(expand = c(0, 0), breaks= c(0,25,50,75))




########## Fig 4C ########## 
category <- c("Pre-existed","TA","RA","PA")

PAT_all <- get_all(immdata_TRD_neg_cutoff)

T1_porp <-  c(PAT_all %>% filter(!is.na(Proportion_t1) & type == "Pre-existed") %>% pull(Proportion_t1) %>% sum()/PAT_all %>% filter(!is.na(Proportion_t1)) %>% pull(Proportion_t1) %>% sum(),
              0,0,0)

T2_porp <- c(PAT_all %>% filter(!is.na(Proportion_t2) & type == "Pre-existed") %>% pull(Proportion_t2) %>% sum()/PAT_all %>% filter(!is.na(Proportion_t2) )%>% pull(Proportion_t2) %>% sum(),
             PAT_all %>% filter(!is.na(Proportion_t2) & type %in% c("TKI","TKI_indi")) %>% pull(Proportion_t2) %>% sum()/PAT_all %>% filter(!is.na(Proportion_t2)) %>% pull(Proportion_t2) %>% sum(),
             0,0
)

T3_porp <- c(PAT_all %>% filter(!is.na(Proportion_t3) & type == "Pre-existed") %>% pull(Proportion_t3) %>% sum()/PAT_all %>% filter(!is.na(Proportion_t3)) %>% pull(Proportion_t3) %>% sum(),
             PAT_all %>% filter(!is.na(Proportion_t3) & type  %in% c("TKI","TKI_indi")) %>% pull(Proportion_t3) %>% sum()/PAT_all %>% filter(!is.na(Proportion_t3)) %>% pull(Proportion_t3) %>% sum(),
             0,0
)

T4_porp <- c(PAT_all %>% filter(!is.na(Proportion_t4) & type == "Pre-existed") %>% pull(Proportion_t4) %>% sum()/PAT_all %>% filter(!is.na(Proportion_t4)) %>% pull(Proportion_t4) %>% sum(),
             PAT_all %>% filter(!is.na(Proportion_t4) & type  %in% c("TKI","TKI_indi")) %>% pull(Proportion_t4) %>% sum()/PAT_all %>% filter(!is.na(Proportion_t4)) %>% pull(Proportion_t4) %>% sum(),
             0,0
)

T5_porp <- c(PAT_all %>% filter(!is.na(Proportion_t5) & type == "Pre-existed") %>% pull(Proportion_t6) %>% sum()/PAT_all %>% filter(!is.na(Proportion_t5)) %>% pull(Proportion_t5) %>% sum(),
             PAT_all %>% filter(!is.na(Proportion_t5) &  type  %in% c("TKI","TKI_indi")) %>% pull(Proportion_t5) %>% sum()/PAT_all %>% filter(!is.na(Proportion_t5)) %>% pull(Proportion_t5) %>% sum(),
             PAT_all %>% filter(!is.na(Proportion_t5) &  type == "Radio") %>% pull(Proportion_t5) %>% sum()/PAT_all %>% filter(!is.na(Proportion_t5)) %>% pull(Proportion_t5) %>% sum(),
             0
)

T6_porp <- c(PAT_all %>% filter(!is.na(Proportion_t6) & type == "Pre-existed") %>% pull(Proportion_t6) %>% sum()/PAT_all %>% filter(!is.na(Proportion_t6)) %>% pull(Proportion_t6) %>% sum(),
             PAT_all %>% filter(!is.na(Proportion_t6) &  type %in% c("TKI","TKI_indi")) %>% pull(Proportion_t6) %>% sum()/PAT_all %>% filter(!is.na(Proportion_t6)) %>% pull(Proportion_t6) %>% sum(),
             PAT_all %>% filter(!is.na(Proportion_t6) &  type == "Radio") %>% pull(Proportion_t6) %>% sum()/PAT_all %>% filter(!is.na(Proportion_t6)) %>% pull(Proportion_t6) %>% sum(),
             PAT_all %>% filter(!is.na(Proportion_t6) &  type == "PD-1") %>% pull(Proportion_t6) %>% sum()/PAT_all %>% filter(!is.na(Proportion_t6)) %>% pull(Proportion_t6) %>% sum()
)


TXQ_flow <- cbind(category,T1_porp,T2_porp,T3_porp,T4_porp,T5_porp,T6_porp) %>% as.data.frame() %>%
  melt(id.vars = 'category') %>%
  mutate(value = as.numeric(value))

TXQ_flow$category<-factor(TXQ_flow$category,levels=category)

library(ggalluvial)

p<-ggplot(TXQ_flow %>% ungroup(), 
    aes(x=variable, y = value, fill = category,label = round(100*value,2),
    stratum= category, alluvium = category)) +
    geom_col(width = 0.4,color=NA)+
    geom_stratum(width=0.4, size = 0.5, alpha = 0.7, color="white") +
    geom_flow(width = 0.4, alpha = 0.5)+
    geom_text(aes(x = variable, y = value, label=round(100*value,2)
    # ifelse(value<0.006, "", round(value*100,2))
  ),
  position=position_stack(vjust=0.5), color="white",size=3)+
    scale_x_discrete(labels = c("T1","T2","T3","T4","T5","T6"))+
  scale_fill_manual(values = c("#B36F6D","#82C1C1","#497F7C","#D5BBA1","#FF5940"))+
  theme_bw()+
  ylab("Vd2+")+
  theme(axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid= element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key.size = unit(25, "pt"),
        legend.title =element_blank(),
        legend.text = element_text(size = 10,  color = "Black"),
       axis.text= element_text(size = 15,  color = "Black"),
      # # legend.position="none",
  )




########################### Fig 4D ###########################
# library(treemap)
library(voronoiTreemap)
library(wesanderson)
library(htmlwidgets)

PAT_all <- get_all(immdata_TRD_neg_cutoff)

### e.g T4 
tree <-  PAT_all  %>% filter(!is.na(Proportion_t4)) %>% 
        mutate( 
          class = case_when(
            type == "Pre-existed" ~ "Pre-existed",
            type == "TKI" ~ "TA",
            type == "Radio" ~ "RA",
            type == "PD-1" ~ "PA"
            )
          ) %>% 
          mutate(tag = paste0(CDR3.aa," ",class)) %>%
        select(CDR3.aa,Proportion_t4,class,tag)

    mycolor <- wes_palette("Darjeeling1", nrow(tree), type = "continuous")
    mycolor <- sample(mycolor, length(mycolor))

    vor <- data.frame(h1 = 'World', 
                  h2 = tree$class,
                  h3 = tree$CDR3.aa,
                  color = mycolor,
                  weight = tree$Proportion_t4,
                  codes = tree$tag
                  )

    vt <- vt_input_from_df(vor)

    # vt_d3(vt_export_json(vt),color_border = "#000000")

    myd3 <- vt_d3(vt_export_json(vt), label = FALSE,seed=5)

    saveWidget(myd3,"treemap/T4.html",selfcontained = TRUE)




############################ Fig 4E ############################
library(ComplexHeatmap)

lt = list(T1 =  immdata_TRD_neg_cutoff$data[[1]]$CDR3.aa %>% sort %>% unlist,
          T2 =  immdata_TRD_neg_cutoff$data[[2]]$CDR3.aa %>% sort %>% unlist,
          T3 =  immdata_TRD_neg_cutoff$data[[3]]$CDR3.aa %>% sort %>% unlist,
          T4 =  immdata_TRD_neg_cutoff$data[[4]]$CDR3.aa %>% sort %>% unlist,
          T5 =  immdata_TRD_neg_cutoff$data[[5]]$CDR3.aa %>% sort %>% unlist,
          t5 =  immdata_TRD_neg_cutoff$data[[6]]$CDR3.aa %>% sort %>% unlist
          )

list_to_matrix(lt)

m = make_comb_mat(lt,mode = "distinct")


pdf('test.pdf', height=5, width=12)

    UpSet(
    m,
    pt_size = unit(5, "mm"), lwd = 3,
    comb_col = c("#4D564D",
    #  "#4D564D",
    #  "#EE0000B2","#C7B458", "#42b540","#0099b4","#925e9f"
    "#02401B","#81A88D","#972D15","#D8B70A","#A2A475"
    )[comb_degree(m)],
    set_order = rev(c("T1", "T2", "T3","T4","T5","T6")),
    # comb_order = order.comb_mat(m, decreasing = TRUE),
    #     c("111111","110000","101000", "100100","100010","100001",
    # "010000", "010010","010001","010011",
    # "001000", "001010","001001","001011",
    # "000100", "000110","000101","000111",
    # "000010","000011","000001"
    # )),
    right_annotation = 
    upset_right_annotation(m, 
        # ylim = c(0, 200),
        gp = gpar(fill = "#008280B2"),
        annotation_name_side = "top",
        axis_param = list(side = "top"))
    )

dev.off()




########## Fig 4F ##########

PAT_all <- get_all(immdata_TRD_neg_cutoff)
# PAT_all <- get_all(immdata_TRD_pos_cutoff)
# PAT_all <- get_all(immdata_TRG_neg_cutoff)
# PAT_all <- get_all(immdata_TRG_pos_cutoff)


NEO <- PAT_all %>% filter(type %in% c("TKI","Radio","PD-1")) 


Neo <- rbind(
# IND %>% select(1,Proportion_t1) %>% rename("Proportion"=3) %>% select(3,4) %>% mutate(Time = "T1"),
NEO %>% select(CDR3.aa,Proportion_t2) %>% rename("Proportion"=2) %>% drop_na(Proportion) %>% mutate(Time = "T2"),
NEO %>% select(CDR3.aa,Proportion_t3) %>% rename("Proportion"=2) %>% drop_na(Proportion) %>% mutate(Time = "T3"),
NEO %>% select(CDR3.aa,Proportion_t4) %>% rename("Proportion"=2) %>% drop_na(Proportion) %>% mutate(Time = "T4"),
NEO %>% select(CDR3.aa,Proportion_t5) %>% rename("Proportion"=2) %>% drop_na(Proportion) %>% mutate(Time = "T5"),
NEO %>% select(CDR3.aa,Proportion_t6) %>% rename("Proportion"=2) %>% drop_na(Proportion) %>% mutate(Time = "t5")
  )%>% select(Proportion, Time) %>% mutate(Proportion = as.numeric(Proportion))


p<-ggplot(Neo, aes(x = Time, y = log(Proportion), color = Time))+
  geom_boxplot(aes(fill = Time), alpha = 0.5,
                outlier.size = 0, outlier.alpha = 0,
               position = position_dodge(0.8),
               width = 0.2,
               color = "black") + 
  geom_violin(aes(fill = Time), color="black", width = 0.6, trim = TRUE, alpha = 0.3)+
  geom_jitter(size = 1)+ 
  ylab("Proportion of Neo Clonotypes (log)")+
  guides(fill=guide_legend(title="Time"))+
  scale_color_npg()+
  scale_fill_npg(alpha = 0.7)+
  #  scale_y_log10()+
  theme_bw()+
  stat_compare_means(method = "kruskal.test", 
        label.x = 1,
         aes(label = paste0("p = ", ..p.format..))
        )+
  stat_compare_means(method = "wilcox.test",
  aes(label = ..p.signif..),size = 6,
  # label.y = -2.2,
  comparisons = list(
  #  c("T5","t5"),
  # c("T2","T5"),c("T3","T5"),c("T4","T5"),
    c("T2","t5"),c("T3","t5"),c("T4","t5"),c("T5","t5")
    ))+
  # annotate("text", x = 1.5, y = -1.8, size = 6,label ="paste(~italic(p), \" < 0.0001\")",parse = TRUE)+
  theme(legend.position="none",
    panel.border = element_rect(size=2,color = "black"),
    axis.text.x=element_text(size=15,angle=0,color="Black"),
    axis.text.y=element_text(size=15,angle=0,color="Black"),
    axis.title.y = element_text(size = 14,  color = "Black"),
    axis.title.x = element_text(size = 15,  color = "Black")
  ) 

