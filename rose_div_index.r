library(tidyverse)
library(data.table)
library(stringr)


dd <- div_trd_pos %>% mutate(Type = str_extract(sample, "([a-z]*[A-Z]+)")) %>% 
  filter(Type != "MTC") %>%
  select(1,Type,Clonality,Gini_coefficient,shannon_entropy,D50) %>% 
  mutate(
      Clonality = Clonality*10,
      Gini_coefficient = scale(Gini_coefficient)+4,
      shannon_entropy = shannon_entropy,
      D50 = log(D50)*2
  ) %>%
    melt() %>% 
    group_by(Type,variable) %>%
    summarize(means = mean(value), SE =  sd(value)/sqrt(length(value))) %>%
    arrange(variable) %>% as.data.frame() %>%
    mutate(s = paste0(Type,"_",variable))


### order ###
l <- c("HC","miPTC","wiPTC","PDTC","ATC")
s <- c("Clonality","Gini_coefficient","shannon_entropy","D50")
dd <- dd %>% slice(match(c(
                          paste0(l,"_D50"),
                          paste0(l,"_shannon_entropy"),
                          paste0(l,"_Gini_coefficient"),
                          paste0(l,"_Clonality"))
                        ,s))
dd$s <- factor(dd$s, levels = dd$s)

myAngle <- seq(-10,-350,length.out =20)



pdf('test.pdf', height=10, width=10)

ggplot(dd, aes(x=variable, y=means, fill= s))+
  geom_bar(color = "white",stat="identity",position="dodge",width = 0.8, alpha = 0.8,size=1)+
  # scale_y_continuous(breaks = 0:nlevels(data$Trajectory)) +
  geom_text(position = position_dodge(.8), 
    aes(y = means+1.1, label = Type),
    alpha = 1, size = 4, vjust = 0.5, hjust = 0,angle = myAngle)+
  geom_errorbar(aes(ymin=means-SE, ymax=means+SE), position=position_dodge(.8),
    width = 0.3,size = 1)+
  coord_polar() +
  ylim(-1,10)+
  scale_fill_manual(values = c("#ECA841","#BD8A5F","#FEB23F","#EFCE71","#EBE08F",
                               "#949F1E","#BCD784","#70B136","#8BC48F","#2f7d3e",
                               "#C1D7FC","#9CAABD","#7D8FAE","#9599C0","#8A8FC9",
                               "#72BBD9","#0BD0C4","#0A9CBD","#11CAF6","#11AAFC"))+
  theme_minimal()+
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.line = element_blank(),
    axis.title = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
  )+
  guides(fill = 'none')

dev.off()



