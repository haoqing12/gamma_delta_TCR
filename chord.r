library(tidyverse)
library(data.table)

##### Fig 1D #####

df <- data.frame()
for (n in 1:nrow(immdata_TRG$meta)){
  print(n)
  df <- rbind(immdata_TRG$data[[n]], df)
}

df <- df %>% select(Clones,V.name,J.name) %>% 
  separate(V.name,into=c('V',NA),sep='\\*') %>%
  separate(J.name,into=c('J',NA),sep='\\*')

head(df)

V <- sort(unique(unlist(strsplit(paste(df$V, collapse=","), ",")))) %>% strsplit(split="\\*00") %>% unlist
J <- sort(unique(unlist(strsplit(paste(df$J, collapse=","), ",")))) %>% strsplit(split="\\*00") %>% unlist

# create and fill matrix of pairwise with sum(Clones)
mat <- matrix(nrow = length(J), ncol = length(V))
colnames(mat) <- V
rownames(mat) <- J
for (j in J) {
  for (v in V) {
  # sub <- df[(stringr::str_detect(df$V, V[[v]], negate = FALSE) & 
  #              stringr::str_detect(df$J, J[[j]], negate = FALSE)),]
  mat[j,v] <- df %>% filter(V == v & J == j) %>% pull(Clones) %>% sum()
  }
}



# convert to proportions / fraction of total
mat <- mat/sum(mat)

# reorder by V and J decreasing order (left to right)
vmax <- colSums(mat)
jmax <- rowSums(mat) 
mat2 <- mat[order(jmax, decreasing = FALSE), order(vmax, decreasing = TRUE)]


library(circlize)
circos.clear()
circos.par("canvas.xlim" = c(-1.25, 1.25), "canvas.ylim" = c(-1.25, 1.25))

pdf(file="1.pdf", width=9, height=5, pointsize=8)
# use circlize
chordDiagram(mat2, annotationTrack = "grid")

# add legends
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, 
              CELL_META$ylim[1], 
              CELL_META$sector.index, 
              facing = "clockwise", 
              niceFacing = TRUE,
              adj = c(-0.25, 0.5))
}, bg.border = NA)

dev.off()


