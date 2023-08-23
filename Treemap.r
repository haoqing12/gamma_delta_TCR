library(immunarch)
library(voronoiTreemap)
library(wesanderson)
library(htmlwidgets)


getcutoff <- function(immdata_TRD){
    immdata_TRD_cutoff <- list()
    for (i in immdata_TRD$meta$Sample){
      print(i)
      cutoff_ <- immdata_TRD$data[[i]]$Clones %>% sum/nrow(immdata_TRD$data[[i]])
      immdata_TRD_cutoff$data[[i]] <- immdata_TRD$data[[i]] %>% filter(Clones >= cutoff_/2)
    }
    immdata_TRD_cutoff$meta <- immdata_TRD$meta

  return(immdata_TRD_cutoff)
}


immdata_TRD_neg_cutoff <- getcutoff(immdata_TRD_neg)
immdata_TRD_pos_cutoff <- getcutoff(immdata_TRD_pos)
immdata_TRG_neg_cutoff <- getcutoff(immdata_TRG_neg)
immdata_TRG_pos_cutoff <- getcutoff(immdata_TRG_pos)




gettreemap <- function(immdata,P_N) {
    
for (i in   immdata$meta$Sample){
    print(i)
    tree <- immdata$data[[i]] %>%
            arrange(desc(Proportion)) %>%
            dplyr::select(CDR3.aa,Proportion) %>% 
            mutate(Proportion = Proportion*100) %>% 
            distinct 

    mycolor <- wes_palette("Darjeeling1", nrow(tree), type = "continuous")
    mycolor <- sample(mycolor, length(mycolor))

    vor <- data.frame(h1 = 'World', 
                  h2 = 'TRX',
                  h3 = tree$CDR3.aa,
                  color = mycolor,
                  weight = tree$Proportion,
                  codes = tree$CDR3.aa
                  )

    vt <- vt_input_from_df(vor)

    myd3 <- vt_d3(vt_export_json(vt),seed = 123, label = FALSE)

    saveWidget(myd3,paste0("treemap2/",P_N,"/",i,"_.html"),selfcontained = TRUE)

}
}


# gettreemap(immdata_TRG_cutoff,"TRG_cutoff")
# gettreemap(immdata_TRD_cutoff,"TRD_cutoff")

gettreemap(immdata_TRG_pos_cutoff,"TRG_pos_cutoff")
gettreemap(immdata_TRG_neg_cutoff,"TRG_neg_cutoff")

gettreemap(immdata_TRD_pos_cutoff,"TRD_pos_cutoff")
gettreemap(immdata_TRD_neg_cutoff,"TRD_neg_cutoff")
