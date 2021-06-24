################################################################################
## draw_histoBGO
## : drawing histogram according to the group of selected GO group
##
## Function Needed
## : subNmodi_scoreONTO.R, sortGO.R, plot_histoBGO.R
##
## Example of Arguments
## : draw_histoBGO("Scerevisiae", "MF")
##
## Writer: KyungJin Sohn
## Date: 06.19.2020
## Last modified: 06.22.2021
################################################################################
draw_histoBGO <- function(species, ONTO) {
                # 1. Install the needed library
        if(!require(dplyr))
                install.packages('dplyr',repos = "http://cran.us.r-project.org")
        library(dplyr)
        source("./subNmodi_scoreONTO.R")
        source("./sortGO.R")
        source("./plot_histoBGO.R")
        
        
                # 2. read the score categorized by ontology
        scoreDatatGO <- subNmodi_scoreONTO(species, ONTO)
        
        
                # 3. Remove GO term that has less than 30 sequences
        scoreDatatGO <- scoreDatatGO %>% group_by(GO) %>% filter(n() >= 30)
        
        
                # 4. select top and bottom GO list
        sortReturn <- sortGO(species, ONTO)
        topGO <- sortReturn[[1]]
        bottomGO <- sortReturn[[2]]
        

                # 5. Draw histogram for each GO group
        plot_histoBGO(species, topGO[,1], scoreDatatGO, 
                      paste0(species, "_", ONTO, "_TOP10"))
        
        plot_histoBGO(species, bottomGO[,1], scoreDatatGO, 
                      paste0(species, "_", ONTO, "_BOTTOM10"))
        
        for (selGO in topGO[,1]) {
                group <- strsplit(selGO, split = '[:]')[[1]][2]
                plot_histoBGO(species, selGO, scoreDatatGO, 
                              paste0(species, "_", ONTO, "_TOP_GO", group))
        }
        
        for (selGO in bottomGO[,1]) {
                group <- strsplit(selGO, split = '[:]')[[1]][2]
                plot_histoBGO(species, selGO, scoreDatatGO, 
                              paste0(species, "_", ONTO, "_BOTTOM_GO", group))
        }
}