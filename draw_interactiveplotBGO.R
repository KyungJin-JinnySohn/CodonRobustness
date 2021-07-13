################################################################################
## change_CToAA.R
## : Draws an interactive line plot based on survival test results.
##
## Additional Function Needed
## : survivalTest_filter.R, combine_survivalTestBGO.R
##
## Example of Arguments
## : draw_interactiveplotBGO("HomoSapiens", "MF", "BLOSUM62", 1234)
##
## Writer: KyungJin Sohn
## Date: 03.18.2021
## Last modified: 07.13.2021
################################################################################
draw_interactiveplotBGO <- function(species, ONTO, subMatrix, seed) {
                # 1. Install the needed library
        if (!requireNamespace("BiocManager", quietly = TRUE)){
                install.packages("BiocManager")
        }
        if (!require("GO.db", character.only = TRUE)){
                BiocManager::install("GO.db")
        }
        if(!require(plotly))
                install.packages('plotly')
        if(!require(lubridate))
                install.packages('lubridate')
        if(!require(RColorBrewer))
                install.packages('RColorBrewer')
        if(!require(htmlwidgets))
                install.packages('htmlwidgets')
        library("GO.db")
        library(plotly)
        library(lubridate)
        library(RColorBrewer)
        library(htmlwidgets)
        
        source("./survivalTest_filter.R")
        source("./combine_survivalTestBGO.R")
        
        
                # 2. Read categorized data that matches several criteria
        survivalData <- survivalTest_filter(species, ONTO)
        
        
                # 3. Select  GO Terms and sort by mean sumScore
        GOData <- survivalData %>% group_by(GO) %>% 
                summarise(avg_sumScore = mean(sumScore)) %>% 
                mutate(term = Term(GO)) %>% 
                arrange(desc(avg_sumScore))
        
        idx <- c(1:26, 33:57)
        idx <- idx[!idx %in% c(15, 16, 38)]
        comResult <- combine_survivalTestBGO(GOData[idx, ]$GO,
                                             species, ONTO, subMatrix, seed)
        survivalData <- comResult[[1]]
        survivedAUC <- comResult[[2]]
        
        
                # 4. Join data
        GOData_survived <- left_join(GOData[idx,], survivedAUC, by = "GO")
        
        allSurvivedData <- left_join(survivalData, GOData_survived, by = "GO")
        
        
                # 5. Draw an interactive plot
        getAllPalette <- rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(nrow(GOData[idx,])))
        
        GOData <- GOData %>% mutate(GO = paste(GO, term)) %>% 
                arrange(desc(avg_sumScore))
        
        p <- allSurvivedData %>%
                mutate(GO = factor(paste(GO, term), levels = GOData[idx,]$GO)) %>% 
                group_by(GO) %>%
                plot_ly(x = ~ MutationNum,
                        text = ~paste('sumScore: ', round(avg_sumScore, 3),
                                      '</br>AUC: ', round(AUC, 3))) %>%
                add_lines(y = ~ survivedRate,
                          color = ~ GO,
                          colors = getAllPalette) %>%
                layout(title = paste0("Survival Test(w/", subMatrix, "): ", 
                                      species, "_", ONTO),
                       xaxis = list(title = "# of Mutation"),
                       yaxis = list(title = "Survival Rate"),
                       autosize = T, margin = 50)
        
        
                # 6. Save a plot
        if (!dir.exists(paste0(getwd(), "/plot/survivalTest/"))) {
                dir.create(paste0("./plot/survivalTest/"))
        }
        
        fileName <- paste0("./plot/survivalTest/", species, "_", ONTO, "_",
                           subMatrix, "_", seed, "_survivalTest.html")  
        
        saveWidget(p, fileName, selfcontained = F, libdir = "lib")
}

