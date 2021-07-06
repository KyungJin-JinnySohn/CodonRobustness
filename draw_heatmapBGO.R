################################################################################
## draw_heatmapBGO
## : draw a heatmpa with categorized score data
## 
## Function Needed
## : subNmodi_scoreONTO.R, plot_heatmapBGO.R
##
## Example of Arguments
## : draw_heatmapBGO(
##      c("HomoSapiens", "MusMusculus", "Dmelanogaster", "Celegans", "Scerevisiae", "Ecoli"), 
##      "BP", "rank") # ONTO: "BP"/"CC"/"MF", indexType: "rank"/"score"
##
## Writer: KyungJin Sohn
## Date: 09.07.2020
## Last modified: 07.06.2021
################################################################################
draw_heatmapBGO <- function(species, ONTO, indexType ) {
                # 1. Throw an error if the length of the species vector is less than 2
        if (length(species) < 2) {
                stop("ERROR: Please select more than one species.")   
        }
        
        
                # 2. Install the needed library
        if(!require(dplyr)){
                install.packages('dplyr',repos = "http://cran.us.r-project.org")
        }
        if(!require(tidyr)){
                install.packages('tidyr',repos = "http://cran.us.r-project.org")
        }
        if (!requireNamespace("BiocManager", quietly = TRUE)){
                install.packages("BiocManager")
        }
        if (!require("GO.db", character.only = TRUE)){
                BiocManager::install("GO.db")
        }
        library(dplyr)
        library(tidyr)
        library("GO.db")
        source("./subNmodi_scoreONTO.R")
        source("./plot_heatmapBGO.R")
        
        
                # 3. Read the categorized score data for each species
        heatmapData <- ""
        for (i in 1:length(species)) {
                        # read the categorized score data
                scoreDatatGO <- subNmodi_scoreONTO(species[i], ONTO)
                
                        # removes GO terms with fewer than 30 sequences 
                        # & summarizes the average of frequency and scores.
                scoreDatatGO <- scoreDatatGO %>% group_by(GO) %>% 
                        filter(n() >= 30) %>% 
                        summarise(GOfreq = n(), 
                                  wR2N = mean(wR2N), CRI2 = mean(CRI2), 
                                  eAI = mean(eAI), sumScore = mean(sumScore))

                        # join all score data
                if (i == 1) {
                        heatmapData <- scoreDatatGO
                } else {
                        heatmapData <- inner_join(heatmapData, scoreDatatGO, 
                                                  by = "GO")
                }
        }
        
        
                # 4. Prepare detail description for each GO term
                        # change data format
        heatmapData <- as.data.frame(heatmapData)
        
                        # concatenate Go freq columns
        heatmapData <- heatmapData %>% 
                unite("GOfreq", seq(2, ncol(heatmapData), by = 5), 
                      sep = ", ", remove = TRUE)
        
                        # add a column of GO term description
        heatmapData$GO_func <- Term(heatmapData$GO)
        
                        # combine the two columns to store them as row names
        heatmapData <- mutate(heatmapData, 
                              GO = paste0(GO, " (", GO_func,")", " [", GOfreq, "]"))
        rownames(heatmapData) <- heatmapData$GO
        heatmapData <- subset(heatmapData, select = -c(GO, GOfreq, GO_func))
        
        
        
                # 5. Draw heatmap for each score
        if (!dir.exists(paste0(getwd(), "/plot/heatmap/"))) {
                dir.create(paste0("./plot/heatmap/"))
                dir.create(paste0("./plot/heatmap/rank/"))
                dir.create(paste0("./plot/heatmap/score/"))
        }

        scores <- c("wR2N", "CRI2", "eAI", "sumScore")
        for (i in 1:length(scores)) {
                scoreData <- heatmapData[, seq(i, ncol(heatmapData), by = 4)]
                colnames(scoreData) <- species
                plot_heatmapBGO(species, ONTO, scores[i], scoreData, indexType)
        }
        
}
