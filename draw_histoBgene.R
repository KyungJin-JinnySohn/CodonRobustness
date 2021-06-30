################################################################################
## draw_histoBgene
## : draw a histogram based on genes
##
## Function Needed
## : subNmodi_score.R, plot_histoBgene.R
##
## Example of Arguments
## : draw_histoBgene("Scerevisiae")
##
## Writer: KyungJin Sohn
## Date: 09.09.2020
## Last modified: 06.29.2021
################################################################################
draw_histoBgene <- function(species) {
                # 1. Install the needed library
        source("./subNmodi_score.R")
        source("./plot_histoBgene.R")
        
        
                # 2. Read score data of species
        scoreData <- subNmodi_score(paste0("Score_", species))
        
        
                # 3. Call function to draw plot for each score
        randomType <- c("GC", "AA")     
        for (type in randomType) {
                plot_histoBgene(scoreData, species, type)
        }
        
}