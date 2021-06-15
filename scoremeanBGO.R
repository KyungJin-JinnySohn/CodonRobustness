################################################################################
## scoremeanBGO
## : Compute score for each GO means
##
## Function Needed
## : subNmodi_scoreONTO.R
##
## Example of Arguments
## : scoremeanBGO("HomoSapiens", "BP")
##
## Writer: KyungJin Sohn
## Date: 07.06.2020
## Last modified: 06.15.2021
################################################################################

scoremeanBGO <- function(species, ONTO) {
                # 1. Install the needed library
        if(!require(dplyr)) 
                install.packages('dplyr',repos = "http://cran.us.r-project.org")
        library(dplyr)
        source("subNmodi_scoreONTO.R")
        
        
        #scoreData <- subNmodi_scoreONTO("HomoSapiens", "BP")
                # 2. read the (species)_Score_(ONTO) data
        scoreData <- subNmodi_scoreONTO(species, ONTO)
        
        
                # 3. Delete rows matching GO terms used less than 30 times
        scoreData <- scoreData %>% group_by(GO) %>% filter(n() >= 30)
        
        
                # 4. Compute score mean for each GO term
        meanScoreData <- scoreData %>% group_by(GO) %>% 
                summarise(wR2Nmean = mean(wR2N), CRI2mean = mean(CRI2),
                          eAImean = mean(eAI), sumScoremean = mean(sumScore))
        
        
                # 5. Save the result
        write.csv(meanScoreData, paste("./summary/GO/", species, "_ScoreMean_", 
                                       ONTO, ".csv", sep = ""))
}
