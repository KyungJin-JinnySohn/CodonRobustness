################################################################################
## sortGO
## : Select top and bottom 10% GO term 
##
## Function Needed
## : .
##
## Example of Arguments
## : sortGO("Scerevisiae", "MF")
##
## Writer: KyungJin Sohn
## Date: 06.19.2020
## Last modified: 06.22.2021
################################################################################
sortGO <- function(species, ONTO) {
                # 1. Install the needed library
        if(!require(dplyr))
                install.packages('dplyr',repos = "http://cran.us.r-project.org")
        library(dplyr)
        
        
                # 2. Read the summary data
        scoreMeanGO <- read.csv(paste0("./summary/GO/", species, "_ScoreMean_", 
                                       ONTO, ".csv"), 
                                header = TRUE, row.names = 1, 
                                stringsAsFactors = FALSE)
        
        
                # 3. Rank the GO group by the score
        scoreMeanGO <- scoreMeanGO %>% mutate(wR2NRank = rank(-wR2Nmean), 
                                              CRI2Rank = rank(-CRI2mean),
                                              eAIRank = rank(-eAImean),
                                              sumScoreRank = rank(-sumScoremean))
        
        
                # 4. Select top and bottom 10% GO term 
        tenPercentage <- round(nrow(scoreMeanGO) * 0.1)
        
        topGO <- scoreMeanGO %>% 
                filter(wR2NRank <= tenPercentage,
                       CRI2Rank <= tenPercentage,
                       eAIRank <= tenPercentage,
                       sumScoreRank <= tenPercentage) %>% 
                select(GO)
        
        bottomGO <- scoreMeanGO %>% 
                filter(wR2NRank > nrow(scoreMeanGO) - tenPercentage,
                       CRI2Rank > nrow(scoreMeanGO) - tenPercentage,
                       eAIRank > nrow(scoreMeanGO) - tenPercentage,
                       sumScoreRank > nrow(scoreMeanGO) - tenPercentage) %>% 
                select(GO)
        
        
        return(list(topGO, bottomGO))
}