################################################################################
## subNmodi_score
## : Subsetting & Modifying score data to make readable way
##
## Function Needed
## : .
##
## Example of Arguments
## : subNmodi_score("HomoSapiens")
##
## Writer: KyungJin Sohn
## Date: 07.06.2020
## Last modified: 06.15.2021
################################################################################
subNmodi_score <- function(fileName) {
        # 1. Install the needed library
        if(!require(dplyr)) 
                install.packages('dplyr',repos = "http://cran.us.r-project.org")
        library(dplyr)
        
        
        # 2. read the (species)_Score_(ONTO) data
        scoreData <- gzfile(paste0("./scoreResult/", fileName, ".csv.gz"))
        scoreData <- read.csv(scoreData, header = TRUE, row.names = 1, 
                              stringsAsFactors = FALSE)
        
        
        # 3. Select the needed columns
        scoreData <- scoreData[, c("Names", "Seqs", "wR2N", "CRI2", "eAI")]
        
        # Delete rows with NA values
        scoreData <- scoreData[complete.cases(scoreData), ]
        
        # Delete duplicated rows
        scoreData <- scoreData[!duplicated(scoreData), ]
        
        # 4. Normalize score:
        #       wR2N = (wR2N - min(wR2N))/(max(wR2N) - min(wR2N))
        scoreData <- scoreData %>% mutate(wR2N = (wR2N + 1)/2)
        
        # Flip and shift the score for convenience
        scoreData$wR2N <- -scoreData$wR2N + 1
        scoreData$CRI2 <- -scoreData$CRI2 + 1
        
        # 5. Add new score column
        scoreData <- scoreData %>% mutate(sumScore = wR2N + CRI2 + eAI)
        
        # Normalize score:
        #       sumScore = (sumScore - min(sumScore))/(max(sumScore) - min(sumScore)
        scoreData <- scoreData %>% mutate(sumScore = sumScore/3)
        
        return(scoreData)
}