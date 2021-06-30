################################################################################
## combine_ONTO_data
## : Combine data categorized by ontology
##
## Function Needed
## : .
##
## Example of Arguments
## : combine_ONTO_data("HomoSapiens")
##
## Writer: KyungJin Sohn
## Date: 12.28.2020
## Last modified: 06.22.2021
################################################################################
combine_ONTO_data <- function(species) {
                # 1. Install the needed library
        if(!require(dplyr))
                install.packages('dplyr',repos = "http://cran.us.r-project.org")
        if(!require(stringr)) 
                install.packages("stringr",repos = "http://cran.us.r-project.org")
        
        library(dplyr)
        library(stringr)
        
                # 2. Read and bind data categorized by ontology
        scoreData <- data.frame()
        for (ONTO in c("BP", "CC", "MF")) {
                openfile <- gzfile(paste0("./scoreResult/", species, "_Score_", 
                                          ONTO, ".csv.gz"))
                data <- read.csv(openfile, header = TRUE, row.names = 1, 
                                 stringsAsFactors = FALSE)
                #close(openfile)
                
                scoreData <- rbind(scoreData, data)
        }
        
                # 3. select the needed column
        scoreData <- scoreData[, !names(scoreData) %in% c("GO", "ONTOLOGY")]

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