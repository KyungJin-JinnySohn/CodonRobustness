################################################################################
## combine_random_data
## : combine random score data w/ diff. seed number
##
## Function Needed
## : subNmodi_score.R
##
## Example of Arguments
## : combine_random_data("Scerevisiae", "GC)
##
## Writer: KyungJin Sohn
## Date: 09.09.2020
## Last modified: 06.29.2021
################################################################################
combine_random_data <- function(species, randomType) {
                # 1. Install the needed libraries and funtions
        source("./subNmodi_score.R")
        
        
                # 2. Read data and bind it
        random_scoreData <- data.frame()
        #seeds <- c(1000, 2000, 3000, 4000, 5000, 1111, 2222, 3333, 4444, 5555)
        # change the seeds vector based on the seed numbers used when creating random data
        seeds <- c(1111)
        for (num in seeds) {
                # 1) Read score data of random sequence
                if (randomType == "AA") {
                        random_file <- paste0("random/Score_", species,
                                              "_randomAA_", num)   
                } else {
                        random_file <- paste0("random/Score_", species,
                                              "_randomGC_", num)
                }
                
                scoreData <- subNmodi_score(random_file)
                
                # 2) Save the result
                random_scoreData <- rbind(random_scoreData, scoreData)
        }
        
        return(random_scoreData)
}