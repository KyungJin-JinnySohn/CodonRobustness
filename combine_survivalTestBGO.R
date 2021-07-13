################################################################################
## combine_survivalTestBGO.R
## : Combine the survival test result data
##
## Additional Function Needed
## : .
##
## Example of Arguments
## : combine_survivalTestBGO(c("GO0044235", "GO05249"), "MF", "BLOSUM62", 1234)
##
## Writer: KyungJin Sohn
## Date: 03.18.2021
## Last modified: 07.13.2021
################################################################################
combine_survivalTestBGO <- function(GOvec, species, ONTO, subMatrix, seed) {
                # 1. Install the needed library
        if(!require(dplyr)) 
                install.packages("dplyr",repos = "http://cran.us.r-project.org")
        library(dplyr)
        
        
                # 2. Combine survival test data
        survivalData <- data.frame()
        survivedAUC <- data.frame()
        for (GOnum in GOvec) {
                openfile <- gzfile(paste0("./survivalResult/", species, "/",
                                          ONTO, "/", sub(":", "", GOnum), "_", 
                                          subMatrix, "_", seed, ".csv.gz"))
                survivalNum_GO <- read.csv(openfile, header = TRUE, row.names = 1, 
                                           stringsAsFactors = FALSE)
                
                survivalNum_GO <- survivalNum_GO %>% group_by(MutationNum) %>% 
                        summarise(survivedRate = sum(survivedNum)/n_distinct(survivalNum_GO$Names))
                
                survivalNum_GO <- cbind(GO = GOnum, survivalNum_GO)
                
                survivalData <- rbind(survivalData, survivalNum_GO)
        }
        
                # 3. Calculate AUC
        survivedAUC <- survivalData %>% group_by(GO) %>% 
                summarise(AUC = sum(survivedRate)/n())
        
        
                # 4. Return the result 
        return(list(survivalData, survivedAUC))
}