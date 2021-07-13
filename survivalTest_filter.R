################################################################################
## survivalTest_filter.R
## : Select GO Terms that meets the conditions from the categorized score data
##
## Additional Function Needed
## : subNmodi_scoreONTO.R
##
## Example of Arguments
## : survivalTest_filter("HomoSapiens", "MF")
##
## Writer: KyungJin Sohn
## Date: 02.25.2021
## Last modified: 07.07.2021
################################################################################
survivalTest_filter <- function(species, ONTO) {
                # 1. Install the needed library
        if(!require(dplyr))
                install.packages('dplyr',repos = "http://cran.us.r-project.org")
        if(!require(stringr)) 
                install.packages("stringr",repos = "http://cran.us.r-project.org")
        library(dplyr)
        library(stringr)
        source("./subNmodi_scoreONTO.R")

        
                # 2. Read the score categorized by ontology
        scoreDatatGO <- subNmodi_scoreONTO(species, ONTO)
        
        
                # 3. Add the length of sequence column
        scoreDatatGO <- scoreDatatGO %>% mutate(seqLength = str_length(Seqs))
        
        
                # 4. Delete and filter data
                        # delete a sequence that is not a multiple of 3 in length
        scoreDatatGO <- scoreDatatGO %>% filter(seqLength %% 3 == 0)
        
                        # find out the mid-20th quantile range.
        q40 <- quantile(scoreDatatGO$seqLength, probs = .40)
        q60 <- quantile(scoreDatatGO$seqLength, probs = .60)
        
                        # filter out GOs that has less than 30 sequences
                        # & select GO only if at least five sequences fall 
                        #   within the range specified above
        survivalDataGO <- scoreDatatGO %>% group_by(GO) %>% 
                filter(n() >= 30, sum(seqLength >= q40) - sum(seqLength > q60) >= 5,
                       seqLength >= q40, seqLength <= q60)
        
        return(survivalDataGO)
}