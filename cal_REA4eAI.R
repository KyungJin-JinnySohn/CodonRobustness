################################################################################
## cal_REA4eAI
## : Calculate the REA for a specific species
##
## Function Needed
## : cal_prob4eAI.R
##
## Example of Arguments
## : "Scerevisiae"
##
## Writer: KyungJin Sohn
## Date: 06.01.2020
## Last modified: 06.07.2021
################################################################################

cal_REA4eAI <- function(species) {
                # 1. Preparing libraries
        if(!require(readxl)) 
                install.packages('readxl',repos = "http://cran.us.r-project.org")
        if(!require(dplyr)) 
                install.packages('dplyr',repos = "http://cran.us.r-project.org")
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
                install.packages("BiocManager", 
                                 repos = "http://cran.us.r-project.org")
                BiocManager::install("Biostrings")
        }
        if(!require(stringr)) 
                install.packages("stringr",repos = "http://cran.us.r-project.org")
        
        library(readxl)
        library(dplyr)
        library(Biostrings)
        library(stringr)
        source("cal_prob4eAI.R")
        
                # 2. Read cost data, codon table, and tWeight data
        cost <- read.table("./value/AAindex_thermo.txt", sep = "\t", 
                           header = TRUE)
        codonTable <- read_excel("./value/codon_table.xlsx")
        tW <- read.csv(paste("./value/RTA/RTA_", species,".csv", sep = ""), 
                       header = TRUE, stringsAsFactors = FALSE)
        
        
                # 3. Join the codon table and tWeight
        mainData <- inner_join(codonTable[, c("codon", "AA")], tW, by = "codon")
        
                # !!!! CAUTION !!!!
                #       Since tWeight is used as a denominator, 
                #       small number are given to prevent error.
                #       May fix, if there is a better way.
        mainData[which(mainData$tWeight == 0), ]$tWeight <- 0.001
        
                # 4. Calculate RTA for each codon
        mainData <- mainData %>% mutate(RTA = tWeight/max(tWeight))
        
        
                # 5. Calculate load using fitness function
        mainData$load <- 0
        for (j in 1:nrow(mainData)) {
                codon <- mainData$codon[j]
                pRTA <- mainData$RTA[j]
                probList <- cal_prob4eAI(codon)
                
                load <- 0
                for (i in 1:nrow(codonTable)) {
                        m_codon <- codonTable$codon[i]
                        
                        subCost <- cost[GENETIC_CODE[[codon]], 
                                        GENETIC_CODE[[m_codon]]]
                        if (is.null(subCost)) {
                                subCost <- max(cost)
                        }
                        
                        prob <- probList[probList$codon == m_codon, ]$prob
                        
                        load <- load + subCost*prob/pRTA
                }
                
                mainData$load[j] <- load
        }
        
        
                # 6. Caluclate REA using the load
        mainData <- mainData %>% group_by(AA) %>% 
                mutate(REA = (1.01*max(load) - load)/(1.01*max(load) - min(load))) %>% 
                arrange(AA, REA, codon)
        mainData$REA <- round(mainData$REA, 3)
        
        
                # 7. Save the result
        write.csv(mainData[, c("codon", "REA")], 
                  file = paste("./value/REA/REA_", species, ".csv", sep = ""), 
                  row.names = FALSE) 
}
