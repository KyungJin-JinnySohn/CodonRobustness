################################################################################
## cal_substitution_score.R
## : Calculate a substitution score after one random mutation occurrence
##
## Additional Function Needed
## : .
##
## Example of Arguments
## : cal_substitution_score("ATGTCCCTTG...", subMatrix, 23)
##
## Writer: KyungJin Sohn
## Date: 02.25.2021
## Last modified: 07.07.2021
################################################################################
cal_substitution_score <- function(seqData, subMatrix, randomSpot) {
                # 1. Install the needed library
        if(!require(stringr)) 
                install.packages("stringr",repos = "http://cran.us.r-project.org")
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
                install.packages("BiocManager", repos = "http://cran.us.r-project.org")
        }
        if (!require("Biostrings", character.only = TRUE)) {
                BiocManager::install("Biostrings")
        }
        library(stringr)
        library(Biostrings)
        
        
                # 2. Prepare val. to check the condition
        quo <- randomSpot %/% 3
        rem <- randomSpot %% 3
        
        
                # 3. Split the sequence by codon
        mutSeq <- strsplit(seqData$Seqs, "(?<=.{3})", perl = TRUE)[[1]]
        
        
                # 4. Identify the original amino acid in the location 
                #    & the codon in the current location and its amino acid
        if (rem == 0) {
                oriAA <- substr(seqData$AA, quo, quo)
                mutC1 <- mutSeq[quo]
        } else {
                oriAA <- substr(seqData$AA, quo+1, quo+1)
                mutC1 <- mutSeq[quo + 1]
        }
        mutAA1 <- GENETIC_CODE[[mutC1]]
        
        
                # 5. Select the nucleotide randomly and replace it
        codon <- c ("A", "C", "G", "T")
        if (rem == 0) {
                codon <- codon[!codon %in% substr(mutC1, 3, 3)]
                newC <- sample(codon, 1, replace = TRUE)
                substr(mutC1, 3, 3) <- newC
        } else {
                codon <- codon[!codon %in% substr(mutC1, rem, rem)]
                newC <- sample(codon, 1, replace = TRUE)
                substr(mutC1, rem, rem) <- newC
        }
        
        
                # 6. Identify the current amino acid in the location after mutation
        mutAA2 <- GENETIC_CODE[[mutC1]]
        
        
                # 7. Calculate the score and return the result
        subScore <- seqData$mutScore - subMatrix[oriAA, mutAA1] + subMatrix[oriAA, mutAA2]
        
        return(list("subScore" = subScore, "newC" = newC))
}