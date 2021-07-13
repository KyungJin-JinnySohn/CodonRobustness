################################################################################
## change_CToAA.R
## : change a sequence of codons into amino acids
##
## Additional Function Needed
## : .
##
## Example of Arguments
## : change_CToAA("ATGTCCCTTG...")
##
## Writer: KyungJin Sohn
## Date: 02.25.2021
## Last modified: 07.07.2021
################################################################################
change_CToAA <- function(seq) {
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
        
        
                # 2. Change codon to amino acid
        dna <- strsplit(seq, "(?<=.{3})", perl = TRUE)[[1]]
        AA <- c()
        for (codon in dna) {
                if (str_length(codon) !=3) {
                        break
                } else {
                        AA <- c(AA, GENETIC_CODE[[codon]])
                }
        }
        AA <- paste(AA, collapse = "")
        
        return(AA)
}
