################################################################################
## cal_prob4eAI
## : Calculate probability p(c'|c)
##
## Function Needed
## : .
##
## Example of Arguments
## : "ATG"
##
## Writer: KyungJin Sohn
## Date: 06.01.2020
## Last modified: 06.07.2021
################################################################################

cal_prob4eAI <- function(codon) {
                # 1. Preparing libraries
        if(!require(readxl)) 
                install.packages('readxl',repos = "http://cran.us.r-project.org")
        
        library(readxl)
        
                # 2. Read the codon table
        codonTable <- read_excel("./value/codon_table.xlsx")
        
                # 3. Calculating probabiity
        codonTable$prob <- 0
        for (i in 1:nrow(codonTable)) {
                m_codon <- codonTable$codon[i]
                
                # convert codon to the vector
                v_codon <- strsplit(codon, "(?<=.{1})", perl = TRUE)[[1]]
                v_m_codon <- strsplit(m_codon, "(?<=.{1})", perl = TRUE)[[1]]
                
                # compare the difference
                com_codon <- v_codon == v_m_codon
                
                prob <- 0
                # if transitions or transversions have occured
                if (sum(com_codon) == 2) {
                        if (!com_codon[3]) {
                                prob <- 1
                        } else if (!com_codon[1]){
                                v_n <- v_codon[1]
                                v_m_c <- v_m_codon[1]
                                if ((v_n == "T" && v_m_c == "C") || 
                                    (v_n == "C" && v_m_c == "T") || 
                                    (v_n == "A" && v_m_c == "G") || 
                                    (v_n == "G" && v_m_c == "A")) {
                                        prob <- 1
                                } else {
                                        prob <- 0.5
                                }
                        } else {
                                v_n <- v_codon[2]
                                v_m_c <- v_m_codon[2]
                                if ((v_n == "T" && v_m_c == "C") || 
                                    (v_n == "C" && v_m_c == "T") || 
                                    (v_n == "A" && v_m_c == "G") || 
                                    (v_n == "G" && v_m_c == "A")) {
                                        prob <- 0.5
                                } else {
                                        prob <- 0.1
                                }
                        }
                }
                
                # save the probability
                codonTable$prob[i] <- prob
        }
        codonTable$prob <- codonTable$prob/sum(codonTable$prob)
        
                # 4. return
        return(codonTable[, c("codon", "prob")])
}
