################################################################################
## cal_robustness_score
## : Calculate the robustness score of genes
##
## Function Needed
## : cal_REA4eAI.R, cal_prob4eAI.R
##
## Example of Arguments
## : list("inputData" = "Scerevisiae.fasta.gz", "outputData" = "Scerevisiae",
##      "species" = "Scerevisiae")
##
## Writer: KyungJin Sohn
## Date: 06.01.2020
## Last modified: 06.11.2021
################################################################################
cal_robustness_score <- function(args) {
        
                # 1. Preparing libraries
        if(!require(seqinr)) 
                install.packages('seqinr',repos = "http://cran.us.r-project.org")
        if(!require(readxl)) 
                install.packages('readxl',repos = "http://cran.us.r-project.org")
        if(!require(dplyr)) 
                install.packages('dplyr',repos = "http://cran.us.r-project.org")
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
                install.packages("BiocManager", 
                                 repos = "http://cran.us.r-project.org")
                BiocManager::install("Biostrings")
        }
        
        library(stringr) 
        library(seqinr)
        library(readxl)
        library(dplyr)
        library(Biostrings)
        source("cal_REA4eAI.R")
        
                # 2. Read a sequence data
        #con <- gzfile(paste("./inSilico/", args$inputData, sep = ""))
        con <- gzfile(paste("./data/", args$inputData, sep = ""))
        entireData <- read.fasta(file = con, as.string = TRUE, 
                                 forceDNAtolower = FALSE, whole.header = TRUE)
        close(con)
        entireData <- data.frame(Names = names(entireData), 
                                 Seqs = unlist(getSequence(entireData, 
                                                           as.string=T)), 
                                 R1N = 0, R2N = 0, wR1N = 0, wR2N = 0, CRI1 = 0, 
                                 CRI2 = 0, eAI = 0, stringsAsFactors = FALSE)
        
                # 3. Read the MD value
        MD <- read_excel("./value/MDvalue.xlsx", col_names = TRUE)
        
                # 4. Read the REA value according to the species
        species <- args$species
        
                # If the REA value does not exist, calculate and make the file.
        if(!file.exists(paste("./value/REA/REA_", species, ".csv", sep = ""))){
                cal_REA4eAI(species)
        }
        
        REA <- read.csv(paste("./value/REA/REA_", species, ".csv", sep = ""), 
                        header = TRUE, stringsAsFactors = FALSE)
        colnames(REA) <- c("Codon", "REA")
        
                # 5. Start calculating robustness score
        for (i in 1:nrow(entireData)) {
                # get the codon usage frequency of a sequence
                seq <- strsplit(entireData[i, 2], "(?<=.{3})", perl = TRUE)[[1]]
                seqFreq <- as.data.frame(table(seq), stringsAsFactors = FALSE)
                colnames(seqFreq) <- c("Codon", "Freq")
                
                # join the two data.frame
                mainData <- inner_join(MD, seqFreq, by = "Codon")
                
                # calculate codon usage fraction in synonmous codon family
                mainData <- mainData %>% group_by(AA) %>% 
                        mutate(Frac = Freq/sum(Freq))
                
                # calculate amino acid usage fraction in a gene
                aaSeq <- translate(DNAString(entireData[i, 2]), 
                                   genetic.code = GENETIC_CODE, 
                                   if.fuzzy.codon="error")
                aaFrac <- as.data.frame(alphabetFrequency(aaSeq, as.prob = TRUE), 
                                        stringsAsFactors = FALSE)
                colnames(aaFrac) <- "aaFrac"
                aaFrac$AA1 <- rownames(aaFrac)
                
                if (nrow(mainData) == 0) {
                        entireData[i, 3:9] <- -2
                } else {
                        ### 1) R1, R2
                        # calculate the correlation between MD and coodn frequencies in synonmous codon family
                        RData <- mainData %>% group_by(AA) %>% 
                                summarise(R1 = cor(MD1, Frac), 
                                          R2 = cor(MD2, Frac), .groups = 'drop')
                        
                        # calculate RN (N<18 if for some aa there is no variance in the MD values or in the freq of its syn codons)
                        entireData$R1N[i] <- mean(RData$R1, na.rm = TRUE)
                        entireData$R2N[i] <- mean(RData$R2, na.rm = TRUE)
                        
                        ### 2) wR1, wR2
                        RData$AA1 <- a(RData$AA)
                        RData <- left_join(RData, aaFrac, by = "AA1")
                        
                        entireData$wR1N[i] <- weighted.mean(RData$R1, 
                                                            RData$aaFrac, 
                                                            na.rm = TRUE)
                        entireData$wR2N[i] <- weighted.mean(RData$R2, 
                                                            RData$aaFrac, 
                                                            na.rm = TRUE)
                        
                        ### 3) CRI
                        # calculate C for each codon using MD1 and MD2
                        CData <- mainData %>% group_by(AA) %>% 
                                mutate(C1 = Frac*((MD1 - min(MD1))/(max(MD1) - min(MD1))), 
                                       C2 = Frac*((MD2 - min(MD2))/(max(MD2) - min(MD2))))
                        
                        # sum C for each amino acid
                        CData <- CData %>% group_by(AA) %>% 
                                summarise(C1m = sum(C1, na.rm = TRUE), 
                                          C2m = sum(C2, na.rm = TRUE), 
                                          .groups = 'drop')
                        
                        CData$AA1 <- a(CData$AA)
                        CData <- left_join(CData, aaFrac, by = "AA1")
                        
                        # calculated weighted mean of C for each gene
                        entireData$CRI1[i] <- weighted.mean(CData$C1m, 
                                                            CData$aaFrac, 
                                                            na.rm = TRUE)
                        entireData$CRI2[i] <- weighted.mean(CData$C2m, 
                                                            CData$aaFrac, 
                                                            na.rm = TRUE)
                        
                        ### 4) eAI
                        mainData <- inner_join(seqFreq, REA, by = "Codon")
                        
                        # multiply each e^(freq.)
                        mainData <- mainData %>% mutate(e = REA^Freq)
                        
                        # calculate eAI for each gene
                        eAI <- 1
                        for (j in mainData$e) {
                                eAI <- eAI * j^(1/sum(mainData$Freq))
                        }
                        
                        entireData$eAI[i] <- eAI
                }
        }
        
        
        gz1 <- gzfile(paste("./scoreResult/Score_", args$fileName, 
                            ".csv.gz", sep = ""), "w")
        write.csv(entireData, gz1)
        close(gz1)
}
