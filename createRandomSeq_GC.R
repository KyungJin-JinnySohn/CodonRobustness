################################################################################
## createRandomSeq_GC
## : Create random sequences with the same GC contents
##
## Function Needed
## : 
##
## Example of Arguments
## : createRandomSeq_GC("Scerevisiae", 1111)
##
## Writer: KyungJin Sohn
## Date: 08.24.2020
## Last modified: 06.17.2021
################################################################################
createRandomSeq_GC <- function(species, seed) {
                # 1. Install the needed library
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
                install.packages("BiocManager")
        }
        if (!require("Biostrings", character.only = TRUE)) {
                BiocManager::install("Biostrings")
        }
        if (!require("universalmotif", character.only = TRUE)) {
                BiocManager::install("universalmotif")
        }
        if(!require(stringr)) 
                install.packages("stringr",repos = "http://cran.us.r-project.org")
        if(!require(seqinr)) 
                install.packages("seqinr",repos = "http://cran.us.r-project.org")
        if(!require(R.utils)) 
                install.packages("R.utils")
        library(Biostrings)
        library(universalmotif)
        library(stringr)
        library(seqinr)
        library(R.utils)
        
        
                # 2. Read the original sequence data
        seqData <- readDNAStringSet(paste0("./data/", species,".fasta.gz"))
        
        
                # 3. Calculate A, C, G, T usage fraction of the original data
        bkg <- get_bkg(seqData, k = 1)
        DNA <- c ("A", "C", "G", "T")
        DNA_prob <- c("A" = round(bkg[which(bkg$klet == "A"), ]$probability, 2),
                      "C" = round(bkg[which(bkg$klet == "C"), ]$probability, 2),
                      "G "= round(bkg[which(bkg$klet == "G"), ]$probability, 2),
                      "T" = round(bkg[which(bkg$klet == "T"), ]$probability, 2))
        
        
                # 4. Prepare variables for the next step
        seqData <- data.frame(Names = names(seqData), seqs = paste(seqData), 
                              stringsAsFactors = FALSE)
        numSeq <- length(seqData$seqs)
        strlen <- str_length(seqData$seqs)
        
        
                # 5. Create file for the random sequence data
        file_out <- paste(species, "_randomGC_", seed, ".fasta", sep = "")
        file.create(file_out)
        
        
                # 6. set seed
        set.seed(seed)
        
        
                # 7. Create sequence and save the result
        for (i in 1:numSeq) {
                seq_length <- strlen[i]
                
                line_header <- paste("Random_", i, sep="")
                
                DNA_sequence <- paste(sample(DNA, seq_length, replace = TRUE, 
                                             prob = DNA_prob), 
                                      collapse = "")
                
                write.fasta(DNA_sequence, 
                            line_header, 
                            file.out = file_out, 
                            open = "a", 
                            nbchar = 60, 
                            as.string = TRUE)
        }
        
                # 8. zip file
        gunzip(file_out, destname = paste0("./data/", file_out, ".gz"), 
               remove = TRUE)
}

