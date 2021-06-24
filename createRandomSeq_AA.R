################################################################################
## createRandomSeq_AA
## : Create random sequences with the same amino acid configuration
##
## Function Needed
## : 
##
## Example of Arguments
## : createRandomSeq_AA("Scerevisiae", 1111)
##
## Writer: KyungJin Sohn
## Date: 09.02.2020
## Last modified: 06.17.2021
################################################################################
createRandomSeq_AA <- function(species, seed) {
                # 1. Install the needed library
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
                install.packages("BiocManager")
        }
        if (!require("Biostrings", character.only = TRUE)) {
                BiocManager::install("Biostrings")
        }
        if(!require(stringr)) 
                install.packages("stringr",repos = "http://cran.us.r-project.org")
        if(!require(seqinr)) 
                install.packages("seqinr",repos = "http://cran.us.r-project.org")
        library(Biostrings)
        library(stringr)
        library(seqinr)
        
        
                # 2. Read the original data
        seqData <- readDNAStringSet(paste("./data/", species,".fasta.gz", sep = ""))
        
        seqData <- data.frame(Names = names(seqData), seqs = paste(seqData), stringsAsFactors = FALSE)
        
        
                # 3. Create file to save the random sequence data
        file_out = paste(species, "_randomAA_", seed, ".fasta", sep = "")
        file.create(file_out)
        
        
                # 4. set seed
        set.seed(seed)
        
        
                # 5. Create random sequences
        randCodonFun <- function(var) {
                # Find synonymous codon
                codon_list <- names(GENETIC_CODE[GENETIC_CODE %in% GENETIC_CODE[[var]]])
                
                # If there are more than one, choose from others
                if (length(codon_list) > 1) {
                        codon_list <- codon_list[!codon_list %in% var]
                }
                
                sample(codon_list, 1, replace = TRUE)
        }
        
        for (i in 1:nrow(seqData)) {
                # split DNA sequence to codon vector
                dna <- strsplit(seqData$seqs[i], "(?<=.{3})", perl = TRUE)[[1]]
                
                # if the length of the last codon is not 3, remove it.
                if (str_length(dna[length(dna)]) !=3 ) {
                        dna <- dna[-length(dna)]
                }
                
                # change codons to synonymous codon
                randCodon <- sapply(as.list(dna), randCodonFun)
                
                # save the result
                line_header <- paste("Random_", i, sep="")
                DNA_sequence <- paste(randCodon, collapse = "")
                
                write.fasta(DNA_sequence, 
                            line_header, 
                            file.out = file_out, 
                            open = "a", 
                            nbchar = 60, 
                            as.string = TRUE)
        }
        
        
                # 6. zip file
        gunzip(file_out, destname = paste0("./data/", file_out, ".gz"), 
               remove = TRUE)
}