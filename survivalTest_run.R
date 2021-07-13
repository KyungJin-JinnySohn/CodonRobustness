################################################################################
## survivalTest_run.R
## : Run survival tests
##
## Additional Function Needed
## : change_CToAA.R, cal_substitution_score.R
##
## Example of Arguments
## : .
##
## Writer: KyungJin Sohn
## Date: 02.25.2021
## Last modified: 07.07.2021
################################################################################
survivalTest_run <- function(matrix, seed, survivalData_GO) {
                # 1. Install the needed library
        if(!require(stringr)) 
                install.packages("stringr",repos = "http://cran.us.r-project.org")
        if(!require(dplyr)) 
                install.packages("dplyr",repos = "http://cran.us.r-project.org")
        library(stringr)
        library(dplyr)
        source("./change_CToAA.R")
        source("./cal_substitution_score.R")
        
        
                # 2. Read a substitution matrix
        subMatrix <- read.table(paste0("./subMatrix/", matrix, ".txt"), 
                                header = TRUE, sep = "", fill = TRUE)
        colnames(subMatrix) <- row.names(subMatrix)
        
        
                # 3. Add a column of AA(Amino Acid)
        survivalData_GO$AA <- ""
        for (i in 1:nrow(survivalData_GO)) {
                survivalData_GO$AA[i] <- change_CToAA(survivalData_GO$Seqs[i])
        }
        
        
                # 4. Add column of max-substitution score
        survivalData_GO$maxScore <- 0
        diagMatrix <- data.frame(Var1 = rownames(subMatrix), score = diag(as.matrix(subMatrix)))
        for (i in 1:nrow(survivalData_GO)) {
                maxAA <- table(strsplit(survivalData_GO$AA[i], "")[[1]])
                maxAA <- data.frame(maxAA)
                maxAA <- left_join(maxAA, diagMatrix, by = "Var1")
                
                survivalData_GO$maxScore[i] <- sum(maxAA$Freq * maxAA$score)
        }
        
        
                # 5. Run the survival test for each seq.
                        # set seed
        set.seed(args$seed)
        
        survivedNum <- ""
        for (i in 1:nrow(survivalData_GO)) {
                        # copy seq 100 times
                idx <- rep(i, 100) 
                mutationData_GO <- survivalData_GO[idx,]
                
                        # prepare values before running the survival test
                mutationData_GO$mutScore <- mutationData_GO$maxScore
                mutationData_GO$survived <- 1
                count <- 1
                
                        # run the survival test
                repeat {
                        for (i in which(mutationData_GO$survived == 1)) {
                                # select a random spot
                                randomSpot <- sample(1:str_length(mutationData_GO$Seqs[i]), 
                                                     1, replace = TRUE)
                                
                                # calculate substitution score after mutaion
                                calResult <- cal_substitution_score(mutationData_GO[i,], subMatrix, randomSpot)
                                mutationData_GO$mutScore[i] <- calResult$subScore
                                newC <- calResult$newC
                                
                                # replace codon into a new randomly chosen codon
                                substr(mutationData_GO$Seqs[i], randomSpot, randomSpot) <- newC
                        }
                        
                                # if the rate of difference is 0.8 or below, the sequence dies
                        if (sum(mutationData_GO$mutScore/mutationData_GO$maxScore < 0.8) > 0) {
                                mutationData_GO[mutationData_GO$mutScore/mutationData_GO$maxScore < 0.8, ]$survived <- 0
                        }
                        
                                
                                # save the result
                        survived <- data.frame(Names = mutationData_GO[1,]$Names,
                                               GO = mutationData_GO[1,]$GO,
                                               MutationNum = count,
                                               survivedNum = 
                                                       sum(mutationData_GO$survived == 1))
                        
                        survivedNum <- rbind(survivedNum, survived)
                        
                                # if all sequences die, stop the survival test
                        if (sum(mutationData_GO$survived == 1) == 0) { break }
                        
                        count <- count + 1
                }
        }
        survivedNum <- survivedNum[-1, ]
        
        return(survivedNum)
}

