################################################################################
## survivalTest_subMatrix.R
## : Select GO Terms that meets the conditions from the categorized score data 
##   of the corresponding species and run survival tests
##
## Additional Function Needed
## : survivalTest_filter.R, survivalTest_run.R
##
## Example of Arguments
## : survivalTest_subMatrix(args = list("name" = "HomoSapiens", "ontology" = "MF", 
##                               "matrix" = "BLOSUM62", "seed" =  1234))
##
## Writer: KyungJin Sohn
## Date: 02.25.2021
## Last modified: 07.13.2021
################################################################################
survivalTest_subMatrix <- function(args) {
                # 1. Install the needed library
        if(!require(dplyr)) 
                install.packages("dplyr",repos = "http://cran.us.r-project.org")
        library(dplyr)
        
        source("./survivalTest_filter.R")
        source("./survivalTest_run.R")

        
                # 2. Read categorized data that matches several criteria
        survivalData <- survivalTest_filter(args$name, args$ontology)
        
        
                # 3. Select GO Terms and sort by mean sumScore
        GOData <- survivalData %>% group_by(GO) %>% 
                summarise(avg_sumScore = mean(sumScore)) %>% 
                arrange(desc(avg_sumScore))
        

                # 4. Run the survival test and save the result
        if (!dir.exists(paste0(getwd(), "/survivalResult/", args$name, "/"))) {
                dir.create(paste0("./survivalResult/", args$name, "/"))
        }
        if (!dir.exists(paste0(getwd(), "/survivalResult/", args$name, "/", 
                               args$ontology, "/"))) {
                dir.create(paste0("./survivalResult/", args$name, "/",
                                  args$ontology, "/"))
        }
        
        for (i in 1:nrow(GOData)) {
                GOnum <- GOData[i,]$GO
                survivalData_GO <- survivalData %>% filter(GO == GOnum) %>% 
                        select(Names, GO, Seqs)
                
                survivedNum <- survivalTest_run(args$matrix, args$seed, survivalData_GO)
                
                gz1 <- gzfile(paste0("./survivalResult/", args$name, "/",
                                     args$ontology, "/", sub(":", "", GOnum), "_", 
                                     args$matrix, "_", args$seed, ".csv.gz"), "w")
                write.csv(survivedNum, gz1)
                close(gz1)
        }

}

