################################################################################
## cateBonto_Celegans
## : Categorize Celegans score data based on ontology. [joined by: Names]
##
## Function Needed
## : .
##
## Example of Arguments
## : cateBonto_Celegans()
##
## Writer: KyungJin Sohn
## Date: 06.25.2020
## Last modified: 06.11.2021
################################################################################
cateBonto_Celegans <- function(species = "Celegans"){
        # 1. Install the needed library
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
                install.packages("BiocManager")
        }
        if (!require("Biostrings", character.only = TRUE)) {
                BiocManager::install("Biostrings")
        }
        if(!require(seqinr)) {
                install.packages('seqinr',repos = "http://cran.us.r-project.org")
        }
        if(!require(dplyr)) {
                install.packages('dplyr',repos = "http://cran.us.r-project.org")
        }
        library(Biostrings)
        library(seqinr)
        library(dplyr)
        
        
                # 2. Read the robustness score data
        scoreData <- gzfile(paste("./scoreResult/Score_", species,".csv.gz", 
                                  sep = ""))
        scoreData <- read.csv(scoreData, header = TRUE, row.names = 1, 
                              stringsAsFactors = FALSE)
        
                # Delete rows with NA values
        scoreData <- scoreData[complete.cases(scoreData), ]
        
                # Filter out the needed Names
        scoreData$Names <- data.frame(do.call('rbind', 
                                              strsplit(scoreData$Names, split = '[=]')),
                                      sstringsAsFactors = FALSE)[, 2]
        
                # 3. Read the GO data
        GOTXT <- read.table("./value/txt/Celegans_GO.gaf", sep = "\t", 
                            stringsAsFactors = FALSE)
        GOTXT <- GOTXT[, c(2, 5, 9)]
        GOTXT <- GOTXT[!duplicated(GOTXT),]
        colnames(GOTXT) <- c("Names", "GO", "ONTOLOGY")
        
        GOTXT[which(GOTXT$ONTOLOGY == 'P'), ]$ONTOLOGY <- "BP"
        GOTXT[which(GOTXT$ONTOLOGY == 'C'), ]$ONTOLOGY <- "CC"
        GOTXT[which(GOTXT$ONTOLOGY == 'F'), ]$ONTOLOGY <- "MF"
        
        
                # 4. Match GO data with scoreData
        DBdata <- inner_join(GOTXT, scoreData, by = "Names")
        
        
                # 5. Seperated DBdata by Ontology
        BP_DBdata <- DBdata %>% filter(ONTOLOGY == "BP")
        CC_DBdata <- DBdata %>% filter(ONTOLOGY == "CC")
        MF_DBdata <- DBdata %>% filter(ONTOLOGY == "MF")
        
        
                # 6. Save the data
        gz1 <- gzfile(paste("./scoreResult/", species, "_Score_BP.csv.gz", sep = ""), "w")
        write.csv(BP_DBdata, gz1)
        close(gz1)
        
        gz1 <- gzfile(paste("./scoreResult/", species, "_Score_CC.csv.gz", sep = ""), "w")
        write.csv(CC_DBdata, gz1)
        close(gz1)
        
        gz1 <- gzfile(paste("./scoreResult/", species, "_Score_MF.csv.gz", sep = ""), "w")
        write.csv(MF_DBdata, gz1)
        close(gz1)
}
