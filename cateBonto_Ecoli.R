################################################################################
## cateBonto_Ecoli
## : Categorize Ecoli score data based on ontology. [joined by: Names]
##
## Function Needed
## : .
##
## Example of Arguments
## : cateBonto_Ecoli()
##
## Writer: KyungJin Sohn
## Date: 08.24.2020
## Last modified: 06.11.2021
################################################################################
cateBonto_Ecoli <- function(species = "Ecoli") {
                # 1. Install the needed library
        if(!require(dplyr)) 
                install.packages('dplyr',repos = "http://cran.us.r-project.org")
        library(dplyr)
        
        
                # 2. Read the robustness score data
        scoreData <- gzfile(paste("./scoreResult/Score_", species,".csv.gz", 
                                  sep = ""))
        scoreData <- read.csv(scoreData, header = TRUE, row.names = 1, 
                              stringsAsFactors = FALSE)
        
                # Delete rows with NA values
        scoreData <- scoreData[complete.cases(scoreData), ]
        
                # Select only gene name from the 'Names' column
        scoreData$Names <- sub(".*gene_symbol:", "", scoreData$Names)
        scoreData$Names <- sub(".desc..*", "", scoreData$Names)
        
                # Delete rows those do not have a gene name
        scoreData <- scoreData[nchar(as.character(scoreData$Names)) <= 10, ]
        
        
                # 3. Read the GO data
        GOTXT <- gzfile("./value/txt/ecocyc.gaf.gz")
        GOTXT <- read.table(GOTXT, fill=TRUE, quote="", skip = 19,
                            sep="\t", stringsAsFactors = FALSE)
        GOTXT <- GOTXT[, c(3, 5, 9)]
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
