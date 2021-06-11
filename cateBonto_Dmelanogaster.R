################################################################################
## cateBonto_Dmelanogaster
## : Categorize Dmelanogaster score data based on ontology. [joined by: Names]
##
## Function Needed
## : .
##
## Example of Arguments
## : cateBonto_Dmelanogaster()
##
## Writer: KyungJin Sohn
## Date: 06.25.2020
## Last modified: 06.11.2021
################################################################################
cateBonto_Dmelanogaster <- function(species = "Dmelanogaster") {
                # 1. Install the needed library
        if (!requireNamespace("BiocManager", quietly = TRUE)){
                install.packages("BiocManager")
                BiocManager::install("AnnotationHub")
        }
        if(!require(dplyr)) 
                install.packages('dplyr',repos = "http://cran.us.r-project.org")
        library(AnnotationHub)
        library(dplyr)
        
        
                # 2. Read the robustness score data
        scoreData <- gzfile(paste("./scoreResult/Score_", species,".csv.gz", 
                                  sep = ""))
        scoreData <- read.csv(scoreData, header = TRUE, row.names = 1, 
                              stringsAsFactors = FALSE)
        
                # Delete rows with NA values
        scoreData <- scoreData[complete.cases(scoreData), ]
        
                # Filter out the needed Names
        for (i in 1:nrow(scoreData)) {
                index <- grep("FlyBase:", strsplit(scoreData$Names[i], "[;]")[[1]])
                splitName <- strsplit(scoreData$Names[i], "[;]")[[1]][index]
                
                index <- grep("FlyBase:", strsplit(splitName, "[,]")[[1]])
                splitName <- strsplit(splitName, "[,]")[[1]][index]
                
                splitName <- strsplit(splitName, "[:]")[[1]][2]
                scoreData$Names[i] <- splitName
        }
        
        
                # 3. Download AnnotationHub and request needed database
        ### CAUTION: Need to find the exact name of the species in OrgDb ###
        hub <- AnnotationHub()
        sub_hub <- query(hub, c("OrgDb", "melanogaster"))
        orgdb <- query(sub_hub, "OrgDb")[[1]]
        
        flyids <- as.vector(scoreData$Names)
        DBdata <- select(orgdb, keys = flyids, columns = "GO", keytype="FLYBASEPROT")
        
                # Select FLYBASEPROT, GO, ONTOLOGY columns
        DBdata <- DBdata[, c("FLYBASEPROT", "GO", "ONTOLOGY")]
        colnames(DBdata) <- c("Names", "GO", "ONTOLOGY")
        
                # Delete rows which have NA values
        DBdata <- DBdata[complete.cases(DBdata), ]
        
                # Delete duplicated rows
        DBdata <- DBdata[!duplicated(DBdata), ]
        
        
                # 5. Seperated DBdata by Ontology
        BP_DBdata <- DBdata %>% filter(ONTOLOGY == "BP")
        CC_DBdata <- DBdata %>% filter(ONTOLOGY == "CC")
        MF_DBdata <- DBdata %>% filter(ONTOLOGY == "MF")
        
        
                # 6. Joing 00_DBdata with score data
        BP_DBdata <- inner_join(BP_DBdata, scoreData, by = "Names")
        CC_DBdata <- inner_join(CC_DBdata, scoreData, by = "Names")
        MF_DBdata <- inner_join(MF_DBdata, scoreData, by = "Names")
        
        
                # 7. Save the data
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
