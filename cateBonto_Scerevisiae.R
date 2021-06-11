################################################################################
## cateBonto_Scerevisiae
## : Categorize Scerevisiae score data based on ontology. [joined by: Names]
##
## Function Needed
## : .
##
## Example of Arguments
## : .
##
## Writer: KyungJin Sohn
## Date: 06.24.2020
## Last modified: 06.11.2021
################################################################################
cateBonto_Scerevisiae <- function(species = "Scerevisiae") {
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
        scoreData$Names <- data.frame(do.call('rbind', 
                                              strsplit(scoreData$Names, split = ' ')),
                                      stringsAsFactors = FALSE)[, 1]
        
        
                # 3. Download AnnotationHub and request needed database
        ### CAUTION: Need to find the exact name of the species in OrgDb ###
        hub <- AnnotationHub()
        sub_hub <- query(hub, c("OrgDb", "cerevisiae")) ## CAUTION: Change species name when needed
        orgdb <- query(sub_hub, "OrgDb")[[1]]
        #keytypes(orgdb)
        
        
                # 4. Match GO data with scoreData
        orfids <- as.vector(scoreData$Names)
        DBdata <- select(orgdb, keys = orfids, columns = "GO", keytype="ORF")
        
                # Select ORF, GO, ONTOLOGY columns
        DBdata <- DBdata[, c("ORF", "GO", "ONTOLOGY")]
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

