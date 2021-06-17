################################################################################
## cateBonto_NucleotideID
## : Categorize HomoSapiens and MusMusculus score data based on ontology.
## [joined by: NucleotideID]
##
## Function Needed
## : .
##
## Example of Arguments
## : cateBonto_NucleotideID("HomoSapiens", "Homo Sapiens")
## : cateBonto_NucleotideID("MusMusculus", "musculus")
##
## Writer: KyungJin Sohn
## Date: 06.24.2020
## Last modified: 06.11.2021
################################################################################
cateBonto_NucleotideID <- function(species, key){
                # 1. Install the needed library
        if (!requireNamespace("BiocManager", quietly = TRUE)){
                install.packages("BiocManager")
        }
        if (!require("AnnotationHub", character.only = TRUE)) {
                BiocManager::install("AnnotationHub")
        }
        if(!require(dplyr)) {
                install.packages('dplyr',repos = "http://cran.us.r-project.org")
        }
        library(AnnotationHub)
        library(dplyr)
        
        
                # 2. Read the robustness score data
        scoreData <- gzfile(paste("./scoreResult/Score_", species,".csv.gz", 
                                  sep = ""))
        scoreData <- read.csv(scoreData, header = TRUE, row.names = 1, 
                              stringsAsFactors = FALSE)
        
                # Delete rows with NA values
        scoreData <- scoreData[complete.cases(scoreData), ]
        
                # Make ccds column
        names <- data.frame(do.call('rbind', 
                                    strsplit(scoreData$Names, split = '|', fixed = TRUE)), 
                            stringsAsFactors = FALSE)
        scoreData$Names <- names[, 1]
        
        
                # 3. Read (species)_NucleotideID.txt and join it with the score data
        IDTXT <- read.table(paste("./value/txt/", species, "_NucleotideID.txt", 
                                  sep = ""), 
                            sep = "\t", header = TRUE, stringsAsFactors = FALSE)
        IDTXT <- IDTXT %>% filter(source == "EBI") 
        IDTXT <- IDTXT[, c("ccds", "nucleotide_ID")]
        IDTXT <- IDTXT[!duplicated(IDTXT),]
        IDTXT$nucleotide_ID <- data.frame(do.call('rbind', 
                                                  strsplit(IDTXT$nucleotide_ID, split = '[.]')), 
                                          sstringsAsFactors = FALSE)[, 1]
        
                # Seperate IDTXT into one matching with only 1 nucleotide_id and others 
        IDTXT_Unique <- IDTXT %>% group_by(ccds) %>% mutate(IDfreq = n()) %>% 
                filter(IDfreq == 1) %>% arrange(ccds)
        IDTXT_Dup <- IDTXT %>% group_by(ccds) %>% mutate(IDfreq = n()) %>% 
                filter(IDfreq > 1) %>% arrange(ccds)
        
        
                # 4. Download AnnotationHub and request needed database
        ### CAUTION: Need to find the exact name of the species in OrgDb ###
        hub <- AnnotationHub()
        sub_hub <- query(hub, c("OrgDb", key))
        orgdb <- query(sub_hub, "OrgDb")[[1]]
        
        
                # 5. Match GO data with scoreData
                # 5-1. Unique
        transids <- as.vector(IDTXT_Unique$nucleotide_ID)
        DBdata <- select(orgdb, keys = transids, columns = "GO", 
                         keytype="ENSEMBLTRANS")
        
                # Select ORF, GO, ONTOLOGY columns
        DBdata <- DBdata[, c("ENSEMBLTRANS", "GO", "ONTOLOGY")]
        colnames(DBdata) <- c("nucleotide_ID", "GO", "ONTOLOGY")
        
                # Delete rows which have NA values
        DBdata <- DBdata[complete.cases(DBdata), ]
        
                # Delete duplicated rows
        DBdata <- DBdata[!duplicated(DBdata), ]
        
        
                # 6. Seperated DBdata by Ontology
        BP_DBdata <- DBdata %>% filter(ONTOLOGY == "BP")
        CC_DBdata <- DBdata %>% filter(ONTOLOGY == "CC")
        MF_DBdata <- DBdata %>% filter(ONTOLOGY == "MF")
        
        
                # 7. Joing 00_DBdata with score data
        BP_DBdata <- inner_join(BP_DBdata, IDTXT_Unique, by = "nucleotide_ID")
        CC_DBdata <- inner_join(CC_DBdata, IDTXT_Unique, by = "nucleotide_ID")
        MF_DBdata <- inner_join(MF_DBdata, IDTXT_Unique, by = "nucleotide_ID")
        
        totalBP_DBdata <- BP_DBdata[, c("ccds", "GO", "ONTOLOGY")]
        totalCC_DBdata <- CC_DBdata[, c("ccds", "GO", "ONTOLOGY")]
        totalMF_DBdata <- MF_DBdata[, c("ccds", "GO", "ONTOLOGY")]
        
                # 5-2. Duplicated
        transids <- as.vector(IDTXT_Dup$nucleotide_ID)
        DBdata <- select(orgdb, keys = transids, columns = "GO", 
                         keytype="ENSEMBLTRANS")
        
                # Select ORF, GO, ONTOLOGY columns
        DBdata <- DBdata[, c("ENSEMBLTRANS", "GO", "ONTOLOGY")]
        colnames(DBdata) <- c("nucleotide_ID", "GO", "ONTOLOGY")
        
                # Delete rows which have NA values
        DBdata <- DBdata[complete.cases(DBdata), ]
        
                # Delete duplicated rows
        DBdata <- DBdata[!duplicated(DBdata), ]
        
        
                # 6. Seperated DBdata by Ontology
        BP_DBdata <- DBdata %>% filter(ONTOLOGY == "BP")
        CC_DBdata <- DBdata %>% filter(ONTOLOGY == "CC")
        MF_DBdata <- DBdata %>% filter(ONTOLOGY == "MF")
        
        
                # 7. Joing 00_DBdata with score data
        BP_DBdata <- inner_join(BP_DBdata, IDTXT_Dup, by = "nucleotide_ID")
        BP_DBdata <- BP_DBdata %>% group_by(ccds, GO) %>% mutate(GOfreq = n()) %>% filter(GOfreq == IDfreq)
        BP_DBdata <- BP_DBdata[, c("ccds", "GO", "ONTOLOGY")]
        BP_DBdata <- BP_DBdata[!duplicated(BP_DBdata), ]
        
        CC_DBdata <- inner_join(CC_DBdata, IDTXT_Dup, by = "nucleotide_ID")
        CC_DBdata <- CC_DBdata %>% group_by(ccds, GO) %>% mutate(GOfreq = n()) %>% filter(GOfreq == IDfreq)
        CC_DBdata <- CC_DBdata[, c("ccds", "GO", "ONTOLOGY")]
        CC_DBdata <- CC_DBdata[!duplicated(CC_DBdata), ]
        
        MF_DBdata <- inner_join(MF_DBdata, IDTXT_Dup, by = "nucleotide_ID")
        MF_DBdata <- MF_DBdata %>% group_by(ccds, GO) %>% mutate(GOfreq = n()) %>% filter(GOfreq == IDfreq)
        MF_DBdata <- MF_DBdata[, c("ccds", "GO", "ONTOLOGY")]
        MF_DBdata <- MF_DBdata[!duplicated(MF_DBdata), ]
        
                # 8. append each DBdata
        totalBP_DBdata <- rbind(totalBP_DBdata, as.data.frame(BP_DBdata))
        totalCC_DBdata <- rbind(totalCC_DBdata, as.data.frame(CC_DBdata))
        totalMF_DBdata <- rbind(totalMF_DBdata, as.data.frame(MF_DBdata))
        colnames(totalBP_DBdata) <- c("Names", "GO", "ONTOLOGY")
        colnames(totalCC_DBdata) <- c("Names", "GO", "ONTOLOGY")
        colnames(totalMF_DBdata) <- c("Names", "GO", "ONTOLOGY")
        
                # 9. Joing 00_DBdata with score data
        totalBP_DBdata <- inner_join(totalBP_DBdata, scoreData, by = "Names")
        totalCC_DBdata <- inner_join(totalCC_DBdata, scoreData, by = "Names")
        totalMF_DBdata <- inner_join(totalMF_DBdata, scoreData, by = "Names")
        
        
                # 10. Save the data
        gz1 <- gzfile(paste("./scoreResult/", species, "_Score_BP.csv.gz", sep = ""), "w")
        write.csv(totalBP_DBdata, gz1)
        close(gz1)
        
        gz1 <- gzfile(paste("./scoreResult/", species, "_Score_CC.csv.gz", sep = ""), "w")
        write.csv(totalCC_DBdata, gz1)
        close(gz1)
        
        gz1 <- gzfile(paste("./scoreResult/", species, "_Score_MF.csv.gz", sep = ""), "w")
        write.csv(totalMF_DBdata, gz1)
        close(gz1)
        
}
