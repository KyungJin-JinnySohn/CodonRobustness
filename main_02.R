################################################################################
## main_02.R
## : The second main function which runs the survival test.
##
## Additional Function Needed
## : survivalTest_subMatrix.R
##
## Example of Arguments
## : .
##
## Writer: KyungJin Sohn
## Date: 02.25.2021
## Last modified: 07.07.2021
################################################################################

#!/home/kjsohn/.Renv/shims/Rscript

        # 1. Preparing libraries
if(!require(getopt)) 
        install.packages("getopt",repos = "http://cran.us.r-project.org")
if(!require(stringr)) 
        install.packages("stringr",repos = "http://cran.us.r-project.org")
library(getopt)
library(stringr)
source("./survivalTest_subMatrix.R")


        # 2. Option information
spec = matrix(c(
        'name', 'n', 1, "character", "Write the name of the species.(default: HomoSapiens)",
        'ontology' , 'o', 1, "character", "Write the ontology type of the data.(default: BP)",
        'matrix', 'm', 1, "character", "Write a name of a substitution matrix you want to use.(default: BLOSUM62)",
        'seed', 's', 1, "integer", "Write a seed number for a survival test.(default: 1234)",
        'help'     , 'h', 0, "logical", "Print this help message"
), byrow = TRUE, ncol = 5);
opt = getopt(spec)


        # 3. Conditions for each argument
species <- c("HomoSapiens", "MusMusculus", "Dmelanogaster", "Celegans", "Scerevisiae", "Ecoli")
if ( !is.null(opt$name) ) {
        if (!opt$name %in% species) {
                write("ERROR: The name of the species is unavailable.\n\t** These are the available species. **", stderr())
                print(species)
                q(status = 1)
        }
}

ONTO <- c("BP", "CC", "MF")
if ( !is.null(opt$ontology) ) {
        if (!opt$ontology %in% ONTO) {
                write("ERROR: The ontology type is unavailable.\n\t** These are the available ontologies. **", stderr())
                print(ONTO)
                q(status = 1)
        }
}

if ( !is.null(opt$matrix) ) {
        if ( ! file.exists(paste("./subMatrix/", opt$matrix, ".txt", sep = "")) ) {
                write("ERROR: You choose the substituation matrix that is not available.\n\t** These are the available substituation matrix. **", stderr())
                print(list.files("./subMatrix/"))
                q(status = 1)
        }
}

if ( !is.null(opt$help) ) {
        cat(getopt(spec, usage = TRUE))
        q(status = 1)
}

if ( is.null(opt$name) ) { opt$name = "HomoSapiens" }
if ( is.null(opt$ontology) ) { opt$ontology = "BP" }
if ( is.null(opt$matrix) ) { opt$matrix = "BLOSUM62" }
if ( is.null(opt$seed) ) { opt$seed = 1234 }


        # 4. Call the function that initiates the survival test
survivalTest_subMatrix(opt)
