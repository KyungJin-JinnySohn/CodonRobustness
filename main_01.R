################################################################################
## main_01.R
## : The first main function to get the inputs needed to start 
##      robustness score calculations.
##
## Additional Function Needed
## : cal_robustness_score.R
##
## Example of Arguments
## : .
##
## Writer: KyungJin Sohn
## Date: 06.01.2020
## Last modified: 06.07.2021
################################################################################

#!/home/kjsohn/.Renv/shims/Rscript

        # 1. Preparing libraries
if(!require(getopt)) 
        install.packages("getopt",repos = "http://cran.us.r-project.org")
if(!require(stringr)) 
        install.packages("stringr",repos = "http://cran.us.r-project.org")

library(getopt)
library(stringr)

        # 2. Option information
spec = matrix(c(
        'inputData', 'i', 1, "character", "Write a file name of the input data.",
        'outputData', 'o', 1, "character", "Write a file name to save the results.(default: Same as input data file name.)",
        'species', 's', 1, "character", "Write a species of the input data.",
        'help', 'h', 0, "logical", "Print this help message"
), byrow = TRUE, ncol = 5);
opt = getopt(spec)

        # 3. Conditions for each argument
if ( !is.null(opt$help) ) {
        cat(getopt(spec, usage = TRUE))
        q(status = 1)
}

if (!is.null(opt$inputData)) {
        if (!file.exists(paste("./data/", opt$inputData, sep = ""))) {
                write("ERROR: The input file does not exist.
                ** These are the available input file. **", stderr())
                print(list.files("./data/"))
                q(status = 1)
        }
}

if (is.null(opt$inputData)) { 
        write("ERROR: Please write a file name of the input data.
              ** These are the available input file. **", stderr())
        print(list.files("./data/"))
        q(status = 1) 
}

if (is.null(opt$outputData)) { 
        opt$outputData = strsplit(opt$inputData, "[.]")[[1]][1]
}

if (!is.null(opt$species)) {
        if (!file.exists(paste("./value/RTA/RTA_", 
                                 opt$species, ".csv", sep = "")) ) {
                write("ERROR: The selected species is not available.
                ** These are the available species. **", stderr())
                speciesFile <- list.files("./value/RTA/")
                species <- sapply(speciesFile, function(fileName) { 
                        name <- strsplit(fileName, "[_]")[[1]][2]
                        strsplit(name, "[.]")[[1]][1]
                })
                print(as.vector(species))
                q(status = 1)
        }
}

if (is.null(opt$species)) { 
        write("ERROR: Please write a species of the input data.
              ** These are the available species. **", stderr())
        speciesFile <- list.files("./value/RTA/")
        species <- sapply(speciesFile, function(fileName) { 
                name <- strsplit(fileName, "[_]")[[1]][2]
                strsplit(name, "[.]")[[1]][1]
        })
        print(as.vector(species))
        q(status = 1) 
}

        # 4. Call the function that initiates the score calculation
source("./cal_robustness_score.R")
cal_robustness_score(opt)
