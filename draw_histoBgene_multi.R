################################################################################
## draw_histoBgene_multi
## : draw a histogram based on genes and draw multiple species at once.
## 
## Function Needed
## : plot_histoBgene_multi.R
##
## Example of Arguments
## : draw_histoBgene_multi(
##      c("HomoSapiens", "MusMusculus", "Dmelanogaster", "Celegans", "Scerevisiae", "Ecoli"), 
##      c("wR2N", "CRI2", "eAI", "sumScore"), "AA") # of "GC"
##
## Writer: KyungJin Sohn
## Date: 09.09.2020
## Last modified: 07.01.2021
################################################################################
draw_histoBgene_multi <- function(species, scores, randomType) {
                # 1. Throw an error if the length of the species vector is less than 2
        if (length(species) < 2) {
                stop("ERROR: Please select more than one species.")
        }
        
        
                # 2. Install the needed library
        source("./plot_histoBgene_multi.R")
        
        
                # 3. Draw histogram for each species
        allNames <- list()
        allplots <- list()
        
        for (speciesName in species) {
                temp <- plot_histoBgene_multi(speciesName, scores, randomType)
                allNames <- append(allNames, temp[["names"]])
                allplots <- append(allplots, temp[["plots"]])
        }
        
                        # title
        title <- textGrob("Histogram of Scores Distributions for diverse species", 
                          gp = gpar(fontsize = 30))
        
                        # legend
        legend <- temp[["legend"]]
        
                        # arrange plot and legend
        girdNames <- arrangeGrob(grobs = allNames, ncol = 1, 
                                 layout_matrix = matrix(1:length(species), 
                                                        ncol = 1))
        
        girdPlots <- arrangeGrob(grobs = allplots, ncol = length(scores), 
                                 layout_matrix = matrix(1:length(allplots), 
                                                        ncol = length(scores), 
                                                        byrow = TRUE))
        
        girdAll <- arrangeGrob(girdNames, girdPlots, legend, 
                               ncol = 1 + 2 * length(scores), 
                               top = title, padding = unit(5, "line"),
                               layout_matrix = rbind(cbind(rep(1, length(species)), 
                                                           matrix(2, ncol = 2 * length(scores), 
                                                                  nrow = length(species))),
                                                     c(NA, rep(3, 2 * length(scores)))))
        
        
                # 4. Save the plot
        if (!dir.exists(paste0(getwd(), "/plot/histogram_all/"))) {
                dir.create(paste0("./plot/histogram_all/"))
        }
        imgHegith <- 15 + 5 * length(species)
        ggsave(paste0(getwd(), "/plot/histogram_all/", 
                      paste0(substring(species, 1, 1), collapse = ""), "_", 
                      paste0(substring(scores, 1, 1), collapse = ""), "_",
                      randomType, "_histo.png"),
               width = 60, height = imgHegith, units = "cm", girdAll)
}

