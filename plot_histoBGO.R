################################################################################
## plot_histoBGO
## : drawing histogram according to the group of selected GO group
##
## Function Needed
## : subNmodi_score.R
##
## Example of Arguments
## : plot_histoBGO("Scerevisiae", vector, df, char_string)
##
## Writer: KyungJin Sohn
## Date: 06.24.2020
## Last modified: 06.22.2021
################################################################################
plot_histoBGO <- function(species, selectedGO, scoreDatatGO, groupName) {
                # 1. Install the needed library
        if(!require(dplyr))
                install.packages('dplyr',repos = "http://cran.us.r-project.org")
        if(!require(ggplot2))
                install.packages('ggplot2')
        if(!require(gridExtra))
                install.packages('gridExtra')
        if(!require(cowplot))
                install.packages('cowplot')
        library(dplyr)
        library(ggplot2)
        library(grid)
        library(gridExtra)
        library(cowplot)
        source("./subNmodi_score.R")
        
        
                # 2. Read score data of random sequence
        seed <- 1111 # change when there are more than one random data
        randomType <- "GC" # change if you want to use diff. type of random data
        random_scoreData <- subNmodi_score(
                paste0("random/Score_", species, "_random", randomType, "_", seed))
        
        
                # 3. Filter out the selected GO term score_data
        selectedGene <- scoreDatatGO %>% filter(GO %in% selectedGO)
        
        
                # 4. Draw plot for each score
        scores <- c("wR2N", "CRI2", "eAI", "sumScore")
        breaks <- seq(from = 0, to = 1, by = 0.025)
        plots <- list() 
        
        for (i in scores) {
                selectedScore <- selectedGene[, c("Names", i)]
                colnames(selectedScore) <- c("Names", "Score")
                
                randomScore <- random_scoreData[, c("Names", i)]
                colnames(randomScore) <- c("Names", "Score")
                
                totalScore <- scoreDatatGO[, c("Names", i)]
                colnames(totalScore) <- c("Names", "Score")
                
                g <- ggplot() +
                        geom_histogram(data = randomScore,
                                       aes(Score, y = ..density.., fill = "Random"),
                                       breaks = breaks, color = "lightgray", alpha = 1) +
                        geom_histogram(data = totalScore,
                                       aes(Score, y = ..density.., fill = "Total"),
                                       breaks = breaks, color = "darkgray", alpha = 0.7) +
                        geom_histogram(data = selectedScore,
                                       aes(Score, y = ..density.., fill = "Selected"),
                                       breaks = breaks, color = "black", alpha = 0.5) +
                        labs (x = i) +
                        theme(axis.text = element_text(size = 10),
                              axis.title = element_text(size = 10)) +
                        theme_classic()
                
                plots[[length(plots) + 1]] <- g +
                        scale_fill_manual("Histogram Legend",
                                          values = c("Random" = "white", 
                                                     "Total" = "lightgray", 
                                                     "Selected" = "black"),
                                          guide = "none")
        }
        
        title <- textGrob(paste0("Histogram of genes containg ", groupName), 
                          gp = gpar(fontsize = 30))
        
        legend_hist <- g +
                scale_fill_manual("Histogram Legend",
                                  values = c("Random" = "white", 
                                             "Total" = "lightgray", 
                                             "Selected" = "black"),
                                  guide = guide_legend(
                                          title.position = "top",
                                          title.theme = element_text(size = 25),
                                          label.position = "bottom",
                                          label.theme = element_text(size = 20),
                                          keywidth = 5,
                                          keyheight = 5,
                                          nrow = 1))
        legend <- get_legend(legend_hist)
        
        girdPlots <- arrangeGrob(plots[[1]], plots[[2]], plots[[3]], plots[[4]],
                                 legend, nrow = 3, top = title, 
                                 padding = unit(5, "line"))
        
        
                # 5. Save the plot
        if (!dir.exists(paste0(getwd(), "/plot/histogram/", species, "/"))) {
                dir.create(paste0("./plot/histogram/", species, "/"))
        }
        ggsave(paste0(getwd(), "/plot/histogram/", species, "/", groupName, "_GOHisto.png"),
               width = 50, height = 40, units = "cm", girdPlots)
}

