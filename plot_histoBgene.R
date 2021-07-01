################################################################################
## plot_histoBgene
## : draw a histogram based on genes
##
## Function Needed
## : combine_random_data.R
##
## Example of Arguments
## : plot_histoBgene(scoreData, "Scerevisiae", "GC")
##
## Writer: KyungJin Sohn
## Date: 09.09.2020
## Last modified: 06.29.2021
################################################################################
plot_histoBgene <- function(scoreData, species, randomType) {
                # 1. Install the needed library
        if(!require(ggplot2))
                install.packages('ggplot2')
        if(!require(gridExtra))
                install.packages('gridExtra')
        if(!require(cowplot))
                install.packages('cowplot')
        library(ggplot2)
        library(gridExtra)
        library(grid)
        library(cowplot)
        source("./combine_random_data.R")
        
        
                # 2. Read score data of random sequences
        random_scoreData <- combine_random_data(species, randomType)
        
        
                # 3. Draw plot for each score
                        # color and description
        if (randomType == "GC") { 
                realC <- "#EFC000FF"
                description <- "* Radom: The GC-content is 100% matched."
        } else {
                realC <- "#0073C2FF"
                description <- "* Radom: The amino acids composition is 100% consistent."
        }
        
        scores <- c("wR2N", "CRI2", "eAI", "sumScore")
        breaks <- seq(from = 0, to = 1, by = 0.025)
        plots <- list() 
        
        for (i in scores) {
                selectedScore <- scoreData[, c("Names", i)]
                colnames(selectedScore) <- c("Names", "Score")
                
                randomScore <- random_scoreData[, c("Names", i)]
                colnames(randomScore) <- c("Names", "Score")
                
                        # t.test
                pvalue <- t.test(selectedScore$Score, randomScore$Score)$p.value
                if (pvalue != 0) {
                        if (pvalue >= 0.05) {
                                pvalue <- paste(formatC(pvalue, format="e", digits = 3), 
                                                "\n(>= 0.05)", sep = "")
                        } else {
                                pvalue <- formatC(pvalue, format="e", digits = 3)
                        }
                }
                pvalue <- paste("p-value = ", pvalue, sep = "")
                
                        # draw plot
                g <- ggplot() +
                        geom_histogram(data = randomScore,
                                       aes(Score, y = ..density.., fill = "Random"),
                                       breaks = breaks, color = "black", alpha = 0.5) +
                        geom_histogram(data = selectedScore,
                                       aes(Score, y = ..density.., fill = "Real  "),
                                       breaks = breaks, color = "black", alpha = 0.5) +
                        annotate("text", label = pvalue, hjust = 1, vjust = 1, 
                                 size = 5, x = 1, 
                                 y = max(density(selectedScore$Score)$y, density(randomScore$Score)$y)) +
                        xlab(i) +
                        theme(axis.text = element_text(size = 10), 
                              axis.title = element_text(size = 10)) +
                        theme_classic()
                
                plots[[length(plots) + 1]] <- g +
                        scale_fill_manual("Histogram Legend",
                                          values = c("Random" = "#868686FF", 
                                                     "Real  " = realC),
                                          guide = "none")
        }
        
                        # title
        title <- textGrob(paste("Histogram of ", species, " genes", sep = ""), 
                          gp = gpar(fontsize = 30))
        
                        # legend
        legend_hist <- g +
                scale_fill_manual(description,
                                  values = c("Random" = "#868686FF", 
                                             "Real  " = realC),
                                  guide = guide_legend(
                                          title.position = "bottom", 
                                          title.theme = element_text(size = 20),
                                          label.position = "top", 
                                          label.theme = element_text(size = 20),
                                          keywidth = 5, keyheight = 5, nrow = 1)) +
                theme(legend.justification = c("left", "top"), 
                      legend.margin = margin(0, 0, 0, 40))
        legend <- get_legend(legend_hist)
        
                        # arrange plot
        girdPlots <- arrangeGrob(plots[[1]], plots[[2]], plots[[3]], plots[[4]], 
                                 legend, ncol = 2, top = title, 
                                 padding = unit(5, "line"))
        
        
                # 3. Save the plot
        if (!dir.exists(paste0(getwd(), "/plot/histogram_all/", species, "/"))) {
                dir.create(paste0("./plot/histogram_all/", species, "/"))
        }
        ggsave(paste0(getwd(), "/plot/histogram_all/", species, "/", 
                      species, "_", randomType, "_histo.png"),
               width = 50, height = 40, units = "cm", girdPlots)
}