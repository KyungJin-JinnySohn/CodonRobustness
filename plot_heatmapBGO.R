################################################################################
## plot_heatmapBGO
## : draw a histogram based on genes for multiple scores at once
##
## Function Needed
## : .
##
## Example of Arguments
## : plot_heatmapBGO(c("Scerevisiae", "Celegans", "Ecoli"), "BP", "wR2N",
##                   scoreData, "rank")
##
## Writer: KyungJin Sohn
## Date: 09.07.2020
## Last modified: 07.06.2021
################################################################################
plot_heatmapBGO <- function(species, ONTO, score, scoreData, indexType) {
                # 1. Install the needed library
        if(!require(RColorBrewer)){
                install.packages('RColorBrewer',repos = "http://cran.us.r-project.org")
        }
        if(!require(pheatmap)){
                install.packages('pheatmap',repos = "http://cran.us.r-project.org")
        }
        library(RColorBrewer)
        library(pheatmap)
        
        
                # 2. Create folder to save the heatmap
        filePath <- paste0(getwd(), "/plot/heatmap/", indexType, "/",
                           paste0(substring(species, 1, 1), collapse = "_"))
        if (!dir.exists(filePath)) {
                dir.create(filePath)
        }

        fileName <- paste0(filePath, "/", 
                           paste(substring(species, 1, 1), collapse = ""), 
                           "_", ONTO, "_", score, "_heatmap.pdf")   
        
        
                # 3. Draw heatmap for 
        if (indexType == "rank") {
                scoreData <- scoreData %>% 
                        mutate(Scerevisiae = dense_rank(-Scerevisiae),
                               Celegans = dense_rank(-Celegans),
                               Ecoli = dense_rank(-Ecoli))
                legend_color <- colorRampPalette(
                        rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
                numFormat <- "%d"
        } else {
                legend_color <- colorRampPalette(
                        brewer.pal(n = 7, name = "RdYlBu"))(100)
                numFormat <- "%.2f"
        }
        
        pheatmap(scoreData, scale = "none", cluster_cols = TRUE, cluster_row = TRUE,  
                 cellwidth = 10, cellheight = 2, border_color = NA,   
                 main = paste0("\nHeatmap of ", score, "(", ONTO, ", ", 
                               nrow(scoreData), "GOs)"),
                 fontsize = 2, display_numbers = TRUE, number_format = numFormat,
                 legend = TRUE, color = legend_color,
                 filename = fileName)
}