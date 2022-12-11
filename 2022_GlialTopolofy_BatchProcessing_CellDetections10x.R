library("dplyr")
library("readr")

process_initial_data <- function(basePath, Cells_Path, filename, resuultsPath) {
  
  Cells_Raw <- read_tsv(paste0(basePath, "/", Cells_Path))
  
  # Convert to data frame
  Cells <- as.data.frame(Cells_Raw) 
  
  # Subset the date set to keep only relevant columns
  Cells <- subset(Cells, select = c(Image, Name, Class, `Centroid X µm`, `Centroid Y µm`))

  # Extract metadata information from image name
  Cells <- cbind(Cells, do.call(rbind , strsplit(Cells$Image , "[_\\.]"))[,1:3])
  colnames(Cells) <- c("Image", "ID", "Class", "X", "Y", "MouseID", "DPI", "Hemisphere")
  Cells <- subset(Cells, select = c(MouseID, Class, DPI, Hemisphere, ID, X, Y))
  
  # Write a .csv file 
  write.csv(Cells, paste0(resultsPath, "/", Cells_Path, filename))
}
basePath <- "D:/Daniel/Project-ECM/2.Images/Exp2-Gfap,NeuN,Iba1_Striatum_10X/QupathProjects"
resultsPath <- "D:/Daniel/Project-ECM/3.DataAnalysis/Exp2-Gfap,NeuN,Iba1_5x/Results/Coordinates_10x"

process_folder <- function(folderPath, filename_suffix) {
  files <- list.files(folderPath, pattern = "_detections.tsv", full.names = FALSE)
  for (file in files) {
    process_initial_data(folderPath, file, filename_suffix, resultsPath)
  }
}

process_folder(paste0(basePath, "/Iba1"), "_Coordinates.csv")
process_folder(paste0(basePath, "/Gfap"), "_Coordinates.csv")
process_folder(paste0(basePath, "/NeuN"), "_Coordinates.csv")
process_folder(paste0(basePath, "/Dapi"), "_Coordinates.csv")




