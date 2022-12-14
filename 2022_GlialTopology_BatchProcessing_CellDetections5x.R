library("dplyr")
library("readr")

process_initial_data <- function(Cells_Path, pattern, coordinates_columns, filename, brain_name) {
  
  Cells_Raw <- list.files(path=Cells_Path, pattern = pattern, full.names=TRUE) %>%        # Create object with all .tsv files in directory
    lapply(read_tsv) %>%                              # Store all files in list
    bind_rows                                         # Combine data sets into one data set 
  
  # Convert to data frame
  Cells <- as.data.frame(Cells_Raw) 
  
  # Transform coordinates from mm to microns for adecuate plotting in BrainRender
  Cells$Z <- (Cells[,c(coordinates_columns[[1]])]*1000)
  Cells$X <- (Cells[,c(coordinates_columns[[2]])]*1000)
  Cells$Y <- (Cells[,c(coordinates_columns[[3]])]*1000)
  
  # Subset the date set to keep only relevant columns
  Cells <- subset(Cells, select = c(Image, Name, Parent, Z, X, Y))
  Cells <- Cells[!(Cells$Name=="Negative" | Cells$Name=="Necrosis"),]
  Cells <- Cells %>% sample_frac(.1)
  
  # Extract metadata information from image name
  Cells <- cbind(Cells, do.call(rbind , strsplit(Cells$Image , "[_\\.]"))[,3:5])
  colnames(Cells) <- c( colnames(Cells[1:3]), paste0("New" , 1:3))
  Cells <- cbind(Cells[c(-2,-3)] , Cells[c(2,3)])
  
  # Rename columns
  colnames(Cells) <- c("Image", "Z", "X", "Y", "MouseID", "DPI", "Section", "ObjectID", "Region")
  
  Cells <- subset(Cells, select = c(MouseID, DPI, Region, Section, ObjectID, Z, X, Y))
  
  # Write a .csv file
  write.csv(Cells, filename)
  
  }

process_brain <- function(basePath, resultsPath, path) {
  
  set.seed(88071)
  
  # Load cells data set
  
  Iba1_Path <- paste0(basePath, "/", path, "/Iba1")
  Iba1_Filename <- paste0(resultsPath, "/", path, "_Iba1_Coordinates.csv")
  Iba1 <- process_initial_data(Iba1_Path, "detections.tsv", coordinates_columns = c(46, 47, 48), filename = Iba1_Filename, brain_name = path)
  
  
  Gfap_Path <- paste0(basePath, "/", path, "/Gfap")
  Gfap_Filename <- paste0(resultsPath, "/", path, "_Gfap_Coordinates.csv")
  process_initial_data(Gfap_Path, "detections.tsv", c(46, 47, 48), Gfap_Filename, path)
  
  
  Neun_Path <- paste0(basePath, "/", path, "/NeuN")
  Neun_Filename <- paste0(resultsPath, "/", path, "_NeuN_Coordinates.csv")
  process_initial_data(Neun_Path, "detections.tsv", c(46, 47, 48), Neun_Filename, path)
}

basePath <- "QupathProjects_5x"
resultsPath <- "ResultsTables/CellCoordinates_5x"

brains <- list.dirs(basePath, full.names = FALSE, recursive = FALSE)

for (brain in brains){
  process_brain(basePath, resultsPath, brain)
}


