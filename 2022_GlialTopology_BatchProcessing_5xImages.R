set.seed(88071)
library(spatstat)
library(ggplot2)
library(dplyr)
library(purrr)

setwd("D:/Research/Project_GlialTopology/3.DataAnalysis/Exp2-Gfap,NeuN,Iba1")

coordinatesPath <- "CellCoordinates_5x"
densityTablesPath <- "ResultsTables"

Cells_Intensity_CSV_Path <- paste0(densityTablesPath, "/Cells_Intensity.csv")
Cells_Intensity_Header <- c("Brain", "Neurons_Intensity", "Astrocytes_Intensity", "Microglia_Intensity")

Tesselation_CSV_Path <- paste0(densityTablesPath, "/Cells_Covariance.csv")
Tesselation_Test_Header <- c("Brain", "AN1", "AN2", "AN3", "MN1", "MN2", "MN3", "AM1", "AM2", "AM3")

# nnndist tables
Neurons_nndist_CSV_Path <- paste0(densityTablesPath, "/Neurons_nndist.csv")
Neurons_nndist_Header <- c("x", "y", "nndist", "Brain")

Astrocytes_nndist_CSV_Path <- paste0(densityTablesPath, "/Astrocytes_nndist.csv")
Astrocytes_nndist_Header <- c("x", "y", "nndist", "Brain")

Microglia_nndist_CSV_Path <- paste0(densityTablesPath, "/Microglia_nndist.csv")
Microglia_nndist_Header <- c("x", "y", "nndist", "Brain")

# Results to generate
Result_Hyperframe <- NULL

# Functions

add_to_hyperframe <- function (...) {
    if (is.null(Result_Hyperframe)){
      Result_Hyperframe <<- hyperframe(...)
    } else {
      Result_Hyperframe <<- rbind(Result_Hyperframe, hyperframe(...))
    }
}


create_empty_table <- function (path, header) {
  df_header <- data.frame(matrix(ncol = length(header), nrow = 0))
  names(df_header) <- header

  write.csv(df_header, path)
}

create_empty_table(Cells_Intensity_CSV_Path, Cells_Intensity_Header)
create_empty_table(Tesselation_CSV_Path, Tesselation_Test_Header)
create_empty_table(Neurons_nndist_CSV_Path, Neurons_nndist_Header)
create_empty_table(Astrocytes_nndist_CSV_Path, Astrocytes_nndist_Header)
create_empty_table(Microglia_nndist_CSV_Path, Microglia_nndist_Header)

# Manipulate coordinates for correct plotting in R
coordinates_manipulation <- function (Raw_Table) {
  Cell_Coor_X <- Raw_Table$Y
  Cell_Coor_Y <- Raw_Table$X

  ## Bind the vectors, rotate and bind to original table
  Coords <- cbind(Cell_Coor_X, Cell_Coor_Y)
  Coords <- secr::rotate(Coords, 180)
  Coords <- as.data.frame(Coords)
  return(cbind(Raw_Table, Coords))
}

# Cretate a point pattern (PPP) object

create_point_pattern <- function(Subset) {
  # We define the limits of the window according to Neuron coordinates
  xlim <- range(Subset$Cell_Coor_X)
  ylim <- range(Subset$Cell_Coor_Y)

  # Create point pattern for neurons
  Cells_PPP <- with(Subset, ppp(x = Subset$Cell_Coor_X, y = Subset$Cell_Coor_Y, xrange = xlim, yrange = ylim))
  unitname(Cells_PPP)  <- list("mm", "mm", 1.3/1000)
  Cells_PPP <- spatstat.geom::rescale (Cells_PPP)
  
  ## We rescale the unit to obtain measurements in mm2
  return(Cells_PPP)
}

define_convex_hull <- function(Neurons_PPP, Cells_PPP) {
  chull <- convexhull(Neurons_PPP)
  Window(Cells_PPP) <- chull
  return(Cells_PPP)
}


tesselation <- function(Cells_Density) {
  ## We define the quantiles for Neurons
  Cells_Quantiles <- c(0, 50, 100, 300)

  ## We define the cutting spots according to quantiles
  Cells_Cut <- cut(Cells_Density, breaks = Cells_Quantiles, labels = 1:3)

  ## We generate the tesselation image
  return(tess(image = Cells_Cut))
}

tesselation_data <- function(Cells_PPP, Cells_Tess) {
  Result <- quadratcount(Cells_PPP, tess = Cells_Tess )
  return(Result)
}

Neurons_Astrocytes_Function_Vector <- c()
Microglia_Neurons_Function_Vector <- c()
Astrocytes_Microglia_Function_Vector <- c()


process_file <- function (basePath, path) {

  Neurons_Raw <- read.csv(file = paste0(basePath, '/', path, '_NeuN_Coordinates.csv'), header = TRUE)
  Astrocytes_Raw <- read.csv(file = paste0(basePath, '/', path, '_Gfap_Coordinates.csv'), header = TRUE)
  Microglia_Raw <- read.csv(file = paste0(basePath, '/', path, '_Iba1_Coordinates.csv'), header = TRUE)

  Neurons <- coordinates_manipulation(Neurons_Raw)
  Astrocytes <- coordinates_manipulation(Astrocytes_Raw)
  Microglia <- coordinates_manipulation(Microglia_Raw)

  # Subset neurons
  Neurons_Subset <- Neurons[(Neurons$Section=="Scene3"),]
  Neurons_Subset <- Neurons_Subset[(Neurons_Subset$Y < 5000),]

  # We subset astrocytes
  Astrocytes_Subset <- Astrocytes[(Astrocytes$Section=="Scene3"),]
  Astrocytes_Subset <- Astrocytes_Subset[(Astrocytes_Subset$Y < 5000),]

  # We subset microglia
  Microglia_Subset <- Microglia[(Microglia$Section=="Scene3"),]
  Microglia_Subset <- Microglia_Subset [(Microglia_Subset$Y < 5000),]

  Neurons_PPP <- create_point_pattern(Neurons_Subset)
  Astrocytes_PPP <- create_point_pattern(Astrocytes_Subset)
  Microglia_PPP <- create_point_pattern(Microglia_Subset)

  Neurons_PPP <- define_convex_hull(Neurons_PPP, Neurons_PPP)
  Astrocytes_PPP <- define_convex_hull(Neurons_PPP, Astrocytes_PPP)
  Microglia_PPP <- define_convex_hull(Neurons_PPP, Microglia_PPP)

  Neurons_Intensity <- summary(Neurons_PPP)$intensity
  Astrocytes_Intensity <- summary(Astrocytes_PPP)$intensity
  Microglia_Intensity <- summary(Microglia_PPP)$intensity

  Intensity_Row <- t(c(path, Neurons_Intensity, Astrocytes_Intensity, Microglia_Intensity))
  write.table(Intensity_Row, Cells_Intensity_CSV_Path, append = TRUE, sep=",", col.names = FALSE)

  Microglia_Density <- density(Microglia_PPP, sigma =0.1, positive=TRUE, equal.ribbon = TRUE, col = topo.colors, main = "")
  Astrocytes_Density <- density(Astrocytes_PPP, sigma =0.1, positive=TRUE, equal.ribbon = TRUE, col = topo.colors, main = "")
  Neurons_Density <- density(Neurons_PPP, sigma =0.1, positive=TRUE, equal.ribbon = TRUE, col = topo.colors, main = "")

  Neurons_Tess <- tesselation(Neurons_Density)
  Astrocytes_Tess <- tesselation(Astrocytes_Density)
  Microglia_Tess <- tesselation(Microglia_Density)

  Astrocytes_in_Neurons <-tesselation_data(Astrocytes_PPP, Neurons_Tess)
  Microglia_in_Neurons <- tesselation_data(Microglia_PPP, Neurons_Tess)
  Astrocytes_in_Microglia <- tesselation_data(Astrocytes_PPP, Microglia_Tess)

  Tesselation_Row <- t(c(path, Astrocytes_in_Neurons, Microglia_in_Neurons, Astrocytes_in_Microglia))
  write.table(Tesselation_Row, Tesselation_CSV_Path, append = TRUE, sep=",", col.names = FALSE)

  Astrocytes_Neurons_Correlation <- ppm(Astrocytes_PPP ~ Neurons_Density)
  Microglia_Neurons_Correlation <- ppm(Microglia_PPP ~ Neurons_Density)
  Astrocytes_Microglia_Correlation <- ppm(Astrocytes_PPP ~ Microglia_Density)

  Neurons_Astrocytes_Function <- effectfun(Astrocytes_Neurons_Correlation)
  Microglia_Neurons_Function <- effectfun(Microglia_Neurons_Correlation)
  Astrocytes_Microglia_Function <- effectfun(Astrocytes_Microglia_Correlation)

  Neurons_Astrocytes_Function_Vector <<- append(Neurons_Astrocytes_Function_Vector, Neurons_Astrocytes_Function)
  Microglia_Neurons_Function_Vector <<- append(Microglia_Neurons_Function_Vector, Microglia_Neurons_Function)
  Astrocytes_Microglia_Function_Vector <<- append(Astrocytes_Microglia_Function, Astrocytes_Microglia_Function)


  #NNdist
  marks(Neurons_PPP) <- nndist(Neurons_PPP, K=5)
  Neurons_nndist <- as.data.frame(Neurons_PPP)

  marks(Astrocytes_PPP) <- nndist(Astrocytes_PPP, K=5)
  Astrocytes_nndist <- as.data.frame(Astrocytes_PPP)

  marks(Microglia_PPP) <- nndist(Microglia_PPP, K=5)
  Microglia_nndist <- as.data.frame(Microglia_PPP)

  Neurons_nndist$Brain <- path
  write.table(Neurons_nndist, Neurons_nndist_CSV_Path, append = TRUE, sep=",", col.names = FALSE)

  Astrocytes_nndist$Brain <- path
  write.table(Astrocytes_nndist, Astrocytes_nndist_CSV_Path, append = TRUE, sep=",", col.names = FALSE)

  Microglia_nndist$Brain <- path
  write.table(Microglia_nndist, Microglia_nndist_CSV_Path, append = TRUE, sep=",", col.names = FALSE)

  fragments <- strsplit(path, "_")[[1]]
  len <- length(fragments)
  mouse <- fragments[1]
  dpi <- fragments[2]
  hemis <- fragments[3]

  add_to_hyperframe(Neurons = Neurons_PPP, Astrocytes = Astrocytes_PPP, Microglia = Microglia_PPP, Neurons_Dens = Neurons_Density, Astrocytes_Dens = Astrocytes_Density, Microglia_Dens = Microglia_Density, Neurons_Tess = Neurons_Tess, Microglia_Tess = Microglia_Tess, ID = mouse, DPI=dpi, Hemis = hemis, stringsAsFactors=TRUE)
}

csv_files <- list.files(coordinatesPath, full.names = FALSE, recursive = FALSE)

brains <- c()

for (csv in csv_files) {
  fragments <- strsplit(csv, "_")[[1]]
  len <- length(fragments)
  brain_name <- paste(fragments[1:(len-2)], collapse="_")
  brains <- append(brains, brain_name)
}

brains <- unique(brains)

for (brain in brains) {
  process_file(coordinatesPath, brain)
}

saveRDS(Result_Hyperframe, "Hyperframes/PointPatterns_5x.rds")
