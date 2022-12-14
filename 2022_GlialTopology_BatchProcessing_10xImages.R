set.seed(88071)
library(spatstat)
library(ggplot2)
library(dplyr)
library(purrr)

coordinatesPath <- "CellCoordinates_10x"

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

# Cretate a point pattern (PPP) object

create_point_pattern <- function(Subset, ReferenceSubset) {
  # We define the limits of the window according to Neuron coordinates
  xlim <- range(ReferenceSubset$X)
  ylim <- range(ReferenceSubset$Y)

  # Create point pattern for neurons
  Cells_PPP <- with(Subset, ppp(x = Subset$X, y = Subset$Y, xrange = xlim, yrange = ylim))
  unitname(Cells_PPP)  <- list("mm", "mm", 3.221/7096)
  Cells_PPP <- spatstat.geom::rescale (Cells_PPP)
  
  ## We rescale the unit to obtain measurements in mm2
  return(Cells_PPP)
  
}

tesselation <- function(Cells_Density, Cells_Quantiles) {
  ## We define the cutting spots according to quantiles
  Cells_Cut <- cut(Cells_Density, breaks = Cells_Quantiles, labels = 1:3)

  ## We generate the tesselation image
  return(tess(image = Cells_Cut))
}

process_file <- function (basePath, path) {

  Dapi_Raw <- read.csv(file = paste0(basePath, '/', path, '_10x_Dapi_Coordinates.csv'), header = TRUE)
  Neurons_Raw <- read.csv(file = paste0(basePath, '/', path, '_10x_NeuN_Coordinates.csv'), header = TRUE)
  Astrocytes_Raw <- read.csv(file = paste0(basePath, '/', path, '_10x_Gfap_Coordinates.csv'), header = TRUE)
  Microglia_Raw <- read.csv(file = paste0(basePath, '/', path, '_10x_Iba1_Coordinates.csv'), header = TRUE)

  # Subset neurons
  Neurons_Subset <- Neurons_Raw[(Neurons_Raw$Class=="Positive"),]
 
  # We subset astrocytes
  Astrocytes_Subset <- Astrocytes_Raw[(Astrocytes_Raw$Class=="Positive"),]
 
  # We subset microglia
  Microglia_Subset <- Microglia_Raw[(Microglia_Raw$Class=="Positive"),]
  
  # We subset Dapi
  Dapi_Subset <- Dapi_Raw[(Dapi_Raw$Class=="Positive"),]

  Neurons_PPP <- create_point_pattern(Neurons_Subset, Dapi_Subset)
  Astrocytes_PPP <- create_point_pattern(Astrocytes_Subset, Dapi_Subset)
  Microglia_PPP <- create_point_pattern(Microglia_Subset, Dapi_Subset)

  Microglia_Density <- density(Microglia_PPP, sigma =0.03, positive = TRUE, equal.ribbon = TRUE, col = topo.colors, main = "")
  Astrocytes_Density <- density(Astrocytes_PPP, sigma =0.03, positive = TRUE, equal.ribbon = TRUE, col = topo.colors, main = "")
  Neurons_Density <- density(Neurons_PPP, sigma =0.03, positive = TRUE, equal.ribbon = TRUE, col = topo.colors, main = "")

  Neurons_Tess <- tesselation(Neurons_Density, c(0, 5, 10, 20))
  Microglia_Tess <- tesselation(Microglia_Density, c(0, 10, 20, 35))

  
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
  brain_name <- paste(fragments[1:3], collapse="_")
  brains <- append(brains, brain_name)
}

brains <- unique(brains)

for (brain in brains) {
  process_file(coordinatesPath, brain)
}

saveRDS(Result_Hyperframe, "Hyperframes/PointPatterns_10x.rds")

