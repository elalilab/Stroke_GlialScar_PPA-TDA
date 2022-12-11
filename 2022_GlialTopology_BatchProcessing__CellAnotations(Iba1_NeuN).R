library("readr")
library("dplyr")

append_annotations <- function(base_path, brain_name, results_path) {
  
  Iba1_csv_path <- paste0(results_path, "/Iba1_Summary.csv")
  Annotations_Path <- paste0(base_path, "/", brain_name, "/Iba1/")
  process_annotation(results_path = Iba1_csv_path, Annotations_Path)
  
  NeuN_csv_path <- paste0(results_path, "/NeuN_Summary.csv")
  Annotations_Path <- paste0(base_path, "/", brain_name, "/NeuN/")
  process_annotation(results_path = NeuN_csv_path, Annotations_Path)
}

process_annotation <- function(results_path, path) {
  
  print (path)
  
  Annotations <- list.files(path = path, pattern = "annotations.tsv", full.names = TRUE) %>% 
    lapply(read_tsv) %>%                              
    bind_rows
 
   print(Annotations)
  
  Annotations <- as.data.frame(Annotations)
  names(Annotations) <- NULL
  
  write.table(Annotations, results_path, append = TRUE, sep=",")
}


basePath <- "D:/Daniel/Project-ECM/2.Images/Exp2-Gfap,NeuN,Iba1_5x/QupathProjects"
resultsPath <- "D:/Daniel/Project-ECM/3.DataAnalysis/Exp2-Gfap,NeuN,Iba1_5x/Results"

Iba1_csv_path <- paste0(resultsPath, "/Iba1_Summary.csv")
NeuN_csv_path <- paste0(resultsPath, "/NeuN_Summary.csv")

Annotations_Header <- c("Image",	"Name",	"Class",	"Parent",	"ROI",	"Centroid X ?m",	"Centroid Y ?m",	"ID",	"Parent ID",	"Side",	"Num Detections",	"Num Negatie", "Num Positive", "Positive %", "Num Positive per mm^2",  "Area ?m^2",	"Perimeter ?m")

df_header <- data.frame(matrix(ncol = 17, nrow = 0))
names(df_header) <- Annotations_Header

write.csv(df_header, Iba1_csv_path)
write.csv(df_header, NeuN_csv_path)

brains <- list.dirs(basePath, full.names = FALSE, recursive = FALSE)


for (brain in brains){
  append_annotations(basePath, brain, resultsPath)
}