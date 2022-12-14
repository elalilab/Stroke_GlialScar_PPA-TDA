library("readr")
library("dplyr")

append_annotations <- function(base_path, brain_name, results_path) {
  
  Gfap_csv_path <- paste0(results_path, "/Gfap_Summary.csv")
  Annotations_Path <- paste0(base_path, "/", brain_name, "/Gfap/")
  process_annotation(results_path = Gfap_csv_path, Annotations_Path)
  
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

basePath <- "QupathProjects_5x"
resultsPath <- "ResultsTables"

Annotations_Header <- c("Image",	"Name",	"Class",	"Parent",	"ROI",	"Centroid X ?m",	"Centroid Y ?m",	"ID",	"Parent ID",	"Side",	"Num Detections",	"Area ?m^2",	"Perimeter ?m")

df_header <- data.frame(matrix(ncol = 13, nrow = 0))
names(df_header) <- Annotations_Header

write.csv(df_header, Gfap_csv_path)

brains <- list.dirs(basePath, full.names = FALSE, recursive = FALSE)

for (brain in brains){
  append_annotations(basePath, brain, resultsPath)
}