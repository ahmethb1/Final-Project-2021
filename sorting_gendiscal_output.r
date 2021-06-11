library(readr)
#Removing the organisms from GenDisCal output file which are not present in matrix

load("../timetree-dist.Rdata")
list1 <- read.delim("C:/Users/Ramal/Documents/GenDisCal_outputs/outputs/list.txt", sep = "\n", stringsAsFactors = F)
list1 <- sort(list1$file)
##copy the paths of files need to be sorted and trimed
paths <- readClipboard()
paths <- gsub('\"', '', paths)

for (j in 1:length(paths)){
  genetic_distance <- read.csv(paths[j], header = T)
  genetic_distance <- genetic_distance[(genetic_distance$File1 %in% colnames(dist)) & (genetic_distance$File2 %in% colnames(dist)), ]
  
  len <- length(genetic_distance$File1)
  name_genetic <- rep(NA, len)
  
  for(i in 1:len){
    p <- sort(c(as.character(genetic_distance$File1[[i]]), as.character(genetic_distance$File2[[i]])))
    name_genetic[i] <- paste(p[1], p[2], sep = "_")
  }
  genetic_distance$name_genetic <- name_genetic
  genetic_distance$File2 <-  genetic_distance$File2[match(sort(genetic_distance$name_genetic), genetic_distance$name_genetic)]
  genetic_distance$File1 <-  genetic_distance$File1[match(sort(genetic_distance$name_genetic), genetic_distance$name_genetic)]
  genetic_distance$Distance <-  genetic_distance$Distance[match(sort(genetic_distance$name_genetic), genetic_distance$name_genetic)]
  genetic_distance$Expected_Relation <-  genetic_distance$Expected_Relation[match(sort(genetic_distance$name_genetic), genetic_distance$name_genetic)]
  genetic_distance$name_genetic <- sort(genetic_distance$name_genetic)
  
  write.csv(genetic_distance, file = paste("C:\\Users\\Ramal\\Documents\\GenDisCal_outputs\\outputs\\sorted_", list1[j], ".csv", sep = ""), row.names = F, quote = F)
}
