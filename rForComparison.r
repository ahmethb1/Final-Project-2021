library(readr)
library(ape)
tree <- read.tree("../all_names.nwk")

names <- read_delim("../names.dmp", "\t", escape_double = FALSE, 
                    col_names = FALSE, trim_ws = TRUE, 
                    col_types = cols(X2=col_skip(), X4=col_skip(), X5=col_skip(),
                                     X6=col_skip(), X7=col_skip(), X8=col_skip()))
colnames(names) <- c("taxid", "txt")

names$txt <- gsub(" ", "_", names$txt)
names$taxid <- paste("t", names$taxid, sep="")

asmbl <- read.delim("../Organisms for evaluation.txt", header = T)
asmbl$taxid <- sub("t", "", asmbl$taxid)
asmbl$species_taxid <- sub("t", "", asmbl$species_taxid)
asmbl <- asmbl[!apply(asmbl=="", 1, all),]

#or as.character(levels(asmbl$species_taxid))[as.integer(asmbl$species_taxid)]
taxid <- unique(c(asmbl$species_taxid, asmbl$taxid))
taxid <- paste("t",unique(c(asmbl$taxid, asmbl$species_taxid)), sep = "")

asmbl$taxid <- paste("t", asmbl$taxid, sep= "")
asmbl$species_taxid <- paste("t",asmbl$species_taxid, sep = "")

names_in_ours <- names$taxid %in% asmbl$species_taxid
names_in_tree <- names$txt %in% tree$tip.label
keep <- subset(names, names_in_ours & names_in_tree)

t_name <- keep$txt
names(t_name) <- keep$taxid

asmbl <- subset(asmbl, taxid %in% keep$taxid | species_taxid %in% keep$taxid)
asmbl$sname <- t_name[asmbl$species_taxid]

tree_in_ours <- tree$tip.label %in% asmbl$sname
tt <- drop.tip(tree, tree$tip.label[!tree_in_ours])

dist <- cophenetic.phylo(tt)

pp <- rep(T, length(asmbl$sname))
ppn <- rep(NA, length(asmbl$sname))
ppi<- 1

for(i in 1:length(asmbl$sname)){
  if(!(T %in%(asmbl$sname[[i]] %in% ppn))){
    ppn[ppi] <- asmbl$sname[[i]]
    pp[i] <- F
    ppi <- ppi +1
  }
}

asmbl <- asmbl[!pp,]

#Changing matrix col and row names to GCF numbers

colnames(dist) <- asmbl$assembly_accession[match(colnames(dist), asmbl$sname)]
rownames(dist) <- asmbl$assembly_accession[match(rownames(dist), asmbl$sname)]

#Removing the organisms from GenDisCal output file which are not present in matrix
gen_dis <- read.csv("../GenDisCal_outputs/output_karl4.csv", header = T)
genetic_distance <- gen_dis[(gen_dis$File1 %in% colnames(dist)) & (gen_dis$File2 %in% colnames(dist)), ]

#Converting matrix to phylogenetic distance list
t = rep(NA,507*506/2)
k <- 1
for (i in 1:506){
  for (j in (i+1):507){
    x <- colnames(dist)[i]
    y <- rownames(dist)[j]
    z <- dist[i,j]
    t[k] <- paste(x,y,z, sep = ",")
    k <- k+1
  }
}

write(t, file = "../read_tree/tree_distances.txt")

#Sorting Gendiscal and Phylogenetic list files
##Takes phylogenetic distances, gives order names and distances which depend on one variable

tree_distance <- read.csv("../read_tree/tree_distances.txt", header = F)
len <- length(tree_distance$V1)
name_tree <- rep(NA, (len-1))

for (i in 1:len){
  p <- sort(c(as.character(tree_distance$V1[[i]]), as.character(tree_distance$V2[[i]])))
  name_tree[i] <- paste(p[1], p[2], sep = "_")
}
tree_distance$name_tree <- name_tree
tree_distance$V2 <-  tree_distance$V2[match(sort(tree_distance$name_tree), tree_distance$name_tree)]
tree_distance$V1 <-  tree_distance$V1[match(sort(tree_distance$name_tree), tree_distance$name_tree)]
tree_distance$V3 <-  tree_distance$V3[match(sort(tree_distance$name_tree), tree_distance$name_tree)]
tree_distance$name_tree <- sort(tree_distance$name_tree)

##Takes genetic distances, gives order names and distances which depend on one variable

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

write.csv(tree_distance, file = "../read_tree/sorted_tree.csv", row.names = F)
write.csv(genetic_distance, file = "../GenDisCal_outputs/sorted_output_karl4.csv", row.names = F)
##Are all sorted names equal to each other according to position

(length(tree_distance$name_tree) == sum(tree_distance$name_tree == genetic_distance$name_genetic)) & (length(genetic_distance$name_genetic) == sum(tree_distance$name_tree == genetic_distance$name_genetic))

#Plotting

plot((genetic_distance$Distance), log(tree_distance$V3), pch = 19, cex = 0.7 
     ,col = ifelse(genetic_distance$Expected_Relation == "Different_Phyla", 'red', 
                   ifelse(genetic_distance$Expected_Relation == "Same_Phylum", "yellow",
                          ifelse(genetic_distance$Expected_Relation == "Same_Class","gray",
                                 ifelse(genetic_distance$Expected_Relation == "Same_Order", "green",
                                        ifelse(genetic_distance$Expected_Relation == "Same_Family", "blue",
                                               ifelse(genetic_distance$Expected_Relation == "Same_Genus", "purple",
                                                      ifelse(genetic_distance$Expected_Relation == "Same_Species", "orange", "black")))))))
)
lines(lowess(log(tree_distance$V3) ~ (genetic_distance$Distance)), lwd = 2, col = "red")
legend("topleft", c("Different Phyla", "Same Phylum","Same Class","Same Order","Same Family","Same Genus", "Same Species"), fill = c("red", "yellow",
"gray", "green", "blue", "purple", "orange"))
#plot(2.718282^(7.9*genetic_distance$Distance + 1.889), (tree_distance$V3), pch = 19, cex = 0.7)
