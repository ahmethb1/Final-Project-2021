genetic_distance <- read.csv("../GenDisCal_outputs/sorted_karl4.csv", header = T)
tree_distance<- read.csv("../read_tree/sorted_tree.csv", header = T)

genetic_distance$range[genetic_distance$Distance<0.1] <- "0-0.1"
genetic_distance$range[genetic_distance$Distance>0.1 & genetic_distance$Distance<0.2] <- "0.1-0.2"
genetic_distance$range[genetic_distance$Distance>0.2 & genetic_distance$Distance<0.3] <- "0.2-0.3"
genetic_distance$range[genetic_distance$Distance>0.3 & genetic_distance$Distance<0.4] <- "0.3-0.4"
genetic_distance$range[genetic_distance$Distance>0.4 & genetic_distance$Distance<0.5] <- "0.4-0.5"
genetic_distance$range[genetic_distance$Distance>0.5 & genetic_distance$Distance<0.6] <- "0.5-0.6"
genetic_distance$range[genetic_distance$Distance>0.6 & genetic_distance$Distance<0.7] <- "0.6-0.7"
genetic_distance$range[genetic_distance$Distance>0.7 & genetic_distance$Distance<0.8] <- "0.7-0.8"
genetic_distance$range[genetic_distance$Distance>0.8 & genetic_distance$Distance<0.9] <- "0.8-0.9"
genetic_distance$range[genetic_distance$Distance>0.9 & genetic_distance$Distance<1] <- "0.9-1"

boxplot(tree_distance$V3 ~ genetic_distance$range, log = "y",  yaxt= "n", ylim = c(10^(-1), 10^4),
        pch = 19, cex = 0.7,ylab = "Phylogenetic Distance", xlab = "Genetic Distance",
        main = " Karlin 4", cex.lab = 1.5 , cex.main =1.5)
axis(2, at = c(0.1, 1,10,100,1000,10000), labels = c(0.1, 1,10,100,1000,10000), las = 1)

