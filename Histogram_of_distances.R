#--------------------------------
distances <- read.csv('output.csv')
par(mar = c(3,3,3,12), xpd = T)
hist(distances$Distance, breaks = length(distances$Distance), freq = F, xlim = c(min(distances$Distance), max(distances$Distance)),
     col = ifelse(distances$Expected_Relation == "Different_Phyla", 'red', 
            ifelse(distances$Expected_Relation == "Same_Phylum", "yellow",
            ifelse(distances$Expected_Relation == "Same_Class","gray",
            ifelse(distances$Expected_Relation == "Same_Order", "green",
            ifelse(distances$Expected_Relation == "Same_Family", "blue",
            ifelse(distances$Expected_Relation == "Same_Genus", "purple", "black")))))))

legend("topright", c("Different_Phyla", "Same_Phylum","Same_Class","Same_Order","Same_Family","Same_Genus"), fill = c("red", "yellow",
                                                                                                     "gray", "green", "blue", "purple"),
                                                                                                      xpd = T, bty = 'n', inset = c(-0.35,0))

