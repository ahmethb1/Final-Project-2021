#--------------------------------
install.packages("tidyverse")
install.packages("ggpubr")
#--------------------------------

library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)

distances <- read.csv('output.csv')


ggplot(distances, aes(x=Distance, color=Expected_Relation,
                      fill=Expected_Relation)) + 
  geom_histogram(position = "identity", alpha=0.5) + theme_bw()

#--------------------------------

ggplot(distances, aes(x=Distance, color=Expected_Relation, 
                      fill=Expected_Relation)) + 
  geom_histogram(position = "identity", alpha=0.5) + 
  scale_color_brewer(palette="Dark2")+
       scale_fill_brewer(palette="Dark2") + theme_bw()

#--------------------------------

karl4 <- read.csv('sorted_output_karl4.csv', stringsAsFactors = F)


q <- ggplot(karl4, aes(x=Distance, color=Expected_Relation,
                  fill=Expected_Relation)) + labs(y= "Frequency", x = "Genetic Distance") +
  geom_histogram(position = "identity", alpha=0.7, bins = 200) + theme_bw() +
  scale_color_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")

w <- karl4 %>% filter(Expected_Relation!="?") %>% 
  ggplot(aes(x=Distance, color=Expected_Relation,fill=Expected_Relation)) + labs(y= "Frequency", x = "Genetic Distance") +
  geom_histogram(position = "identity", alpha=0.7, bins = 200) + theme_bw() +
  scale_color_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")

dstrbtn <- plot_grid(q, w + rremove("x.text"), 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)

annotate_figure(dstrbtn,
                top = text_grob("Distribution of Distances - Karlin 4", color = "black", face = "bold", size = 14))

#--------------------------------
RColorBrewer::display.brewer.all()
#--------------------------------
treedist <- read.csv('sorted_tree.csv')

bacteria_treedist <- read.csv('sorted_bac_tree.csv')


freq4 <- read.csv('sorted_freq4_PaSiT.csv', stringsAsFactors = F)
karl4 <- read.csv('sorted_karl4_PaSiT.csv', stringsAsFactors = F)
exp4 <- read.csv('sorted_exp4_hamming.csv', stringsAsFactors = F)
GC <- read.csv('sorted_GC_manhattan.csv', stringsAsFactors = F)
len <- read.csv('sorted_len_pearson.csv', stringsAsFactors = F)
bacteria_karl4 <- read.csv('sorted_bac_karl4.sdb.csv', stringsAsFactors = F)


dffreq4 <- data.frame(freq4, Divergence.time=treedist$V3)
dfkarl4 <- data.frame(karl4, Divergence.time=treedist$V3)
dfexp4 <- data.frame(exp4, Divergence.time=treedist$V3)
dfGC <- data.frame(GC, Divergence.time=treedist$V3)
dflen <- data.frame(len, Divergence.time=treedist$V3)
dfbactkarl4 <- data.frame(bacteria_karl4, Divergence.time=bacteria_treedist$V3)
#--------------------------------
correct_order <- c("Different_Phyla", "Same_Phylum", "Same_Class","Same_Order","Same_Family", "Same_Genus")
dffreq4$Expected_Relation <- factor(dffreq4$Expected_Relation, levels = correct_order)
dfkarl4$Expected_Relation <- factor(dfkarl4$Expected_Relation, levels = correct_order)
dfexp4$Expected_Relation <- factor(dfexp4$Expected_Relation, levels = correct_order)
dfGC$Expected_Relation <- factor(dfGC$Expected_Relation, levels = correct_order)
dflen$Expected_Relation <- factor(dflen$Expected_Relation, levels = correct_order)
dfbactkarl4$Expected_Relation <- factor(dfbactkarl4$Expected_Relation, levels = correct_order)
#--------------------------------


##rlevel


 
a <- dffreq4 %>% filter(Expected_Relation!="?") %>% 
  ggplot(aes(x=Distance,y=Divergence.time, 
                   color=Expected_Relation)) + labs(y= "Divergence Time", x = "Genetic Distance") +
  geom_point() + scale_y_log10() +ggtitle("Frequency 4") + scale_colour_discrete("Expected Relations") +
  theme(plot.title = element_text(hjust = 0.5))

b <- dfkarl4 %>% filter(Expected_Relation!="?") %>% 
  ggplot(aes(x=Distance,y=Divergence.time, 
             color=Expected_Relation)) + labs(y= "Divergence Time", x = "Genetic Distance") +
  geom_point() + scale_y_log10() +ggtitle("Karlin 4") + scale_colour_discrete("Expected Relations") +
  theme(plot.title = element_text(hjust = 0.5)) 

c <- dfexp4 %>% filter(Expected_Relation!="?") %>% 
  ggplot(aes(x=Distance,y=Divergence.time, 
             color=Expected_Relation)) + labs(y= "Divergence Time", x = "Genetic Distance") +
  geom_point() + scale_y_log10() +ggtitle("Markov 4") + scale_colour_discrete("Expected Relations") +
  theme(plot.title = element_text(hjust = 0.5)) 

d <- dfGC %>% filter(Expected_Relation!="?") %>% 
  ggplot(aes(x=Distance,y=Divergence.time, 
             color=Expected_Relation)) + labs(y= "Divergence Time", x = "Genetic Distance") +
  geom_point() + scale_y_log10() +ggtitle("GC") + scale_colour_discrete("Expected Relations") +
  theme(plot.title = element_text(hjust = 0.5))

ggarrange(a, b, c, d + rremove("x.text"), 
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")

#--------------------------------

y <- dfbactkarl4 %>% filter(Divergence.time>100, Distance>0.5) %>% 
  ggplot(aes(x=Distance, y=Divergence.time, color=Expected_Relation)) + labs(y= "Divergence Time", x = "Genetic Distance") +
  geom_point() + scale_y_log10() + scale_colour_discrete("Expected Relations") +
  theme(plot.title = element_text(hjust = 0.5))



dfbactkarl4$Expected_Relation %>% table() -> w
dfbactkarl4 %>% filter(Divergence.time>100, Distance>0.5) %>% lm(log10(Divergence.time)~Distance, data=., weights = 1/w[as.character(Expected_Relation)]) -> model
summary(model)$r.squared

z <- dfbactkarl4 %>% filter(Divergence.time>100, Distance>0.5) %>% 
  mutate(Divergence=jitter(Divergence.time,factor = 1000)) %>% 
  ggplot(aes(x=Distance, y=Divergence, color=Expected_Relation)) + labs(y= "Divergence Time", x = "Genetic Distance") +
  geom_point(alpha=0.3) + ggtitle("",subtitle = paste("Adj R2 =", signif(summary(model)$r.squared,5)))  + scale_y_log10() +
  geom_line(mapping = aes(x=Distance ,y=10^predict(model)), inherit.aes=FALSE) + scale_colour_discrete("Expected Relations") + 
  
  theme(plot.title = element_text(hjust = 0.5))

figure <- plot_grid(y, z + rremove("x.text"), 
          labels = c("A", "B"),
          ncol = 1, nrow = 2, common.legend = TRUE, legend = "bottom")
annotate_figure(figure,
                top = text_grob("Karlin 4- Bacteria Only", color = "black", face = "bold", size = 14))


#--------------------------------
