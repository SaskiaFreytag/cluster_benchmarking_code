## Overview

library(ggplot2)
library(reshape2)
library(readxl)
library(gridExtra)
library(cowplot)

overview <- read_xlsx("Overview.xlsx")
overview <- overview[!apply(overview, 1, function(x) all(is.na(x))),]

overview$`Stability Cells 1` <- factor(overview$`Stability Cells 1`, 
                                          levels=c("poor", "fair", "good"))
overview$`Stability Cells 2` <- factor(overview$`Stability Cells 2`, 
                                       levels=c("poor", "fair", "good"))
overview$`Stability Genes` <- factor(overview$`Stability Genes`, 
                                       levels=c("poor", "fair", "good"))
overview$`Stability Aligners` <- factor(overview$`Stability Aligners`, 
                                       levels=c("poor", "fair", "good"))


overview11 <- melt(overview[, c(1,2,3,4,9,10)])
colnames(overview11) <- c("Method", "Test", "Rank")

overview11$Method <- factor(overview11$Method, 
                            levels=sort(unique(overview11$Method), decreasing=TRUE))

g1 <- ggplot(overview11, aes(y=Method, x=Test)) + geom_tile(aes(fill = Rank),
colour = "white") + scale_fill_viridis_c(trans = 'reverse') +
  theme_minimal() + theme (axis.title.x = element_blank(), axis.line.y = element_blank(),
                axis.text.x = element_text(angle = 45, hjust = 1)) + guides(fill=FALSE)

g1a <- ggplot(overview11, aes(y=Method, x=Test)) + geom_tile(aes(fill = Rank),
                                                            colour = "white") + scale_fill_viridis_c(trans = 'reverse') +
  theme_minimal() + theme (axis.title.x = element_blank(), axis.line.y = element_blank(),
                           axis.text.x = element_text(angle = 45, hjust = 1))                                                                                                

overview12 <- melt(overview[, c(1,5,6,7,8)], id.vars="X__1")
colnames(overview12) <- c("Method", "Test", "Performance")

overview12$Method <- factor(overview12$Method, 
                            levels=sort(unique(overview12$Method), decreasing=TRUE))
overview12$Performance <- factor(overview12$Performance, levels=c("good", "fair", "poor"))


g2 <- ggplot(overview12, aes(y=Method, x=Test)) + geom_tile(aes(fill = Performance),
colour = "white") + scale_fill_viridis_d(begin=1, end=0, na.value = "grey50",
                                         option="C") +
  theme_minimal() + theme (axis.title.x = element_blank(), axis.line.y = element_blank(),
                           axis.text.x = element_text(angle = 45, hjust = 1),
                           axis.text.y = element_blank(), axis.title.y = element_blank()) +
  guides(fill=FALSE)

g2a <- ggplot(overview12, aes(y=Method, x=Test)) + geom_tile(aes(fill = Performance),
    colour = "white") + scale_fill_viridis_d(begin=1, end=0, na.value = "grey50", 
                                             option="C") +
  theme_minimal() + theme (axis.title.x = element_blank(), axis.line.y = element_blank(),
                           axis.text.x = element_text(angle = 45, hjust = 1),
                           axis.text.y = element_blank(), axis.title.y = element_blank()) 


overview_all <- data.frame(Method = overview$X__1,
  `Total Score` = rank (apply(overview[,c(2,3,4,9,10)],1 , function(x) mean(x))),
  `Total Stability` = order(apply(overview[, c(5,6,7,8)], 1, 
              function(x) sum(x=="good") -sum(x=="poor")), decreasing=T),
  check.names = FALSE)

overview_all$`Total Score` <- as.numeric(overview_all$`Total Score`)
overview_all$`Total Stability` <- as.numeric(overview_all$`Total Stability`)

overview_all$Rank <- rank(apply(overview_all, 1, function(x) sum(as.numeric(x[2:3]))))

overview_all$Method <- factor(overview_all$Method, 
                            levels=sort(unique(overview_all$Method), decreasing=TRUE))
overview_all$x <- "       Total Score " 

g3 <- ggplot(overview_all, aes(y=Method, x=x)) + 
  geom_tile(aes(fill = Rank), colour = "white") + 
  scale_fill_viridis_c(trans = 'reverse') +
  theme_minimal() + theme (axis.title.x = element_blank(), axis.line.y = element_blank(),
                           axis.text.x = element_text(angle = 45, hjust = 1),
                           axis.text.y = element_blank(), axis.title.y = element_blank()) + 
  guides(fill=FALSE)

g3a <- ggplot(overview_all, aes(y=Method, x=x)) + 
  geom_tile(aes(fill = Rank), colour = "white") + 
  scale_fill_viridis_c(trans = 'reverse') +
  theme_minimal() + theme (axis.title.x = element_blank(), axis.line.y = element_blank(),
                           axis.text.x = element_text(angle = 45, hjust = 1),
                           axis.text.y = element_blank(), axis.title.y = element_blank())

gg2 <- cowplot::get_legend(g2a)
gg1 <- cowplot::get_legend(g1a)

gg <- grid.arrange(g1, g2, g3, gg2, gg1, ncol=4, widths=c(5,3,1,2),
             layout_matrix = rbind(c(1, 2, 3, 4),
                                   c(1, 2, 3, 5),
                                   c(1, 2, 3, NA)))
ggsave("Overview.pdf", gg, width=6, height=5)
