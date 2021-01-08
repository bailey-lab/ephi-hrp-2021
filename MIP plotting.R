#E Reichert
#Nov 2020

##Background##

#This script visualizes MIP sequencing results, by sample, on a dichotomous scale (i.e., did the MIP bind its target or not?)
#using a rectangle plot
#These visualizations are used for the MIP-WGS Comparison figure in our manuscript Feleke et al, which characterizes
#different deletion breakpoint patterns alongs chr8/13 for 2017-18 samples collected in Ethiopia.

library(readr)
library(dplyr)
library(ggplot2)
library(readxl)
library(tidyr)
library(stringr)

## Chromosome 8 (pfhrp2) ##

#Read in Chr8 (pfhrp2) MIP data
chr8 <- read_excel("mip_0630_2020.xlsx", sheet = "chr8_combined")

#read in long format of data with log10 coverage produced in "EPHI MIPS Analysis.Rmd"
chr8_vis <- read_csv("MIPS_chr8_visual")

#clean up barcodes
#'Gene' column contains barcode + region information
chr8 <- chr8 %>% 
  mutate(Gene = str_replace_all(Gene, "-","")) %>%
  mutate(Gene = str_replace(Gene, "^E18","")) %>%
  mutate(Gene = str_replace(Gene, "HRP1$","")) %>%
  mutate(Gene = str_replace(Gene, "HRP2$",""))

chr8$BarcodeNumber <- stringr::str_extract(chr8$Gene, "^.{4}")

#separate out barcode number + region identifiers
chr8 <- chr8 %>%
  separate(Gene, c("BarcodeNumber", "Reg"), sep = 4,
           convert = FALSE, extra = "warn", fill = "warn")

#define region labels
chr8$Reg[chr8$Reg=="DT"] <- "Tg"
chr8$Reg[chr8$Reg=="T"] <- "Tg"
chr8$Reg[chr8$Reg=="DA"] <- "Am"
chr8$Reg[chr8$Reg=="A"] <- "Am"
chr8$Reg[chr8$Reg=="DG"] <- "GA"
chr8$Reg[chr8$Reg=="G"] <- "GA"

#extract probes genomic location from this dataset
#pull out rows 2-4 which have probe start location, end location, and mean location
chr8_probes <- chr8[c(2:4),c(1,3:ncol(chr8))]
chr8_probes <- t(chr8_probes)
rownames(chr8_probes) <- c()
colnames(chr8_probes) <- c("Start", "End", "Gene")
chr8_probes <- chr8_probes[c(2:nrow(chr8_probes)),]
chr8_probes <- as.data.frame(chr8_probes)
chr8_probes <- data.frame(lapply(chr8_probes, as.character), stringsAsFactors=FALSE)
chr8_probes <- data.frame(lapply(chr8_probes, as.numeric), stringsAsFactors=FALSE)
chr8_probes$Gene <- round(chr8_probes$Gene, 0)

#set which sample to focus on
#in this case we are analysing barcode 2463, replace all to look at new barcode
mip2463 <- subset(chr8.long2, BarcodeNumber == 2463)
mip2463 <- mip2463 %>% select(BarcodeNumber, Gene, log10_coverage)

mip2463$Gene <- as.character(mip2463$Gene)

#change MIP target names covering hrp2 to their mean location
mip2463$Gene[mip2463$Gene == "HRP2_1"] <- 1374360
mip2463$Gene[mip2463$Gene == "HRP2_2"] <- 1374409
mip2463$Gene[mip2463$Gene == "HRP2_3"] <- 1374528
mip2463$Gene[mip2463$Gene == "HRP2_4"] <- 1374564
mip2463$Gene[mip2463$Gene == "HRP2_5"] <- 1374712
mip2463$Gene[mip2463$Gene == "HRP2_6"] <- 1374726
mip2463$Gene[mip2463$Gene == "HRP2_7"] <- 1374874
mip2463$Gene[mip2463$Gene == "HRP2_8"] <- 1374959

mip2463$Gene <- as.numeric(mip2463$Gene)

#merge with probe data by 'gene' (probe) to get start and end location
mip2463 <- left_join(mip2463, chr8_probes, by = "Gene")

#set coverage for each probe to be neg
mip2463$positive <- c(rep('No', nrow(mip2463)))

#replace with 'yes' if coverage > 0.1 (any coverage)
mip2463$positive[mip2463$log10_coverage > -1] <- 'Yes'
mip2463$ylim <- c(rep(0, nrow(mip2463)))
mip2463$ymax <- c(rep(0.1, nrow(mip2463)))

#Generate figure
#blue = MIP did not bind target, orange = MIP bound target
#this visualizes all of the probes along chr8 with true spacing
#pdf(file = "mip2463.pdf", width = 11, height = 1)
ggplot() + 
  geom_rect(data=mip2463, mapping=aes(xmin=Start, xmax=End, ymin= ylim, ymax=ymax, fill=positive)) + ylab('') + xlab('Chromosome 8 Position') +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 5, angle = 90),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = 'none') + 
  scale_fill_manual(values = c("blue4", "orange")) + scale_x_continuous(limits = c(1326000, 1420600), breaks = c(seq(1345000, 1402000, 2000))) + labs(fill = "MIP Positive")
#dev.off()

## Chromosome 13 (pfhrp3) ##

#Read in Chr13 (HRP3) data
chr13 <- read_excel("mip_0630_2020.xlsx", sheet = "chr13_combined")

#read in long format of data with log10 coverage produced in "EPHI MIPS Analysis.Rmd"
chr13_vis <- read_csv("MIPS_chr13_visual")

#clean up barcodes
chr13 <- chr13 %>% 
  mutate(Gene = str_replace_all(Gene, "-","")) %>%
  mutate(Gene = str_replace(Gene, "^E18","")) %>%
  mutate(Gene = str_replace(Gene, "HRP1$","")) %>%
  mutate(Gene = str_replace(Gene, "HRP2$",""))

chr13$BarcodeNumber <- stringr::str_extract(chr13$Gene, "^.{4}")

#separate out barcode number + region identifiers
chr13 <- chr13 %>%
  separate(Gene, c("BarcodeNumber", "Reg"), sep = 4,
           convert = FALSE, extra = "warn", fill = "warn")

#define region labels
chr13$Reg[chr13$Reg=="DT"] <- "Tg"
chr13$Reg[chr13$Reg=="T"] <- "Tg"
chr13$Reg[chr13$Reg=="DA"] <- "Am"
chr13$Reg[chr13$Reg=="A"] <- "Am"
chr13$Reg[chr13$Reg=="DG"] <- "GA"
chr13$Reg[chr13$Reg=="G"] <- "GA"

#extract probes genomic location from this dataset
#pull out rows 2-4 which have probe start location, end location, and mean location
chr13_probes <- chr13[c(2:4),c(1,3:ncol(chr13))]
chr13_probes <- t(chr13_probes)
rownames(chr13_probes) <- c()
colnames(chr13_probes) <- c("Start", "End", "Gene")
chr13_probes <- chr13_probes[c(2:nrow(chr13_probes)),]
chr13_probes <- as.data.frame(chr13_probes)

chr13_probes <- data.frame(lapply(chr13_probes, as.character), stringsAsFactors=FALSE)
chr13_probes <- data.frame(lapply(chr13_probes, as.numeric), stringsAsFactors=FALSE)
chr13_probes$Gene <- round(chr13_probes$Gene, 0)

#select which sample to focus on
mip2463 <- subset(chr13.long2, BarcodeNumber == "2463")

mip2463 <- mip2463 %>% select(BarcodeNumber, Gene, log10_coverage)
mip2463$Gene <- as.character(mip2463$Gene)

#rename MIPs covering hrp3 to their mean location
mip2463$Gene[mip2463$Gene == "HRP3_1"] <- 2840751
mip2463$Gene[mip2463$Gene == "HRP3_2"] <- 2840852
mip2463$Gene[mip2463$Gene == "HRP3_3"] <- 2840904
mip2463$Gene[mip2463$Gene == "HRP3_4"] <- 2841107
mip2463$Gene[mip2463$Gene == "HRP3_5"] <- 2841240
mip2463$Gene[mip2463$Gene == "HRP3_6"] <- 2841341
mip2463$Gene[mip2463$Gene == "HRP3_7"] <- 2841528
mip2463$Gene[mip2463$Gene == "HRP3_8"] <- 2841558
mip2463$Gene[mip2463$Gene == "HRP3_9"] <- 2841762

mip2463$Gene <- as.numeric(mip2463$Gene)

#merge in MIP probes gene location info
mip2463 <- left_join(mip2463, chr13_probes, by = "Gene")

mip2463$positive <- c(rep(0, nrow(mip2463)))
mip2463$positive[mip2463$log10_coverage > -1] <- 1
mip2463$ylim <- c(rep(0, nrow(mip2463)))
mip2463$ymax <- c(rep(0.1, nrow(mip2463)))

#produce plot
#pdf(file = "mip2463_hrp3.pdf", width = 10, height = 0.8)
ggplot() + 
  geom_rect(data=mip2463, mapping=aes(xmin=Start, xmax=End, ymin= ylim, ymax=ymax, fill=factor(positive))) + ylab('') + 
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 5, angle = 90),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") + 
  scale_fill_manual(values = c("blue4", "orange")) + scale_x_continuous(limits = c(2776000, 2890000), breaks = c(seq(2774000, 2900000, 2000))) + labs(fill = "MIP Positive")
#dev.off()
