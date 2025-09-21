library(ggplot2)
library(superb)
library(dplyr)
setwd("~/Downloads")
Bdepth <- read.delim("L_dussii_B.genomemapped.depth.100b.bed", header = F, sep = "\t")
Cdepth <- read.delim("L_dussii_C.genomemapped.depth.100b.bed", header = F, sep = "\t")
Ddepth <- read.delim("L_dussii_D.genomemapped.depth.100b.bed", header = F, sep = "\t")
Gdepth <- read.delim("L_dussii_G.genomemapped.depth.100b.bed", header = F, sep = "\t")
JC1depth <- read.delim("L_dussii_JC1.genomemapped.depth.100b.bed", header = F, sep = "\t")
JC2depth <- read.delim("L_dussii_JC2.genomemapped.depth.100b.bed", header = F, sep = "\t")

colnames(Bdepth) <- c("SampleChrom","StartBase","StopBase" ,"SampleDepth")
colnames(Cdepth) <- c("SampleChrom","StartBase","StopBase" , "SampleDepth")
colnames(Ddepth) <- c("SampleChrom","StartBase", "StopBase" ,"SampleDepth")
colnames(Gdepth) <- c("SampleChrom","StartBase", "StopBase" ,"SampleDepth")
colnames(JC1depth) <- c("SampleChrom","StartBase", "StopBase" ,"SampleDepth")
colnames(JC2depth) <- c("SampleChrom","StartBase", "StopBase" ,"SampleDepth")

Bdepth$SampleDepth <- as.numeric(Bdepth$SampleDepth)
Cdepth$SampleDepth <- as.numeric(Cdepth$SampleDepth)
Ddepth$SampleDepth <- as.numeric(Ddepth$SampleDepth)
Gdepth$SampleDepth <- as.numeric(Gdepth$SampleDepth)
JC1depth$SampleDepth <- as.numeric(JC1depth$SampleDepth)
JC2depth$SampleDepth <- as.numeric(JC2depth$SampleDepth)

HsubsetBdepth <- read.delim("L_dussii_H.Bsubset.genomemapped.depth.100b.bed", header = F, sep = "\t")
HsubsetCdepth <- read.delim("L_dussii_H.Csubset.genomemapped.depth.100b.bed", header = F, sep = "\t")
HsubsetDdepth <- read.delim("L_dussii_H.Dsubset.genomemapped.depth.100b.bed", header = F, sep = "\t")
HsubsetGdepth <- read.delim("L_dussii_H.Gsubset.genomemapped.depth.100b.bed", header = F, sep = "\t")
HsubsetJC1depth <- read.delim("L_dussii_H.JC1subset.genomemapped.depth.100b.bed", header = F, sep = "\t")
HsubsetJC2depth <- read.delim("L_dussii_H.JC2subset.genomemapped.depth.100b.bed", header = F, sep = "\t")

colnames(HsubsetBdepth) <- c("RefChrom","RefStart", "RefStop", "RefDepth")
colnames(HsubsetCdepth) <- c("RefChrom","RefStart", "RefStop","RefDepth")
colnames(HsubsetDdepth) <- c("RefChrom","RefStart", "RefStop","RefDepth")
colnames(HsubsetGdepth) <- c("RefChrom","RefStart", "RefStop","RefDepth")
colnames(HsubsetJC1depth) <- c("RefChrom","RefStart", "RefStop","RefDepth")
colnames(HsubsetJC2depth) <- c("RefChrom","RefStart", "RefStop","RefDepth")

HsubsetBdepth$RefDepth <- as.numeric(HsubsetBdepth$RefDepth)
HsubsetCdepth$RefDepth <- as.numeric(HsubsetCdepth$RefDepth)
HsubsetDdepth$RefDepth <- as.numeric(HsubsetDdepth$RefDepth)
HsubsetGdepth$RefDepth <- as.numeric(HsubsetGdepth$RefDepth)
HsubsetJC1depth$RefDepth <- as.numeric(HsubsetJC1depth$RefDepth)
HsubsetJC2depth$RefDepth <- as.numeric(HsubsetJC2depth$RefDepth)

Bmerged <- cbind(Bdepth, HsubsetBdepth)
#BmergedSample <- sample_n(Bmerged, 10000000)
#rm(Bdepth, HsubsetBdepth, Bmerged)

Cmerged <- cbind(Cdepth, HsubsetCdepth)
#CmergedSample <- sample_n(Cmerged, 10000000)
#rm(Cdepth, HsubsetCdepth, Cmerged)

Dmerged <- cbind(Ddepth, HsubsetDdepth)
#DmergedSample <- sample_n(Dmerged, 10000000)
#rm(Ddepth, HsubsetDdepth, Dmerged)

Gmerged <- cbind(Gdepth, HsubsetGdepth)
#GmergedSample <- sample_n(Gmerged, 10000000)
#rm(Gdepth, HsubsetGdepth, Gmerged)

JC1merged <- cbind(JC1depth, HsubsetJC1depth)
#JC1mergedSample <- sample_n(JC1merged, 10000000)
#rm(JC1depth, HsubsetJC1depth, JC1merged)

JC2merged <- cbind(JC2depth, HsubsetJC2depth)
#JC2mergedSample <- sample_n(JC2merged, 10000000)
#rm(JC2depth, HsubsetJC2depth, JC2merged)

Bmerged_rmmissing <- Bmerged[Bmerged$RefDepth != "0", ]
Cmerged_rmmissing <- Cmerged[Cmerged$RefDepth != "0", ]
Dmerged_rmmissing <- Dmerged[Dmerged$RefDepth != "0", ]
Gmerged_rmmissing <- Gmerged[Gmerged$RefDepth != "0", ]
JC1merged_rmmissing <- JC1merged[JC1merged$RefDepth != "0", ]
JC2merged_rmmissing <- JC2merged[JC2merged$RefDepth != "0", ]

Bmerged_rmNA <- Bmerged_rmmissing[Bmerged_rmmissing$RefDepth != "NaN", ]
Cmerged_rmNA <- Cmerged_rmmissing[Cmerged_rmmissing$RefDepth != "NaN", ]
Dmerged_rmNA <- Dmerged_rmmissing[Dmerged_rmmissing$RefDepth != "NaN", ]
Gmerged_rmNA <- Gmerged_rmmissing[Gmerged_rmmissing$RefDepth != "NaN", ]
JC1merged_rmNA <- JC1merged_rmmissing[JC1merged_rmmissing$RefDepth != "NaN", ]
JC2merged_rmNA <- JC2merged_rmmissing[JC2merged_rmmissing$RefDepth != "NaN", ]

Bmerged_rmOutliers <- Bmerged_rmmissing[which((Bmerged_rmmissing$SampleDepth < 100 ) & (Bmerged_rmmissing$RefDepth < 100)), ]
Cmerged_rmOutliers <- Cmerged_rmmissing[which((Cmerged_rmmissing$SampleDepth < 100 ) & (Cmerged_rmmissing$RefDepth < 100)), ]
Dmerged_rmOutliers <- Dmerged_rmmissing[which((Dmerged_rmmissing$SampleDepth < 100 ) & (Dmerged_rmmissing$RefDepth < 100)), ]
Gmerged_rmOutliers <- Gmerged_rmmissing[which((Gmerged_rmmissing$SampleDepth < 100 ) & (Gmerged_rmmissing$RefDepth < 100)), ]
JC1merged_rmOutliers <- JC1merged_rmmissing[which((JC1merged_rmmissing$SampleDepth < 100 ) & (JC1merged_rmmissing$RefDepth < 100)), ]
JC2merged_rmOutliers <- JC2merged_rmmissing[which((JC2merged_rmmissing$SampleDepth < 100 ) & (JC2merged_rmmissing$RefDepth < 100)), ]

#rm(Bmerged_rmmissing, Cmerged_rmmissing, Dmerged_rmmissing, Gmerged_rmmissing, JC1merged_rmmissing, JC2merged_rmmissing)

BnetDepth <- data.frame(Chrom = Bmerged_rmOutliers$SampleChrom, Base = Bmerged_rmOutliers$StartBase, NetDepth = (Bmerged_rmOutliers$SampleDepth - Bmerged_rmOutliers$RefDepth))
CnetDepth <- data.frame(Chrom = Cmerged_rmOutliers$SampleChrom, Base = Cmerged_rmOutliers$StartBase, NetDepth = (Cmerged_rmOutliers$SampleDepth - Cmerged_rmOutliers$RefDepth))
DnetDepth <- data.frame(Chrom = Dmerged_rmOutliers$SampleChrom, Base = Dmerged_rmOutliers$StartBase, NetDepth = (Dmerged_rmOutliers$SampleDepth - Dmerged_rmOutliers$RefDepth))
GnetDepth <- data.frame(Chrom = Gmerged_rmOutliers$SampleChrom, Base = Gmerged_rmOutliers$StartBase, NetDepth = (Gmerged_rmOutliers$SampleDepth - Gmerged_rmOutliers$RefDepth))
JC1netDepth <- data.frame(Chrom = JC1merged_rmOutliers$SampleChrom, Base = JC1merged_rmOutliers$StartBase, NetDepth = (JC1merged_rmOutliers$SampleDepth - JC1merged_rmOutliers$RefDepth))
JC2netDepth <- data.frame(Chrom = JC2merged_rmOutliers$SampleChrom, Base = JC2merged_rmOutliers$StartBase, NetDepth = (JC2merged_rmOutliers$SampleDepth - JC2merged_rmOutliers$RefDepth))

#BnetDepth <- data.frame(Chrom = Bmerged_rmNA$SampleChrom, Base = Bmerged_rmNA$StartBase, NetDepth = (as.numeric(Bmerged_rmNA$SampleDepth) - as.numeric(Bmerged_rmNA$RefDepth)))
#CnetDepth <- data.frame(Chrom = Cmerged_rmNA$SampleChrom, Base = Cmerged_rmNA$StartBase, NetDepth = (as.numeric(Cmerged_rmNA$SampleDepth) - as.numeric(Cmerged_rmNA$RefDepth)))
#DnetDepth <- data.frame(Chrom = Dmerged_rmNA$SampleChrom, Base = Dmerged_rmNA$StartBase, NetDepth = (as.numeric(Dmerged_rmNA$SampleDepth) - as.numeric(Dmerged_rmNA$RefDepth)))
#GnetDepth <- data.frame(Chrom = Gmerged_rmNA$SampleChrom, Base = Gmerged_rmNA$StartBase, NetDepth = (as.numeric(Gmerged_rmNA$SampleDepth) - as.numeric(Gmerged_rmNA$RefDepth)))
#JC1netDepth <- data.frame(Chrom = JC1merged_rmNA$SampleChrom, Base = JC1merged_rmNA$StartBase, NetDepth = (as.numeric(JC1merged_rmNA$SampleDepth) - as.numeric(JC1merged_rmNA$RefDepth)))
#JC2netDepth <- data.frame(Chrom = JC2merged_rmNA$SampleChrom, Base = JC2merged_rmNA$StartBase, NetDepth = (as.numeric(JC2merged_rmNA$SampleDepth) - as.numeric(JC2merged_rmNA$RefDepth)))

#BnetDepthNorm = data.frame(Chrom = BnetDepth$Chrom, Base = BnetDepth$Base, NormDepth = (BnetDepth$NetDepth-min(BnetDepth$NetDepth))/(max(BnetDepth$NetDepth)-min(BnetDepth$NetDepth)))
#CnetDepthNorm = data.frame(Chrom = CnetDepth$Chrom, Base = CnetDepth$Base, NormDepth = (CnetDepth$NetDepth-min(CnetDepth$NetDepth))/(max(CnetDepth$NetDepth)-min(CnetDepth$NetDepth)))
#DnetDepthNorm = data.frame(Chrom = DnetDepth$Chrom, Base = DnetDepth$Base, NormDepth = (DnetDepth$NetDepth-min(DnetDepth$NetDepth))/(max(DnetDepth$NetDepth)-min(DnetDepth$NetDepth)))
#GnetDepthNorm = data.frame(Chrom = GnetDepth$Chrom, Base = GnetDepth$Base, NormDepth = (GnetDepth$NetDepth-min(GnetDepth$NetDepth))/(max(GnetDepth$NetDepth)-min(GnetDepth$NetDepth)))
#JC1netDepthNorm = data.frame(Chrom = JC1netDepth$Chrom, Base = JC1netDepth$Base, NormDepth = (JC1netDepth$NetDepth-min(JC1netDepth$NetDepth))/(max(JC1netDepth$NetDepth)-min(JC1netDepth$NetDepth)))
#JC2netDepthNorm = data.frame(Chrom = JC2netDepth$Chrom, Base = JC2netDepth$Base, NormDepth = (JC2netDepth$NetDepth-min(JC2netDepth$NetDepth))/(max(JC2netDepth$NetDepth)-min(JC2netDepth$NetDepth)))

#BnetDepthScaff <- subset(BnetDepthNorm, BnetDepthNorm$Chrom %in% c("Ledus.S1","Ledus.S2","Ledus.S3","Ledus.S4","Ledus.S5","Ledus.S6"))
#CnetDepthScaff <- subset(CnetDepthNorm, CnetDepthNorm$Chrom %in% c("Ledus.S1","Ledus.S2","Ledus.S3","Ledus.S4","Ledus.S5","Ledus.S6"))
#DnetDepthScaff <- subset(DnetDepthNorm, DnetDepthNorm$Chrom %in% c("Ledus.S1","Ledus.S2","Ledus.S3","Ledus.S4","Ledus.S5","Ledus.S6"))
#GnetDepthScaff <- subset(GnetDepthNorm, GnetDepthNorm$Chrom %in% c("Ledus.S1","Ledus.S2","Ledus.S3","Ledus.S4","Ledus.S5","Ledus.S6"))
#JC1netDepthScaff <- subset(JC1netDepthNorm, JC1netDepthNorm$Chrom %in% c("Ledus.S1","Ledus.S2","Ledus.S3","Ledus.S4","Ledus.S5","Ledus.S6"))
#JC2netDepthScaff <- subset(JC2netDepthNorm, JC2netDepthNorm$Chrom %in% c("Ledus.S1","Ledus.S2","Ledus.S3","Ledus.S4","Ledus.S5","Ledus.S6"))

BnetDepthScaff <- subset(BnetDepth, BnetDepth$Chrom %in% c("Ledus.S1","Ledus.S2","Ledus.S3","Ledus.S4","Ledus.S5","Ledus.S6"))
CnetDepthScaff <- subset(CnetDepth, CnetDepth$Chrom %in% c("Ledus.S1","Ledus.S2","Ledus.S3","Ledus.S4","Ledus.S5","Ledus.S6"))
DnetDepthScaff <- subset(DnetDepth, DnetDepth$Chrom %in% c("Ledus.S1","Ledus.S2","Ledus.S3","Ledus.S4","Ledus.S5","Ledus.S6"))
GnetDepthScaff <- subset(GnetDepth, GnetDepth$Chrom %in% c("Ledus.S1","Ledus.S2","Ledus.S3","Ledus.S4","Ledus.S5","Ledus.S6"))
JC1netDepthScaff <- subset(JC1netDepth, JC1netDepth$Chrom %in% c("Ledus.S1","Ledus.S2","Ledus.S3","Ledus.S4","Ledus.S5","Ledus.S6"))
JC2netDepthScaff <- subset(JC2netDepth, JC2netDepth$Chrom %in% c("Ledus.S1","Ledus.S2","Ledus.S3","Ledus.S4","Ledus.S5","Ledus.S6"))

BnetDepthFinal <- cbind("B", BnetDepthScaff)
CnetDepthFinal <- cbind("C", CnetDepthScaff)
DnetDepthFinal <- cbind("D", DnetDepthScaff)
GnetDepthFinal <- cbind("G", GnetDepthScaff)
JC1netDepthFinal <- cbind("JC1", JC1netDepthScaff)
JC2netDepthFinal <- cbind("JC2", JC2netDepthScaff)

colnames(BnetDepthFinal) <- c("Sample", "Chrom", "Base", "NormDepth")
colnames(CnetDepthFinal) <- c("Sample", "Chrom", "Base", "NormDepth")
colnames(DnetDepthFinal) <- c("Sample", "Chrom", "Base", "NormDepth")
colnames(GnetDepthFinal) <- c("Sample", "Chrom", "Base", "NormDepth")
colnames(JC1netDepthFinal) <- c("Sample", "Chrom", "Base", "NormDepth")
colnames(JC2netDepthFinal) <- c("Sample", "Chrom", "Base", "NormDepth")

combined <- bind_rows(BnetDepthFinal, CnetDepthFinal, DnetDepthFinal, JC1netDepthFinal, JC2netDepthFinal)

ggplot(data = combined, aes(x = Chrom, y = NormDepth, fill = Sample)) +
  geom_boxplot(outliers = F, outlier.shape = ".", outlier.alpha = 0.5) +
  xlab("") +
  #ylim(-30,30)+
  ylab("Net Read Depth") +
  scale_fill_manual(values=c("#648FFF", "#AFC5FF", "lightblue", "#b95657", "#8c383a", "#4e0d14")) +
  theme_bw()

# Isolate just sex chromosome for line plot
BS6.netDepthFinal <- subset(BnetDepthFinal, BnetDepthFinal$Chrom %in% c("Ledus.S6"))
CS6.netDepthFinal <- subset(CnetDepthFinal, CnetDepthFinal$Chrom %in% c("Ledus.S6"))
DS6.netDepthFinal <- subset(DnetDepthFinal, DnetDepthFinal$Chrom %in% c("Ledus.S6"))
GS6.netDepthFinal <- subset(GnetDepthFinal, GnetDepthFinal$Chrom %in% c("Ledus.S6"))
JC1S6.netDepthFinal <- subset(JC1netDepthFinal, JC1netDepthFinal$Chrom %in% c("Ledus.S6"))
JC2S6.netDepthFinal <- subset(JC2netDepthFinal, JC2netDepthFinal$Chrom %in% c("Ledus.S6"))

combinedS6male <- bind_rows(BS6.netDepthFinal, CS6.netDepthFinal, DS6.netDepthFinal)
combinedS6female <- bind_rows(JC1S6.netDepthFinal, JC2S6.netDepthFinal)
combinedS6malefemale <- bind_rows(BS6.netDepthFinal, CS6.netDepthFinal, DS6.netDepthFinal, JC1S6.netDepthFinal, JC2S6.netDepthFinal)

ggplot(data = combinedS6malefemale, aes(x = Base, y = NormDepth, color = Sample)) +
  geom_point(alpha = 0.5, show.legend = F) +
  geom_line(show.legend = F, stat = "smooth", span = 0.5) +
  #xlab("") +
 # ylim(-25,5)+
  ylab("Net Read Depth") +
  scale_color_manual(values=c("#648FFF", "#AFC5FF", "lightblue", "#b95657", "#8c383a")) +
  #geom_smooth(method = loess, ) +
  theme_bw()

# Separate male and female plot
leduss6male <- ggplot(data = combinedS6male, aes(x = Base, y = NormDepth, color = Sample)) +
  geom_line(show.legend = T) +
  #xlab("") +
  ylim(-25,5)+
  ylab("Net Read Depth") +
  scale_color_manual(values=c("#648FFF", "#AFC5FF", "lightblue")) +
  theme_bw()

leduss6female <- ggplot(data = combinedS6female, aes(x = Base, y = NormDepth, color = Sample)) +
  geom_line(show.legend = T) +
  #xlab("") +
  ylim(-25,5)+
  ylab("Net Read Depth") +
  scale_color_manual(values=c("#b95657", "#8c383a", "#4e0d14")) +
  theme_bw()

ggarrange(leduss6male, leduss6female, ncol = 1)


#save(Bmerged_rmOutliers, BmergedSample, BnetDepth, BnetDepthFinal, BnetDepthNorm, BnetDepthScaff, 
#     Cmerged_rmOutliers, CmergedSample, CnetDepth, CnetDepthFinal, CnetDepthNorm, CnetDepthScaff,
#     Dmerged_rmOutliers, DmergedSample, DnetDepth, DnetDepthFinal, DnetDepthNorm, DnetDepthScaff,
#     Gmerged_rmOutliers, GmergedSample, GnetDepth, GnetDepthFinal, GnetDepthNorm, GnetDepthScaff,
#     JC1merged_rmOutliers, JC1mergedSample, JC1netDepth, JC1netDepthFinal, JC1netDepthNorm, JC1netDepthScaff,
#     JC2merged_rmOutliers, JC2mergedSample, JC2netDepth, JC2netDepthFinal, JC2netDepthNorm, JC2netDepthScaff,
#     file = "leioReadMapping.Rda")
#load("leioReadMapping.Rda")
#rm(Bmerged_rmOutliers, BmergedSample, BnetDepth, BnetDepthFinal, BnetDepthNorm, BnetDepthScaff, 
#     Cmerged_rmOutliers, CmergedSample, CnetDepth, CnetDepthFinal, CnetDepthNorm, CnetDepthScaff,
#     Dmerged_rmOutliers, DmergedSample, DnetDepth, DnetDepthFinal, DnetDepthNorm, DnetDepthScaff,
#     Gmerged_rmOutliers, GmergedSample, GnetDepth, GnetDepthFinal, GnetDepthNorm, GnetDepthScaff,
#     JC1merged_rmOutliers, JC1mergedSample, JC1netDepth, JC1netDepthFinal, JC1netDepthNorm, JC1netDepthScaff,
#     JC2merged_rmOutliers, JC2mergedSample, JC2netDepth, JC2netDepthFinal, JC2netDepthNorm, JC2netDepthScaff,
#     combined)

BnetDepthFinalSample <- sample_n(BnetDepthFinal, 1000)
CnetDepthFinalSample <- sample_n(CnetDepthFinal, 1000)
DnetDepthFinalSample <- sample_n(DnetDepthFinal, 1000)
GnetDepthFinalSample <- sample_n(GnetDepthFinal, 1000)
JC1netDepthFinalSample <- sample_n(JC1netDepthFinal, 1000)
JC2netDepthFinalSample <- sample_n(JC2netDepthFinal, 1000)

aovB <- aov(BnetDepthFinalSample$NormDepth ~ BnetDepthFinalSample$Chrom)
TukeyHSD(x=aovB, 'BnetDepthFinalSample$Chrom', conf.level=0.95)

aovC <- aov(CnetDepthFinalSample$NormDepth ~ CnetDepthFinalSample$Chrom)
TukeyHSD(x=aovC, 'CnetDepthFinalSample$Chrom', conf.level=0.95)

aovD <- aov(DnetDepthFinalSample$NormDepth ~ DnetDepthFinalSample$Chrom)
TukeyHSD(x=aovD, 'DnetDepthFinalSample$Chrom', conf.level=0.95)

aovG <- aov(GnetDepthFinalSample$NormDepth ~ GnetDepthFinalSample$Chrom)
TukeyHSD(x=aovG, 'GnetDepthFinalSample$Chrom', conf.level=0.95)

aovJC1 <- aov(JC1netDepthFinalSample$NormDepth ~ JC1netDepthFinalSample$Chrom)
TukeyHSD(x=aovJC1, 'JC1netDepthFinalSample$Chrom', conf.level=0.95)

aovJC2 <- aov(JC2netDepthFinalSample$NormDepth ~ JC2netDepthFinalSample$Chrom)
TukeyHSD(x=aovJC2, 'JC2netDepthFinalSample$Chrom', conf.level=0.95)


pairwise.wilcox.test(BnetDepthFinalSample$NormDepth, BnetDepthFinalSample$Chrom, p.adjust.method = "bonf")
pairwise.wilcox.test(CnetDepthFinalSample$NormDepth, CnetDepthFinalSample$Chrom, p.adjust.method = "bonf")
pairwise.wilcox.test(DnetDepthFinalSample$NormDepth, DnetDepthFinalSample$Chrom, p.adjust.method = "bonf")
pairwise.wilcox.test(JC1netDepthFinalSample$NormDepth, JC1netDepthFinalSample$Chrom, p.adjust.method = "bonf")
pairwise.wilcox.test(JC2netDepthFinalSample$NormDepth, JC2netDepthFinalSample$Chrom, p.adjust.method = "bonf")

#################################################################################################
# Phymatoceros
Onedepth <- read.delim("PhymatoCJR5469_1.genomeMapped.bam.depth.100b.bed", header = F, sep = "\t")
Twodepth <- read.delim("PhymatoCJR5469_2.genomeMapped.bam.depth.100b.bed", header = F, sep = "\t")
Threedepth <- read.delim("PhymatoCJR5469_3.genomeMapped.bam.depth.100b.bed", header = F, sep = "\t")
Fivedepth <- read.delim("PhymatoCJR5469_5.genomeMapped.bam.depth.100b.bed", header = F, sep = "\t")
Sixdepth <- read.delim("PhymatoCJR5469_6.genomeMapped.bam.depth.100b.bed", header = F, sep = "\t")

colnames(Onedepth) <- c("SampleChrom","SampleBase", "StopBase", "SampleDepth")
colnames(Twodepth) <- c("SampleChrom","SampleBase", "StopBase", "SampleDepth")
colnames(Threedepth) <- c("SampleChrom","SampleBase", "StopBase", "SampleDepth")
colnames(Fivedepth) <- c("SampleChrom","SampleBase", "StopBase", "SampleDepth")
colnames(Sixdepth) <- c("SampleChrom","SampleBase", "StopBase",  "SampleDepth")

Onedepth$SampleDepth <- as.numeric(Onedepth$SampleDepth)
Twodepth$SampleDepth <- as.numeric(Twodepth$SampleDepth)
Threedepth$SampleDepth <- as.numeric(Threedepth$SampleDepth)
Fivedepth$SampleDepth <- as.numeric(Fivedepth$SampleDepth)
Sixdepth$SampleDepth <- as.numeric(Sixdepth$SampleDepth)

HsubsetOnedepth <- read.delim("Phymatoceros_phymatodes_H40.subset1.genomeMapped.bam.depth.100b.bed", header = F, sep = "\t")
HsubsetTwodepth <- read.delim("Phymatoceros_phymatodes_H40.subset2.genomeMapped.bam.depth.100b.bed", header = F, sep = "\t")
HsubsetThreedepth <- read.delim("Phymatoceros_phymatodes_H40.subset3.genomeMapped.bam.depth.100b.bed", header = F, sep = "\t")
HsubsetFivedepth <- read.delim("Phymatoceros_phymatodes_H40.subset5.genomeMapped.bam.depth.100b.bed", header = F, sep = "\t")
HsubsetSixdepth <- read.delim("Phymatoceros_phymatodes_H40.subset6.genomeMapped.bam.depth.100b.bed", header = F, sep = "\t")

colnames(HsubsetOnedepth) <- c("RefChrom","RefBase", "RefStop", "RefDepth")
colnames(HsubsetTwodepth) <- c("RefChrom","RefBase", "RefStop", "RefDepth")
colnames(HsubsetThreedepth) <- c("RefChrom","RefBase", "RefStop", "RefDepth")
colnames(HsubsetFivedepth) <- c("RefChrom","RefBase", "RefStop", "RefDepth")
colnames(HsubsetSixdepth) <- c("RefChrom","RefBase",  "RefStop", "RefDepth")

HsubsetOnedepth$RefDepth <- as.numeric(HsubsetOnedepth$RefDepth)
HsubsetTwodepth$RefDepth <- as.numeric(HsubsetTwodepth$RefDepth)
HsubsetThreedepth$RefDepth <- as.numeric(HsubsetThreedepth$RefDepth)
HsubsetFivedepth$RefDepth <- as.numeric(HsubsetFivedepth$RefDepth)
HsubsetSixdepth$RefDepth <- as.numeric(HsubsetSixdepth$RefDepth)

Onemerged <- cbind(Onedepth, HsubsetOnedepth)
#OnemergedSample <- sample_n(Onemerged, 10000000)
#rm(Onedepth, HsubsetOnedepth, Onemerged)

Twomerged <- cbind(Twodepth, HsubsetTwodepth)
#TwomergedSample <- sample_n(Twomerged, 10000000)
#rm(Twodepth, HsubsetTwodepth, Twomerged)

Threemerged <- cbind(Threedepth, HsubsetThreedepth)
#ThreemergedSample <- sample_n(Threemerged, 10000000)
#rm(Threedepth, HsubsetThreedepth, Threemerged)

Fivemerged <- cbind(Fivedepth, HsubsetFivedepth)
#FivemergedSample <- sample_n(Fivemerged, 10000000)
#rm(Fivedepth, HsubsetFivedepth, Fivemerged)

Sixmerged <- cbind(Sixdepth, HsubsetSixdepth)
#SixmergedSample <- sample_n(Sixmerged, 10000000)
#rm(Sixdepth, HsubsetSixdepth, Sixmerged)


Onemerged_rmNA <- Onemerged[Onemerged$RefDepth != "NaN", ]
Twomerged_rmNA <- Twomerged[Twomerged$RefDepth != "NaN", ]
Threemerged_rmNA <- Threemerged[Threemerged$RefDepth != "NaN", ]
Fivemerged_rmNA <- Fivemerged[Fivemerged$RefDepth != "NaN", ]
Sixmerged_rmNA <- Sixmerged[Sixmerged$RefDepth != "NaN", ]

Onemerged_rmmissing <- Onemerged_rmNA[Onemerged_rmNA$RefDepth != 0.0, ]
Twomerged_rmmissing <- Twomerged_rmNA[Twomerged_rmNA$RefDepth != 0.0, ]
Threemerged_rmmissing <- Threemerged_rmNA[Threemerged_rmNA$RefDepth != 0.0, ]
Fivemerged_rmmissing <- Fivemerged_rmNA[Fivemerged_rmNA$RefDepth != 0.0, ]
Sixmerged_rmmissing <- Sixmerged_rmNA[Sixmerged_rmNA$RefDepth != 0.0, ]

Onemerged_rmOutliers <- Onemerged_rmmissing[which((Onemerged_rmmissing$SampleDepth < 100 ) & (Onemerged_rmmissing$RefDepth < 100)), ]
Twomerged_rmOutliers <- Twomerged_rmmissing[which((Twomerged_rmmissing$SampleDepth < 100 ) & (Twomerged_rmmissing$RefDepth < 100)), ]
Threemerged_rmOutliers <- Threemerged_rmmissing[which((Threemerged_rmmissing$SampleDepth < 100 ) & (Threemerged_rmmissing$RefDepth < 100)), ]
Fivemerged_rmOutliers <- Fivemerged_rmmissing[which((Fivemerged_rmmissing$SampleDepth < 100 ) & (Fivemerged_rmmissing$RefDepth < 100)), ]
Sixmerged_rmOutliers <- Sixmerged_rmmissing[which((Sixmerged_rmmissing$SampleDepth < 100 ) & (Sixmerged_rmmissing$RefDepth < 100)), ]

#rm(Onemerged_rmmissing, Twomerged_rmmissing, Threemerged_rmmissing, Fivemerged_rmmissing, Sixmerged_rmmissing)

OnenetDepth <- data.frame(Chrom = Onemerged_rmOutliers$SampleChrom, Base = Onemerged_rmOutliers$SampleBase, NetDepth = (Onemerged_rmOutliers$SampleDepth - Onemerged_rmOutliers$RefDepth))
TwonetDepth <- data.frame(Chrom = Twomerged_rmOutliers$SampleChrom, Base = Twomerged_rmOutliers$SampleBase, NetDepth = (Twomerged_rmOutliers$SampleDepth - Twomerged_rmOutliers$RefDepth))
ThreenetDepth <- data.frame(Chrom = Threemerged_rmOutliers$SampleChrom, Base = Threemerged_rmOutliers$SampleBase, NetDepth = (Threemerged_rmOutliers$SampleDepth - Threemerged_rmOutliers$RefDepth))
FivenetDepth <- data.frame(Chrom = Fivemerged_rmOutliers$SampleChrom, Base = Fivemerged_rmOutliers$SampleBase, NetDepth = (Fivemerged_rmOutliers$SampleDepth - Fivemerged_rmOutliers$RefDepth))
SixnetDepth <- data.frame(Chrom = Sixmerged_rmOutliers$SampleChrom, Base = Sixmerged_rmOutliers$SampleBase, NetDepth = (Sixmerged_rmOutliers$SampleDepth - Sixmerged_rmOutliers$RefDepth))

#OnenetDepthNorm = data.frame(Chrom = OnenetDepth$Chrom, Base = OnenetDepth$Base, NormDepth = (OnenetDepth$NetDepth-min(OnenetDepth$NetDepth))/(max(OnenetDepth$NetDepth)-min(OnenetDepth$NetDepth)))
#TwonetDepthNorm = data.frame(Chrom = TwonetDepth$Chrom, Base = TwonetDepth$Base, NormDepth = (TwonetDepth$NetDepth-min(TwonetDepth$NetDepth))/(max(TwonetDepth$NetDepth)-min(TwonetDepth$NetDepth)))
#ThreenetDepthNorm = data.frame(Chrom = ThreenetDepth$Chrom, Base = ThreenetDepth$Base, NormDepth = (ThreenetDepth$NetDepth-min(ThreenetDepth$NetDepth))/(max(ThreenetDepth$NetDepth)-min(ThreenetDepth$NetDepth)))
#FivenetDepthNorm = data.frame(Chrom = FivenetDepth$Chrom, Base = FivenetDepth$Base, NormDepth = (FivenetDepth$NetDepth-min(FivenetDepth$NetDepth))/(max(FivenetDepth$NetDepth)-min(FivenetDepth$NetDepth)))
#SixnetDepthNorm = data.frame(Chrom = SixnetDepth$Chrom, Base = SixnetDepth$Base, NormDepth = (SixnetDepth$NetDepth-min(SixnetDepth$NetDepth))/(max(SixnetDepth$NetDepth)-min(SixnetDepth$NetDepth)))

OnenetDepthScaff <- subset(OnenetDepth, OnenetDepth$Chrom %in% c("Phphy.S1","Phphy.S2","Phphy.S3","Phphy.S4","Phphy.S5"))
TwonetDepthScaff <- subset(TwonetDepth, TwonetDepth$Chrom %in% c("Phphy.S1","Phphy.S2","Phphy.S3","Phphy.S4","Phphy.S5"))
ThreenetDepthScaff <- subset(ThreenetDepth, ThreenetDepth$Chrom %in% c("Phphy.S1","Phphy.S2","Phphy.S3","Phphy.S4","Phphy.S5"))
FivenetDepthScaff <- subset(FivenetDepth, FivenetDepth$Chrom %in% c("Phphy.S1","Phphy.S2","Phphy.S3","Phphy.S4","Phphy.S5"))
SixnetDepthScaff <- subset(SixnetDepth, SixnetDepth$Chrom %in% c("Phphy.S1","Phphy.S2","Phphy.S3","Phphy.S4","Phphy.S5"))

OnenetDepthFinal <- cbind("1", OnenetDepthScaff)
TwonetDepthFinal <- cbind("2", TwonetDepthScaff)
ThreenetDepthFinal <- cbind("3", ThreenetDepthScaff)
FivenetDepthFinal <- cbind("5", FivenetDepthScaff)
SixnetDepthFinal <- cbind("6", SixnetDepthScaff)

colnames(OnenetDepthFinal) <- c("Sample", "Chrom", "Base", "NormDepth")
colnames(TwonetDepthFinal) <- c("Sample", "Chrom", "Base", "NormDepth")
colnames(ThreenetDepthFinal) <- c("Sample", "Chrom", "Base", "NormDepth")
colnames(FivenetDepthFinal) <- c("Sample", "Chrom", "Base", "NormDepth")
colnames(SixnetDepthFinal) <- c("Sample", "Chrom", "Base", "NormDepth")

combined <- bind_rows(OnenetDepthFinal, ThreenetDepthFinal, TwonetDepthFinal, FivenetDepthFinal, SixnetDepthFinal)
combined$Sample <- factor(combined$Sample, levels = c("1", "3", "2", "5", "6"))


ggplot(data = combined, aes(x = Chrom, y = NormDepth, fill = Sample)) +
  geom_boxplot(outliers = F, outlier.alpha = 0.1, outlier.shape = ".") +
  #ylim(0.35, 0.65) +
  xlab("") +
  scale_fill_manual(values=c("#648FFF", "#AFC5FF", "#b95657", "#8c383a", "#4e0d14")) +
  ylab("Net Read Depth") +
  theme_bw()

One.S5.netDepth <- subset(OnenetDepthFinal, OnenetDepthFinal$Chrom %in% "Phphy.S5")
Two.S5.netDepth <- subset(TwonetDepthFinal, TwonetDepthFinal$Chrom %in% "Phphy.S5")
Three.S5.netDepth <- subset(ThreenetDepthFinal, ThreenetDepthFinal$Chrom %in% "Phphy.S5")
Five.S5.netDepth <- subset(FivenetDepthFinal, FivenetDepthFinal$Chrom %in% "Phphy.S5")
Six.S5.netDepth <- subset(SixnetDepthFinal, SixnetDepthFinal$Chrom %in% "Phphy.S5")

combinedS5male <- bind_rows(One.S5.netDepth, Three.S5.netDepth)
combinedS5female <- bind_rows(Two.S5.netDepth, Five.S5.netDepth, Six.S5.netDepth)
combinedS5malefemale <- bind_rows(One.S5.netDepth, Three.S5.netDepth,Two.S5.netDepth, Five.S5.netDepth, Six.S5.netDepth )

ggplot(data = combinedS5malefemale, aes(x = Base, y = NormDepth, color = Sample)) +
  geom_point(alpha = 0.5, show.legend = F) +
  geom_line(show.legend = F, stat = "smooth", span = 0.5) +
  #xlab("") +
  #ylim(-55,5)+
  ylab("Net Read Depth") +
  scale_color_manual(values=c("#648FFF", "#b95657", "#AFC5FF", "#8c383a", "#4e0d14")) +
  theme_bw()

phphyS5male <- ggplot(data = combinedS5male, aes(x = Base, y = NormDepth, color = Sample)) +
  geom_line(show.legend = T) +
  #xlab("") +
  ylim(-30,12)+
  ylab("Net Read Depth") +
  scale_color_manual(values=c("#648FFF", "#AFC5FF")) +
  theme_bw()
phphyS5female <- ggplot(data = combinedS5female, aes(x = Base, y = NormDepth, color = Sample)) +
  geom_line(show.legend = T) +
  #xlab("") +
  ylim(-30,12)+
  ylab("Net Read Depth") +
  scale_color_manual(values=c("#b95657", "#8c383a", "#4e0d14")) +
  theme_bw()

ggarrange(phphyS5female, phphyS5male, ncol = 1)

save(Onemerged_rmOutliers, OnemergedSample, OnenetDepth, OnenetDepthFinal, OnenetDepthNorm, OnenetDepthScaff, 
     Twomerged_rmOutliers, TwomergedSample, TwonetDepth, TwonetDepthFinal, TwonetDepthNorm, TwonetDepthScaff,
     Threemerged_rmOutliers, ThreemergedSample, ThreenetDepth, ThreenetDepthFinal, ThreenetDepthNorm, ThreenetDepthScaff,
     Fivemerged_rmOutliers, FivemergedSample, FivenetDepth, FivenetDepthFinal, FivenetDepthNorm, FivenetDepthScaff,
     Sixmerged_rmOutliers, SixmergedSample, SixnetDepth, SixnetDepthFinal, SixnetDepthNorm, SixnetDepthScaff,
     file = "phymatocerosReadMapping.Rda")
load("phymatocerosReadMapping.Rda")
rm(Onemerged_rmOutliers, OnemergedSample, OnenetDepth, OnenetDepthFinal, OnenetDepthNorm, OnenetDepthScaff, 
   Twomerged_rmOutliers, TwomergedSample, TwonetDepth, TwonetDepthFinal, TwonetDepthNorm, TwonetDepthScaff,
   Threemerged_rmOutliers, ThreemergedSample, ThreenetDepth, ThreenetDepthFinal, ThreenetDepthNorm, ThreenetDepthScaff,
   Fivemerged_rmOutliers, FivemergedSample, FivenetDepth, FivenetDepthFinal, FivenetDepthNorm, FivenetDepthScaff,
   Sixmerged_rmOutliers, SixmergedSample, SixnetDepth, SixnetDepthFinal, SixnetDepthNorm, SixnetDepthScaff,
   combined)

OnenetDepthFinalSample <- sample_n(OnenetDepthFinal, 10000)
TwonetDepthFinalSample <- sample_n(TwonetDepthFinal, 10000)
ThreenetDepthFinalSample <- sample_n(ThreenetDepthFinal, 10000)
FivenetDepthFinalSample <- sample_n(FivenetDepthFinal, 10000)
SixnetDepthFinalSample <- sample_n(SixnetDepthFinal, 10000)

aov1 <- aov(OnenetDepthFinalSample$NormDepth ~ OnenetDepthFinalSample$Chrom)
TukeyHSD(x=aov1, 'OnenetDepthFinalSample$Chrom', conf.level=0.95)

aov2 <- aov(TwonetDepthFinalSample$NormDepth ~ TwonetDepthFinalSample$Chrom)
TukeyHSD(x=aov2, 'TwonetDepthFinalSample$Chrom', conf.level=0.95)

aov3 <- aov(ThreenetDepthFinalSample$NormDepth ~ ThreenetDepthFinalSample$Chrom)
TukeyHSD(x=aov3, 'ThreenetDepthFinalSample$Chrom', conf.level=0.95)

aov5 <- aov(FivenetDepthFinalSample$NormDepth ~ FivenetDepthFinalSample$Chrom)
TukeyHSD(x=aov5, 'FivenetDepthFinalSample$Chrom', conf.level=0.95)

aov6 <- aov(SixnetDepthFinalSample$NormDepth ~ SixnetDepthFinalSample$Chrom)
TukeyHSD(x=aov6, 'SixnetDepthFinalSample$Chrom', conf.level=0.95)


######################################################################################
# Phaeomegaceros
depth52 <- read.delim("Phaeom_14765_2.genomeMapped.bam.depth.100b.bed", header = F, sep = "\t")
depth54 <- read.delim("Phaeom_14765_4.genomeMapped.bam.depth.100b.bed", header = F, sep = "\t")
depth58 <- read.delim("Phaeom_14765_8.genomeMapped.bam.depth.100b.bed", header = F, sep = "\t")
depth62 <- read.delim("Phaeom_14766_2.genomeMapped.bam.depth.100b.bed", header = F, sep = "\t")
depth66 <- read.delim("Phaeom_14766_6.genomeMapped.bam.depth.100b.bed", header = F, sep = "\t")
depth67 <- read.delim("Phaeom_14766_7.genomeMapped.bam.depth.100b.bed", header = F, sep = "\t")

colnames(depth52) <- c("SampleChrom","StartBase","StopBase" ,"SampleDepth")
colnames(depth54) <- c("SampleChrom","StartBase","StopBase" ,"SampleDepth")
colnames(depth58) <- c("SampleChrom","StartBase","StopBase" ,"SampleDepth")
colnames(depth62) <- c("SampleChrom","StartBase","StopBase" ,"SampleDepth")
colnames(depth66) <- c("SampleChrom","StartBase","StopBase" ,"SampleDepth")
colnames(depth67) <- c("SampleChrom","StartBase","StopBase" ,"SampleDepth")

depth52$SampleDepth <- as.numeric(depth52$SampleDepth)
depth54$SampleDepth <- as.numeric(depth54$SampleDepth)
depth58$SampleDepth <- as.numeric(depth58$SampleDepth)
depth62$SampleDepth <- as.numeric(depth62$SampleDepth)
depth66$SampleDepth <- as.numeric(depth66$SampleDepth)
depth67$SampleDepth <- as.numeric(depth67$SampleDepth)

subsetdepth52 <- read.delim("Phaeomegaceros_14765-5.652.genomeMapped.bam.depth.100b.bed", header = F, sep = "\t")
subsetdepth54 <- read.delim("Phaeomegaceros_14765-5.654subset.genomeMapped.bam.depth.100b.bed", header = F, sep = "\t")
subsetdepth58 <- read.delim("Phaeomegaceros_14765-5.658subset.genomeMapped.bam.depth.100b.bed", header = F, sep = "\t")
subsetdepth62 <- read.delim("Phaeomegaceros_14765-5.662subset.genomeMapped.bam.depth.100b.bed", header = F, sep = "\t")
subsetdepth66 <- read.delim("Phaeomegaceros_14765-5.666subset.genomeMapped.bam.depth.100b.bed", header = F, sep = "\t")
subsetdepth67 <- read.delim("Phaeomegaceros_14765-5.667subset.genomeMapped.bam.depth.100b.bed", header = F, sep = "\t")

colnames(subsetdepth52) <- c("RefChrom","RefStart", "RefStop", "RefDepth")
colnames(subsetdepth54) <- c("RefChrom","RefStart", "RefStop", "RefDepth")
colnames(subsetdepth58) <- c("RefChrom","RefStart", "RefStop", "RefDepth")
colnames(subsetdepth62) <- c("RefChrom","RefStart", "RefStop", "RefDepth")
colnames(subsetdepth66) <- c("RefChrom","RefStart", "RefStop", "RefDepth")
colnames(subsetdepth67) <- c("RefChrom","RefStart", "RefStop", "RefDepth")

subsetdepth52$RefDepth <- as.numeric(subsetdepth52$RefDepth)
subsetdepth54$RefDepth <- as.numeric(subsetdepth54$RefDepth)
subsetdepth58$RefDepth <- as.numeric(subsetdepth58$RefDepth)
subsetdepth62$RefDepth <- as.numeric(subsetdepth62$RefDepth)
subsetdepth66$RefDepth <- as.numeric(subsetdepth66$RefDepth)
subsetdepth67$RefDepth <- as.numeric(subsetdepth67$RefDepth)

merged52 <- cbind(depth52, subsetdepth52)
#merged52Sample <- sample_n(merged52, 10000000)
#rm(depth52, subsetdepth52, merged52)

merged54 <- cbind(depth54, subsetdepth54)
#merged54Sample <- sample_n(merged54, 10000000)
#rm(depth54, subsetdepth54, merged54)

merged58 <- cbind(depth58, subsetdepth58)
#merged58Sample <- sample_n(merged58, 10000000)
#rm(depth58, subsetdepth58, merged58)

merged62 <- cbind(depth62, subsetdepth62)
#merged62Sample <- sample_n(merged62, 10000000)
#rm(depth62, subsetdepth62, merged62)

merged66 <- cbind(depth66, subsetdepth66)
#merged66Sample <- sample_n(merged66, 10000000)
#rm(depth66, subsetdepth66, merged66)

merged67 <- cbind(depth67, subsetdepth67)
#merged67Sample <- sample_n(merged67, 10000000)
#rm(depth67, subsetdepth67, merged67)

merged52_rmmissing <- merged52[merged52$RefDepth != "0", ]
merged54_rmmissing <- merged54[merged54$RefDepth != "0", ]
merged58_rmmissing <- merged58[merged58$RefDepth != "0", ]
merged62_rmmissing <- merged62[merged62$RefDepth != "0", ]
merged66_rmmissing <- merged66[merged66$RefDepth != "0", ]
merged67_rmmissing <- merged67[merged67$RefDepth != "0", ]

merged52_rmNA <- merged52_rmmissing[merged52_rmmissing$RefDepth != "NaN", ]
merged54_rmNA <- merged54_rmmissing[merged54_rmmissing$RefDepth != "NaN", ]
merged58_rmNA <- merged58_rmmissing[merged58_rmmissing$RefDepth != "NaN", ]
merged62_rmNA <- merged62_rmmissing[merged62_rmmissing$RefDepth != "NaN", ]
merged66_rmNA <- merged66_rmmissing[merged66_rmmissing$RefDepth != "NaN", ]
merged67_rmNA <- merged67_rmmissing[merged67_rmmissing$RefDepth != "NaN", ]

merged52_rmOutliers <- merged52_rmmissing[which((merged52_rmmissing$SampleDepth < 100 ) & (merged52_rmmissing$RefDepth < 100)), ]
merged54_rmOutliers <- merged54_rmmissing[which((merged54_rmmissing$SampleDepth < 100 ) & (merged54_rmmissing$RefDepth < 100)), ]
merged58_rmOutliers <- merged58_rmmissing[which((merged58_rmmissing$SampleDepth < 100 ) & (merged58_rmmissing$RefDepth < 100)), ]
merged62_rmOutliers <- merged62_rmmissing[which((merged62_rmmissing$SampleDepth < 100 ) & (merged62_rmmissing$RefDepth < 100)), ]
merged66_rmOutliers <- merged66_rmmissing[which((merged66_rmmissing$SampleDepth < 100 ) & (merged66_rmmissing$RefDepth < 100)), ]
merged67_rmOutliers <- merged67_rmmissing[which((merged67_rmmissing$SampleDepth < 100 ) & (merged67_rmmissing$RefDepth < 100)), ]

#rm(merged52_rmmissing, merged54_rmmissing, merged58_rmmissing, merged62_rmmissing, merged66_rmmissing, merged67_rmmissing)

netDepth52 <- data.frame(Chrom = merged52_rmOutliers$SampleChrom, Base = merged52_rmOutliers$StartBase, NetDepth = (merged52_rmOutliers$SampleDepth - merged52_rmOutliers$RefDepth))
netDepth54 <- data.frame(Chrom = merged54_rmOutliers$SampleChrom, Base = merged54_rmOutliers$StartBase, NetDepth = (merged54_rmOutliers$SampleDepth - merged54_rmOutliers$RefDepth))
netDepth58 <- data.frame(Chrom = merged58_rmOutliers$SampleChrom, Base = merged58_rmOutliers$StartBase, NetDepth = (merged58_rmOutliers$SampleDepth - merged58_rmOutliers$RefDepth))
netDepth62 <- data.frame(Chrom = merged62_rmOutliers$SampleChrom, Base = merged62_rmOutliers$StartBase, NetDepth = (merged62_rmOutliers$SampleDepth - merged62_rmOutliers$RefDepth))
netDepth66 <- data.frame(Chrom = merged66_rmOutliers$SampleChrom, Base = merged66_rmOutliers$StartBase, NetDepth = (merged66_rmOutliers$SampleDepth - merged66_rmOutliers$RefDepth))
netDepth67 <- data.frame(Chrom = merged67_rmOutliers$SampleChrom, Base = merged67_rmOutliers$StartBase, NetDepth = (merged67_rmOutliers$SampleDepth - merged67_rmOutliers$RefDepth))

#netDepth52 <- data.frame(Chrom = merged52_rmOutliers$SampleChrom, Base = merged52_rmOutliers$SampleBase, NetDepth = (merged52_rmOutliers$SampleDepth - merged52_rmOutliers$RefDepth))
#netDepth54 <- data.frame(Chrom = merged54_rmOutliers$SampleChrom, Base = merged54_rmOutliers$SampleBase, NetDepth = (merged54_rmOutliers$SampleDepth - merged54_rmOutliers$RefDepth))
#netDepth58 <- data.frame(Chrom = merged58_rmOutliers$SampleChrom, Base = merged58_rmOutliers$SampleBase, NetDepth = (merged58_rmOutliers$SampleDepth - merged58_rmOutliers$RefDepth))
#netDepth62 <- data.frame(Chrom = merged62_rmOutliers$SampleChrom, Base = merged62_rmOutliers$SampleBase, NetDepth = (merged62_rmOutliers$SampleDepth - merged62_rmOutliers$RefDepth))
#netDepth66 <- data.frame(Chrom = merged66_rmOutliers$SampleChrom, Base = merged66_rmOutliers$SampleBase, NetDepth = (merged66_rmOutliers$SampleDepth - merged66_rmOutliers$RefDepth))
#netDepth67 <- data.frame(Chrom = merged67_rmOutliers$SampleChrom, Base = merged67_rmOutliers$SampleBase, NetDepth = (merged67_rmOutliers$SampleDepth - merged67_rmOutliers$RefDepth))

#netDepth52Norm = data.frame(Chrom = netDepth52$Chrom, Base = netDepth52$Base, NormDepth = (netDepth52$NetDepth-min(netDepth52$NetDepth))/(max(netDepth52$NetDepth)-min(netDepth52$NetDepth)))
#netDepth54Norm = data.frame(Chrom = netDepth54$Chrom, Base = netDepth54$Base, NormDepth = (netDepth54$NetDepth-min(netDepth54$NetDepth))/(max(netDepth54$NetDepth)-min(netDepth54$NetDepth)))
#netDepth58Norm = data.frame(Chrom = netDepth58$Chrom, Base = netDepth58$Base, NormDepth = (netDepth58$NetDepth-min(netDepth58$NetDepth))/(max(netDepth58$NetDepth)-min(netDepth58$NetDepth)))
#netDepth62Norm = data.frame(Chrom = netDepth62$Chrom, Base = netDepth62$Base, NormDepth = (netDepth62$NetDepth-min(netDepth62$NetDepth))/(max(netDepth62$NetDepth)-min(netDepth62$NetDepth)))
#netDepth66Norm = data.frame(Chrom = netDepth66$Chrom, Base = netDepth66$Base, NormDepth = (netDepth66$NetDepth-min(netDepth66$NetDepth))/(max(netDepth66$NetDepth)-min(netDepth66$NetDepth)))
#netDepth67Norm = data.frame(Chrom = netDepth67$Chrom, Base = netDepth67$Base, NormDepth = (netDepth67$NetDepth-min(netDepth67$NetDepth))/(max(netDepth67$NetDepth)-min(netDepth67$NetDepth)))

netDepth52Scaff <- subset(netDepth52, netDepth52$Chrom %in% c("Phfim.S1","Phfim.S2","Phfim.S3","Phfim.S4","Phfim.S5", "Phfim.S6"))
netDepth54Scaff <- subset(netDepth54, netDepth54$Chrom %in% c("Phfim.S1","Phfim.S2","Phfim.S3","Phfim.S4","Phfim.S5", "Phfim.S6"))
netDepth58Scaff <- subset(netDepth58, netDepth58$Chrom %in% c("Phfim.S1","Phfim.S2","Phfim.S3","Phfim.S4","Phfim.S5", "Phfim.S6"))
netDepth62Scaff <- subset(netDepth62, netDepth62$Chrom %in% c("Phfim.S1","Phfim.S2","Phfim.S3","Phfim.S4","Phfim.S5", "Phfim.S6"))
netDepth66Scaff <- subset(netDepth66, netDepth66$Chrom %in% c("Phfim.S1","Phfim.S2","Phfim.S3","Phfim.S4","Phfim.S5", "Phfim.S6"))
netDepth67Scaff <- subset(netDepth67, netDepth67$Chrom %in% c("Phfim.S1","Phfim.S2","Phfim.S3","Phfim.S4","Phfim.S5", "Phfim.S6"))

netDepth52Final <- cbind("14765-2", netDepth52Scaff)
netDepth54Final <- cbind("14765-4", netDepth54Scaff)
netDepth58Final <- cbind("14765-8", netDepth58Scaff)
netDepth62Final <- cbind("14766-2", netDepth62Scaff)
netDepth66Final <- cbind("14766-6", netDepth66Scaff)
netDepth67Final <- cbind("14766-7", netDepth67Scaff)

colnames(netDepth52Final) <- c("Sample", "Chrom", "Base", "NormDepth")
colnames(netDepth54Final) <- c("Sample", "Chrom", "Base", "NormDepth")
colnames(netDepth58Final) <- c("Sample", "Chrom", "Base", "NormDepth")
colnames(netDepth62Final) <- c("Sample", "Chrom", "Base", "NormDepth")
colnames(netDepth66Final) <- c("Sample", "Chrom", "Base", "NormDepth")
colnames(netDepth67Final) <- c("Sample", "Chrom", "Base", "NormDepth")

combined <- bind_rows(netDepth52Final, netDepth54Final, netDepth58Final, netDepth62Final, netDepth66Final, netDepth67Final)

ggplot(data = combined, aes(x = Chrom, y = NormDepth, fill = Sample)) +
  geom_boxplot(outliers = F, outlier.alpha = 0.1, outlier.shape = ".") +
  #ylim(-30, 30) +
  xlab("") +
  ylab("Net Read Depth") +
  scale_fill_manual(values=c("#ff8585", "#b95657", "lightblue", "#AFC5FF", "#648FFF", "blue", "darkblue")) +
  theme_bw()


S2netDepth52Scaff <- subset(netDepth52Final, netDepth52Final$Chrom %in% c("Phfim.S2"))
S2netDepth54Scaff <- subset(netDepth54Final, netDepth54Final$Chrom %in% c("Phfim.S2"))
S2netDepth58Scaff <- subset(netDepth58Final, netDepth58Final$Chrom %in% c("Phfim.S2"))
S2netDepth62Scaff <- subset(netDepth62Final, netDepth62Final$Chrom %in% c("Phfim.S2"))
S2netDepth66Scaff <- subset(netDepth66Final, netDepth66Final$Chrom %in% c("Phfim.S2"))
S2netDepth67Scaff <- subset(netDepth67Final, netDepth67Final$Chrom %in% c("Phfim.S2"))

S2netDepth52Final <- cbind("1", S2netDepth52Scaff)
S2netDepth54Final <- cbind("1", S2netDepth54Scaff)
S2netDepth58Final <- cbind("2", S2netDepth58Scaff)
S2netDepth62Final <- cbind("2", S2netDepth62Scaff)
S2netDepth66Final <- cbind("2", S2netDepth66Scaff)
S2netDepth67Final <- cbind("2", S2netDepth67Scaff)

colnames(S2netDepth52Final) <- c("Sex","Sample", "Chrom", "Base", "NetDepth")
colnames(S2netDepth54Final) <- c("Sex","Sample", "Chrom", "Base", "NetDepth")
colnames(S2netDepth58Final) <- c("Sex","Sample", "Chrom", "Base", "NetDepth")
colnames(S2netDepth62Final) <- c("Sex","Sample", "Chrom", "Base", "NetDepth")
colnames(S2netDepth66Final) <- c("Sex","Sample", "Chrom", "Base", "NetDepth")
colnames(S2netDepth67Final) <- c("Sex","Sample", "Chrom", "Base", "NetDepth")

combinedS2male <- bind_rows(S2netDepth52Final, S2netDepth54Final )
combinedS2female <- bind_rows(S2netDepth58Final, S2netDepth62Final, S2netDepth66Final, S2netDepth67Final)
combinedS2malefemale<- bind_rows(S2netDepth52Final, S2netDepth54Final, S2netDepth58Final, S2netDepth62Final, S2netDepth66Final, S2netDepth67Final)

ggplot(data = combinedS2malefemale, aes(x = Base, y = NetDepth, color = Sample)) +
  geom_point(alpha = 0.5, show.legend = F) +
  geom_line(show.legend = F, stat = "smooth", span = 0.5) +
  #xlab("") +
  #ylim(-55,5)+
  ylab("Net Read Depth") +
  scale_color_manual(values=c("#ff8585", "#b95657","lightblue", "#AFC5FF", "#648FFF", "blue", "darkblue")) +
  theme_bw()

phfimS2male <- ggplot(data = combinedS2male, aes(x = Base, y = NetDepth, color = Sample)) +
  geom_line(show.legend = T) +
  #xlab("") +
  ylim(-55,5)+
  ylab("Net Read Depth") +
  scale_color_manual(values=c("#648FFF", "#AFC5FF")) +
  theme_bw()
phfimS2female <- ggplot(data = combinedS2female, aes(x = Base, y = NetDepth, color = Sample)) +
  geom_line(show.legend = T) +
  #xlab("") +
  ylim(-55,5)+
  ylab("Net Read Depth") +
  scale_color_manual(values=c("#ff8585", "#b95657", "#8c383a", "#4e0d14")) +
  theme_bw()

ggarrange(phfimS2female, phfimS2male, ncol = 1)

save(merged52_rmOutliers, merged52Sample, netDepth52, netDepth52Final, netDepth52Norm, netDepth52Scaff, 
     merged54_rmOutliers, merged54Sample, netDepth54, netDepth54Final, netDepth54Norm, netDepth54Scaff,
     merged58_rmOutliers, merged58Sample, netDepth58, netDepth58Final, netDepth58Norm, netDepth58Scaff,
     merged62_rmOutliers, merged62Sample, netDepth62, netDepth62Final, netDepth62Norm, netDepth62Scaff,
     merged66_rmOutliers, merged66Sample, netDepth66, netDepth66Final, netDepth66Norm, netDepth66Scaff,
     merged67_rmOutliers, merged67Sample, netDepth67, netDepth67Final, netDepth67Norm, netDepth67Scaff,
     file = "phaeomegacerosReadMapping.Rda")
load("phaeomegacerosReadMapping.Rda")
rm(merged52_rmOutliers, merged52Sample, netDepth52, netDepth52Final, netDepth52Norm, netDepth52Scaff, 
   merged54_rmOutliers, merged54Sample, netDepth54, netDepth54Final, netDepth54Norm, netDepth54Scaff,
   merged58_rmOutliers, merged58Sample, netDepth58, netDepth58Final, netDepth58Norm, netDepth58Scaff,
   merged62_rmOutliers, merged62Sample, netDepth62, netDepth62Final, netDepth62Norm, netDepth62Scaff,
   merged66_rmOutliers, merged66Sample, netDepth66, netDepth66Final, netDepth66Norm, netDepth66Scaff,
   merged67_rmOutliers, merged67Sample, netDepth67, netDepth67Final, netDepth67Norm, netDepth67Scaff,
   combined)

netDepth52FinalSample <- sample_n(netDepth52Final, 10000)
netDepth54FinalSample <- sample_n(netDepth54Final, 10000)
netDepth58FinalSample <- sample_n(netDepth58Final, 10000)
netDepth62FinalSample <- sample_n(netDepth62Final, 10000)
netDepth66FinalSample <- sample_n(netDepth66Final, 10000)
netDepth67FinalSample <- sample_n(netDepth67Final, 10000)

combinedFinalSample <- bind_rows(netDepth52FinalSample, netDepth54FinalSample, netDepth58FinalSample, netDepth62FinalSample, netDepth66FinalSample, netDepth67FinalSample)

ggplot(data = combinedFinalSample, aes(x = Chrom, y = NormDepth, fill = Sample)) +
  geom_boxplot(outlier.alpha = 0.1, outlier.shape = ".") +
  ylim(-30, 30) +
  xlab("") +
  ylab("Net Read Depth") +
  scale_fill_manual(values=c("#648FFF", "#AFC5FF", "#ff8585", "#b95657", "#8c383a", "#4e0d14")) +
  theme_bw()


aov52 <- aov(netDepth52FinalSample$NormDepth ~ netDepth52FinalSample$Chrom)
TukeyHSD(x=aov52, 'netDepth52FinalSample$Chrom', conf.level=0.95)

aov54 <- aov(netDepth54FinalSample$NormDepth ~ netDepth54FinalSample$Chrom)
TukeyHSD(x=aov54, 'netDepth54FinalSample$Chrom', conf.level=0.95)

aov58 <- aov(netDepth58FinalSample$NormDepth ~ netDepth58FinalSample$Chrom)
TukeyHSD(x=aov58, 'netDepth58FinalSample$Chrom', conf.level=0.95)

aov62 <- aov(netDepth62FinalSample$NormDepth ~ netDepth62FinalSample$Chrom)
TukeyHSD(x=aov62, 'netDepth62FinalSample$Chrom', conf.level=0.95)

aov66 <- aov(netDepth66FinalSample$NormDepth ~ netDepth66FinalSample$Chrom)
TukeyHSD(x=aov66, 'netDepth66FinalSample$Chrom', conf.level=0.95)

aov67 <- aov(netDepth67FinalSample$NormDepth ~ netDepth67FinalSample$Chrom)
TukeyHSD(x=aov67, 'netDepth67FinalSample$Chrom', conf.level=0.95)

kruskal.test(netDepth52FinalSample$NormDepth ~ netDepth52FinalSample$Chrom)
pairwise.wilcox.test(netDepth52FinalSample$NormDepth, netDepth52FinalSample$Chrom, p.adjust.method = "bonf")

kruskal.test(netDepth54FinalSample$NormDepth ~ netDepth54FinalSample$Chrom)
pairwise.wilcox.test(netDepth54FinalSample$NormDepth, netDepth54FinalSample$Chrom, p.adjust.method = "bonf")

kruskal.test(netDepth58FinalSample$NormDepth ~ netDepth58FinalSample$Chrom)
pairwise.wilcox.test(netDepth58FinalSample$NormDepth, netDepth58FinalSample$Chrom, p.adjust.method = "bonf")

kruskal.test(netDepth62FinalSample$NormDepth ~ netDepth62FinalSample$Chrom)
pairwise.wilcox.test(netDepth62FinalSample$NormDepth, netDepth62FinalSample$Chrom, p.adjust.method = "bonf")

kruskal.test(netDepth66FinalSample$NormDepth ~ netDepth66FinalSample$Chrom)
pairwise.wilcox.test(netDepth66FinalSample$NormDepth, netDepth66FinalSample$Chrom, p.adjust.method = "bonf")

