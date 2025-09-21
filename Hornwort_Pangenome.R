library(UpSetR)
library(viridis)
library(RColorBrewer)
library(wesanderson)
library(ggplot2)
library(ggpubr)

paDF <- read.delim("/Volumes/WD14TB-2/projects/hornwort_genomes/analysis/orthofinder/hornworts-only/OrthoFinder/Results_Apr12/Orthogroups/Orthogroups.PresenceAbsence.tsv", header = T, sep = "\t")

colnames(paDF) <- c("Orthogroup","AnagrBONN", "Anang","AnagrOXF", "Anfus", "Anpun", "Ledus", "Mefla", "Noorb", "Papea", "Phcar", "Phsp", "Phchi", "Phphy")
upset(paDF, nsets = 14, nintersects = 50, order.by = "freq", keep.order = T, sets = rev(c("Ledus", "AnagrOXF", "Anfus", "Anpun", "Noorb", "Papea", "Phcar", "Phsp", "Mefla", "Phchi", "Phphy")))
upset(paDF, nsets = 14, nintersects = 12, order.by = "degree", keep.order = T, sets = rev(c("Ledus", "AnagrOXF", "Anfus", "Anpun", "Noorb", "Papea", "Phcar", "Phsp", "Mefla", "Phchi", "Phphy")))


coreDF <- read.delim("/Volumes/WD14TB-2/projects/hornwort_genomes/analysis/orthofinder/hornworts-only/OrthoFinder/Results_Apr12/Orthogroups/Orthogroups.GeneCount.CoreGenes.no_Anang_AnagrBONN.tsv", header = T, sep ="\t")
periphDF <- read.delim("/Volumes/WD14TB-2/projects/hornwort_genomes/analysis/orthofinder/hornworts-only/OrthoFinder/Results_Apr12/Orthogroups/Orthogroups.GeneCount.PeripheralGenes.no_Anang_AnagrBONN.tsv", header = T, sep ="\t")
privateDF <- read.delim("/Volumes/WD14TB-2/projects/hornwort_genomes/analysis/orthofinder/hornworts-only/OrthoFinder/Results_Apr12/Orthogroups/Orthogroups.GeneCount.PrivateGenes.no_Anang_AnagrBONN.tsv", header = T, sep ="\t")


species <- c(rep("AnagrOXF",4), rep("Anfus", 4), rep("Anpun", 4), rep("Ledus", 4), rep("Mefla", 4), rep("Noorb", 4), rep("Papea", 4), rep("Phcar", 4), rep("Phsp", 4), rep("Phchi", 4), rep("Phphy", 4))
geneType <- c(rep(c("Core", "Peripheral", "Dispensible", "Private"), 11))
perc <- as.numeric(c(45,10,23,22,41,10,38,11,43,9,43,4,61,8,20,11,44,9,30,18,43,8,34,15,35,8,44,14,38,9,46,7,33,8,50,8,55,12,29,5,37,9,36,18))
abs <- as.numeric(c(9478,2147,4761,4716,8720,2027,8025,2330,9249,1995,9173,863,8030,994,2662,1489,8949,1859,6066,3637,9146,1760,7299,3194,9393,2097,11781,3794,9652,2273,11580,1837,9555,2374,14262,2411,8798,1850,4580,787,9661,2243,9432,4695))
geneCountsDF<- data.frame(species,geneType,perc,abs)
colnames(geneCountsDF) <- c("Species","Category","Perc","Count")
geneCountsDF$Species <- factor(geneCountsDF$Species, levels = c("Phphy","Phchi","Mefla","Phsp", "Phcar", "Papea", "Noorb", "Anpun", "Anfus", "AnagrOXF", "Ledus"))
geneCountsDF$Category <- factor(geneCountsDF$Category, levels = c("Private", "Dispensible", "Peripheral", "Core"))

percPlot <-  ggplot(data = geneCountsDF, aes(x = Species, y = Count, fill = Category)) +
  geom_bar(position='fill', stat="identity", show.legend = F) +
  coord_flip() +
  theme_bw() +
  labs(x = "Proportion") +
  scale_fill_manual(values = wes_palette("Moonrise3", type = 'discrete')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_blank())

absPlot <- ggplot(data = geneCountsDF, aes(x = Species, y = Count, fill = Category)) +
  geom_bar(position='stack', stat="identity", show.legend = T) +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = wes_palette("Moonrise3", type = 'discrete')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())

ggarrange(percPlot, absPlot, nrow=1)

sum(coreDF$Anthoceros_agrestis_Oxford_PROT_primary_transcripts)
sum(coreDF$Anthoceros_fusiformis_PROT_primary_transcripts)
sum(coreDF$Anthoceros_punctatus_PROT_primary_transcripts)
sum(coreDF$Leiosporoceros_dussii_PROT_primary_transcripts)
sum(coreDF$Megaceros_flagellaris_PROT_primary_transcripts)
sum(coreDF$Notothylas_orbicularis_PROT_primary_transcripts)
sum(coreDF$Paraphymatoceros_pearsonii_PROT_primary_transcripts)
sum(coreDF$Phaeoceros_carolinianus_PROT_primary_transcripts)
sum(coreDF$Phaeoceros_sp_PROT_primary_transcripts)
sum(coreDF$Phaeomegaceros_fimbriatus_PROT_primary_transcripts)
sum(coreDF$Phymatoceros_phymatodes_PROT_primary_transcripts)

sum(periphDF$Anthoceros_agrestis_Oxford_PROT_primary_transcripts)
sum(periphDF$Anthoceros_fusiformis_PROT_primary_transcripts)
sum(periphDF$Anthoceros_punctatus_PROT_primary_transcripts)
sum(periphDF$Leiosporoceros_dussii_PROT_primary_transcripts)
sum(periphDF$Megaceros_flagellaris_PROT_primary_transcripts)
sum(periphDF$Notothylas_orbicularis_PROT_primary_transcripts)
sum(periphDF$Paraphymatoceros_pearsonii_PROT_primary_transcripts)
sum(periphDF$Phaeoceros_carolinianus_PROT_primary_transcripts)
sum(periphDF$Phaeoceros_sp_PROT_primary_transcripts)
sum(periphDF$Phaeomegaceros_fimbriatus_PROT_primary_transcripts)
sum(periphDF$Phymatoceros_phymatodes_PROT_primary_transcripts)

sum(privateDF$Anthoceros_agrestis_Oxford_PROT_primary_transcripts)
sum(privateDF$Anthoceros_fusiformis_PROT_primary_transcripts)
sum(privateDF$Anthoceros_punctatus_PROT_primary_transcripts)
sum(privateDF$Leiosporoceros_dussii_PROT_primary_transcripts)
sum(privateDF$Megaceros_flagellaris_PROT_primary_transcripts)
sum(privateDF$Notothylas_orbicularis_PROT_primary_transcripts)
sum(privateDF$Paraphymatoceros_pearsonii_PROT_primary_transcripts)
sum(privateDF$Phaeoceros_carolinianus_PROT_primary_transcripts)
sum(privateDF$Phaeoceros_sp_PROT_primary_transcripts)
sum(privateDF$Phaeomegaceros_fimbriatus_PROT_primary_transcripts)
sum(privateDF$Phymatoceros_phymatodes_PROT_primary_transcripts)

upset(paDF, nsets = 14, nintersects = 50, order.by = "freq", keep.order = T, sets = rev(c("Anang", "Phphy")))


### 

divTime<-c(32.08,10.82,296.9,257,257,257,257,257,257,257,32.08,296.9,257,257,257,257,257,257,257,296.9,257,257,257,257,257,257,257,296.9,296.9,296.9,296.9,296.9,296.9,296.9,141,141,141,141,54.9,99.4,97.3,97.3,97.3,141,141,53.3,53.3,141,141,4.8,141,141,141,141,99.4)
sharedOGs<-c(10903,12251,8489,9094,9029,9427,9447,9545,9197,9430,10445,8464,9032,8955,9390,9359,9439,9149,9306,8359,8943,8855,9183,9192,9262,9051,9232,8330,8143,8419,8395,8438,8355,8339,9064,9451,9442,9531,9603,9472,9625,9683,9799,9125,9333,10640,10792,9557,9956,12809,9565,9917,9607,10029,9574)
clades<-c(1,1,3,3,3,3,3,3,3,3,1,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
 
 
divOGDF <- data.frame(divTime, sharedOGs, as.character(clades))
colnames(divOGDF) <- c("DivTime","OGs","Clades")
ggplot(data = divOGDF, aes(x = DivTime, y = OGs)) +
  geom_point() +
  #scale_color_manual(values = c("red", "blue", "black")) +
  theme_bw() +
  xlab("Divergence Time (Ma)") +
  ylab("Shared Orthogroups") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_smooth(method='lm', formula = y~log(x), color = "black") +
#  geom_smooth(method="lm", formula = y~I(1/x), color = "red")
  annotate(geom = "text", label="R-squared:  0.82", x =270, y = 13000, size= 3) +
  scale_x_log10()
  #scale_y_continuous(trans='log10', breaks = c(8000,9000,10000,11000,12000,13000))

lm <- lm(divOGDF$OGs ~ log(divOGDF$DivTime))
summary(lm)
lm.rep <- lm(divOGDF$OGs ~ I(1/divOGDF$DivTime))
summary(lm.rep)
mean(c(397,998,346,216,988,1500,1188,676,1935,540,1275))


### Syntenic Block Size vs Div. Time
divTime2<-c(3.44,32.08,10.82,296.9,257,257,257,257,257,257,32.08,10.82,296.9,257,257,257,257,257,257,32.08,296.9,257,257,257,257,257,257,296.9,257,257,257,257,257,257,296.9,296.9,296.9,296.9,296.9,296.9,141,141,141,54.9,99.4,97.3,97.3,141,141,53.3,141,141,141,141,99.4)
meanSynBlock<-c(183.1154,42.75814,57.90055,11.42133,10.63231,14.16245,14.49216,14.77303,12.74571,13.81356,41.67901,59.41117,11.56296,10.71543,14.16554,14.33621,14.63363,12.42308,13.90506,45.79812,11.99476,10.85128,14.27891,14.7791,15.43533,12.1525,14.66667,11.7602,10.71392,14.10345,14.03478,14.79128,12.09719,13.93038,9.963277,11.18612,12.16867,11.7659,10.47475,11.77287,12.35075,13.22646,13.28111,17.92678,14.09567,29.46058,29.81356,15.03191,17.12378,43.99519,16.17918,19.86687,16.28713,18.64516,19.41737)
sharedOGs2<-c(13839,10181,11219,8009,8583,8513,8910,8879,8688,8886,10903,12251,8489,9094,9029,9427,9447,9197,9430,10445,8464,9032,8955,9390,9359,9149,9306,8359,8943,8855,9183,9192,9051,9232,8330,8143,8419,8395,8355,8339,9064,9451,9442,9603,9472,9625,9683,9125,9333,10640,9557,9956,9565,9917,9574)
synBlockSum<-c(14283,9193,10480,4283,3817,3923,4623,4491,4461,4075,10128,11704,4683,4029,4193,4989,4873,4845,4394,9755,4582,4232,4198,4951,4893,4861,4444,4610,4157,4090,4842,4748,4730,4402,3527,3546,4040,4071,4148,3732,4965,5899,5764,8569,6188,7100,7036,5652,5257,9151,6682,6417,6580,6358,6932)
propSynGenes <-c(53.99,39.17,44.49,21.96,16.47,16.61,17.48,17.55,21.32,15.71,42.05,48.42,23.27,16.94,17.30,18.43,18.60,22.49,16.55,46.03,26.74,20.34,19.75,20.56,21.07,26.19,18.86,26.76,19.89,19.17,20.03,20.37,25.37,18.61,20.94,20.51,20.08,21.14,28.42,19.04,23.69,24.80,25.14,46.92,26.59,29.30,30.11,30.21,22.17,34.92,31.02,24.17,31.82,24.75,32.97)


divSynOGdf <- data.frame(divTime2, 
                         meanSynBlock, 
                         synBlockSum,
                         propSynGenes)
hwtSynDF <- cbind(divSynOGdf, "Hornworts")

colnames(hwtSynDF) <- c("DivTime", "SynBlockSize", "SynBlockSum", "PercSynGenes", "Phylum")

hwtNumSynGenes <- ggplot(data = hwtSynDF, aes(x = DivTime, y = synBlockSum)) +
  geom_point() +
  #scale_color_manual(values = c("red", "blue", "black")) +
  theme_bw() +
  xlab("Divergence Time (Ma)") +
  ylab("No. Syntenic Genes") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_smooth(method='lm', formula = y~log(x), color = "red") +
  ylim(0,25000) +
  xlim(0,300) +
  #geom_smooth(method='nls', formula = y~a*exp(-b*x), color = "black", method.args = list(start = st), se = FALSE) +
  #geom_smooth(method="lm", formula = y~I(1/x), color = "black") +
  annotate(geom = "text", label="R-squared:  0.96", x =270, y = 14000, size= 3)

lm <- lm(hwtSynDF$SynBlockSum ~ log(hwtSynDF$DivTime))
summary(lm)

fm0 <- lm(log(SynBlockSize) ~ DivTime, hwtSynDF)
st <- list(a = exp(coef(fm0)[[1]]), b = -coef(fm0)[[2]])
m.exp <- nls(SynBlockSize ~ a * exp(-b * DivTime), data = hwtSynDF, start = st)
summary(m.exp)

scatterplot(1/SynBlockSize ~ DivTime, data = hwtSynDF)
fm1 <- lm((1/SynBlockSize)~DivTime, data = hwtSynDF)
summary(fm1)


hwtAvgSynBlock <- ggplot(data = hwtSynDF, aes(y = meanSynBlock, x = DivTime)) +
  geom_point() +
  #scale_color_manual(values = c("red", "blue", "black")) +
  theme_bw() +
  xlab("Divergence Time (Ma") +
  ylab("Mean Syntenic Block Size") +
  ylim(0,375) +
  xlim(0,300) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_smooth(method='lm', formula = y~log(x), color = "red")

ggplot(data = hwtSynDF, aes(x = DivTime, y = PercSynGenes)) +
  geom_point() +
  #scale_color_manual(values = c("red", "blue", "black")) +
  theme_bw() +
  xlab("Divergence Time (Ma)") +
  ylab("% Syntenic Genes") +
  geom_smooth(method='lm', formula = y~log(x), color = "red") +
  #geom_smooth(method='nls', formula = y~a*exp(-b*x), color = "black", method.args = list(start = st), se = FALSE) +
  #geom_smooth(method="lm", formula = y~I(1/x), color = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

propSynGenes_lm <- lm(PropSynGenes ~ log(DivTime), data = hwtSynDF)
summary(lm)

### Angiosperm
angOldDivTime <- c(247,247,247,247,247,247,247,247,247,247,247,247,247,247,247,247,247,245.7,242.9,242.9,242.9,242.9,241.8,241.8,242.9,242.9,242.9,242.9,242.9,242.9,241.8,241.8,241.8,245.7,245.7,245.7,245.7,245.7,245.7,245.7,245.7,245.7,245.7,245.7,245.7,245.7,245.7,245.7,236.4,236.4,230.2,242.9,242.9,236.4,236.4,236.4,230.3,236.4,236.4,242.9,242.9,242.9,56.3,236.4,242.9,242.9,231.7,33.5,231.7,236.4,231.7,231.7,242.9,242.9,242.9,236.4,242.9,242.9,231.7,56.3,231.7,236.4,231.7,231.7,242.9,242.9,242.9,242.9,242.9,236.4,236.4,236.4,158.2,236.4,236.4,242.9,242.9,242.9,100.3,242.9,242.9,242.9,242.9,242.9,242.9,96.1,96.1,100.3,242.9,242.9,242.9,242.9,242.9,242.9,100.3,100.3,61.2,231.7,114.2,236.4,114.2,114.2,242.9,242.9,242.9,231.7,236.4,231.7,231.7,242.9,242.9,242.9,236.4,98.3,98.3,242.9,242.9,242.9,236.4,236.4,242.9,242.9,242.9,37.3,242.9,242.9,242.9,242.9,242.9,242.9,64.9,100.3,100.3)
angYoungDivTime <- c(154,154,154,154,154,154,154,154,154,154,154,154,154,154,154,154,154,153.4,152.1,152.1,152.1,152.1,151.7,151.7,152.1,152.1,152.1,152.1,152.1,152.1,151.7,151.7,151.7,153.4,153.4,153.4,153.4,153.4,153.4,153.4,153.4,153.4,153.4,153.4,153.4,153.4,153.4,153.4,148.6,148.6,144.9,152.1,152.1,148.6,148.6,148.6,144.9,148.6,148.6,152.1,152.1,152.1,34.1,148.6,152.1,152.1,146.3,19.2,146.3,148.6,146.3,146.3,152.1,152.1,152.1,148.6,152.1,152.1,146.3,34.1,146.3,148.6,146.3,146.3,152.1,152.1,152.1,152.1,152.1,148.6,148.6,148.6,97.4,148.6,148.6,152.1,152.1,152.1,68,152.1,152.1,152.1,152.1,152.1,152.1,66,66,68,152.1,152.1,152.1,152.1,152.1,152.1,68,68,39,146.3,80.5,148.6,80.5,80.5,152.1,152.1,152.1,146.3,148.6,146.3,146.3,152.1,152.1,152.1,148.6,71.3,71.3,152.1,152.1,152.1,148.6,148.6,152.1,152.1,152.1,29.4,152.1,152.1,152.1,152.1,152.1,152.1,42.5,68,68)
angMeanSynBlockSize <- c(13.10902,13.06545,9.535135,8.5,7.939394,8.297872,8.883562,7.428571,8.987654,8.009259,8.404762,6.921053,9.357798,7.969231,8.933884,9.351351,8.065359,10.13158,9.5,7.662791,7.055556,8.04329,8.547529,7.150538,8.425656,7.792135,7.949853,6.673913,8.676923,7.386861,8.517241,8.318182,7.898113,8.6139,7.867925,6.738095,7.367742,8.547872,7.807018,7.878136,7.663866,7.666667,6.548387,8.206186,7.367521,8.278481,8.138889,8.045918,11.48458,8.97205,12.54382,7.768519,7.741935,12.93266,11.9,13.23288,7.712177,15.40323,10.25397,7.366337,7.776316,7.556604,45.16095,8.984848,7.434783,6.615385,9.548632,355.5185,9.784703,7.198473,11.20367,7.72381,6.952381,6.192308,6.425532,7.90942,6.307692,5.75,8.174603,48.32764,7.997868,6.851064,9.136126,7.382353,6.307692,5.75,6.25,6.880597,8.5,9.781701,9.240618,10.19092,9.460477,11.92917,8.044898,7.086957,7.090909,6.886364,41.64659,7.186441,7.170732,7.103448,6.75,7.552846,6.526316,67.87097,55.90777,91.4434,7.125,7.727273,6.2,6,7.347826,6,31.61076,27.77113,38.65147,9.709627,24.58333,7.125604,29.47817,11.22977,7.10465,6.891892,7.06383,10.06065,7.056338,11.85,8.041096,6.969697,6.931034,6.54902,7.508772,44.51892,11.63565,6.8875,6.835821,6.922078,7.795652,6.945455,6.333333,6.142857,7.428571,16.08486,7.284091,7.328571,7.028037,6.421053,7.1,6.214286,52.32773,51.24157,41.70468)
angNumSynBlocks <-c(269,272,185,99,39,137,153,53,403,287,178,55,240,257,96,260,101,41,147,188,57,455,320,502,103,28,388,455,46,12,282,11,3,62,9,256)
angSumSynBlocks <- c(3487,3593,1764,816,262,1170,1297,416,2184,865,2118,263,2040,1036,1081,1038,1234,4235,2774,1318,381,1858,2248,665,2890,1387,2695,307,2820,1012,1976,1647,2093,2231,834,283,1142,1607,445,2198,912,2070,203,2388,862,1308,1172,1577,5214,2889,6297,839,240,7682,5355,7728,2090,6685,3230,744,591,801,17116,4151,342,86,6283,19198,6908,943,6106,1622,292,161,302,2183,82,23,3090,16963,3751,322,3490,502,82,23,50,461,85,6094,4186,6512,7540,5726,1971,326,234,303,10370,848,294,618,54,929,124,14728,11517,19386,171,85,124,6,169,12,9989,7887,14417,6253,19470,1475,14857,8602,611,510,664,7133,1002,6162,1761,230,201,334,1712,16472,10379,551,458,533,1793,382,38,43,52,12321,641,513,752,122,71,87,12454,18242,14263)
angPercSynGenes <- c(13.81947885,15.41233244,5.958050461,3.257550051,0.830981002,3.806363459,4.005868273,1.537807515,4.87010815,3.442923101,5.388627401,0.703123956,5.554043017,3.31472268,3.801251846,3.55522066,2.399751079,16.30162824,8.595287155,4.755375956,1.114181691,5.562124863,6.414792832,2.237701057,6.082737863,4.990914161,6.42102379,0.766140101,7.157996269,2.983402612,6.352778537,5.169004802,3.869584108,7.350058478,3.233059389,0.876826076,3.627181629,4.851467214,1.600834592,4.821074104,3.525250768,5.168345755,0.532096144,6.371993116,2.693665823,4.481831109,3.914103463,3.022896959,16.24779919,7.490277418,16.66799015,2.128442229,0.703967148,14.80553521,16.64853101,16.67457817,4.702388318,15.27266912,8.434411354,2.097015136,1.630907209,1.370097327,50.32267549,12.49492046,0.981038983,0.291179956,13.27529924,69.53907453,16.53086375,2.364119535,15.57116809,4.807635307,0.944326763,0.508207071,0.560239679,5.49860205,0.198352705,0.063863166,5.7426405,49.76383959,7.771194166,0.694451933,7.637931412,1.248212843,0.219245475,0.060273326,0.08280202,1.136882082,0.241316167,11.49442632,12.57208073,13.71611517,16.54361349,12.75221594,4.999175681,0.890467085,0.62619586,0.508440447,28.13118845,1.551508055,0.841550858,1.258232977,0.114368011,1.996068025,0.301952954,38.5050785,29.52471288,31.65914083,0.34664153,0.287070028,0.28316644,0.014323228,0.410039912,0.033575825,30.34003068,23.41606793,25.78723785,13.19114824,31.61535464,2.471369809,25.17751529,16.06844312,1.204724254,0.990762596,0.900937572,17.03891264,2.507350641,15.68417838,5.208133087,0.74203123,0.632981152,0.618747684,3.161967734,30.80663562,21.62584517,1.219647166,0.997050212,0.781972095,3.476928745,0.828759248,0.087807471,0.097649597,0.078478128,27.12771228,1.504624196,1.183104438,1.146586162,0.328606252,0.187409265,0.144733449,35.51335244,31.83928509,24.56914)
angAvgNumGenes <- c(25232.5,23312.5,29607,25049.5,31529,30738,32377.5,27051.5,44845,25124,39305,37404.5,36730,31254.5,28438,29196.5,51422,25979,32273.5,27716,34195.5,33404.5,35044,29718,47511.5,27790.5,41971.5,40071,39396.5,33921,31104.5,31863,54088.5,30353.5,25796,32275.5,31484.5,33124,27798,45591.5,25870.5,40051.5,38151,37476.5,32001,29184.5,29943,52168.5,32090.5,38570,37779,39418.5,34092.5,51886,32165,46346,44445.5,43771,38295.5,35479,36237.5,58463,34012.5,33221.5,34861,29535,47328.5,27607.5,41788.5,39888,39213.5,33738,30921.5,31680,53905.5,39701,41340.5,36014.5,53808,34087,48268,46367.5,45693,40217.5,37401,38159.5,60385,40549.5,35223.5,53017,33296,47477,45576.5,44902,39426.5,36610,37368.5,59594,36863,54656.5,34935.5,49116.5,47216,46541.5,41066,38249.5,39008,61233.5,49330.5,29609.5,43790.5,41890,41215.5,35740,32923.5,33682,55907.5,47403,61584,59683.5,59009,53533.5,50717,51475.5,73701,41863,39962.5,39288,33812.5,30996,31754.5,53980,54143.5,53469,47993.5,45177,45935.5,68161,51568.5,46093,43276.5,44035,66260.5,45418.5,42602,43360.5,65586,37126.5,37885,60110.5,35068.5,57294,58052.5)

ggplot(data = NULL, aes(x = angAvgNumGenes, y = angPercSynGenes)) +
  geom_point()
ggplot(data = NULL, aes(x = angSumSynBlocks, y = angPercSynGenes)) +
  geom_point() +
  geom_smooth(method="lm")
ggplot(data = NULL, aes(x = angSumSynBlocks, y = angAvgNumGenes)) +
  geom_point()


angSynDF <- data.frame(angOldDivTime, angMeanSynBlockSize, angSumSynBlocks, angPercSynGenes)
angSynDF <- cbind(angSynDF, "Angiosperm")
colnames(angSynDF) <- c("DivTime", "SynBlockSize", "SynBlockSum", "PercSynGenes", "Phylum")

tmp1 <- ggplot(data = angSynDF, aes(x = YoungDivTime, y = MeanSynBlockSize)) +
  geom_point() +
  #scale_color_manual(values = c("red", "blue", "black")) +
  xlab("Divergence Time (Ma)") +
  ylab("Mean Syntenic Block Size") +
  theme_bw() +
  geom_smooth(method='lm', formula = y~log(x), color = "red") +
  ylim(0,375) +
  xlim(0,300) +
  #geom_smooth(method='nls', formula = y~a*exp(-b*x), color = "black", method.args = list(start = st), se = FALSE) +
  geom_smooth(method="lm", formula = y~I(1/x), color = "black") +
  #annotate(geom = "text", label="R-squared:  0.96", x =270, y = 14000, size= 3)
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

tmp2 <- ggplot(data = angSynDF, aes(x = YoungDivTime, y = NumSynGenes)) +
  geom_point() +
  #scale_color_manual(values = c("red", "blue", "black")) +
  xlab("Divergence Time (Ma)") +
  ylab("No. Syntenic Genes") +
  geom_smooth(method='lm', formula = y~log(x), color = "red") +
  ylim(0,25000) +
  xlim(0,300) +
  #geom_smooth(method='nls', formula = y~a*exp(-b*x), color = "black", method.args = list(start = st), se = FALSE) +
  geom_smooth(method="lm", formula = y~I(1/x), color = "black") +
  #annotate(geom = "text", label="R-squared:  0.96", x =270, y = 14000, size= 3)
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggarrange(tmp1, tmp2, hwtAvgSynBlock, hwtNumSynGenes, nrow = 2, ncol = 2)

loglm <- lm(MeanSynBlockSize ~ log(YoungDivTime), data = angSynDF)
summary(loglm)
replm<- lm(I(1/MeanSynBlockSize)~YoungDivTime, data = angSynDF)
summary(replm)

#### Mosses
mossDivTime <- c(420,276,403,420,420,403)
mossMeanSynBlockSize <- c(6.428571,7.933718,6,6.333333,6.8,7)
mossSynBlockSum <- c(45,2753,12,38,34,7)
mossPercSynGenes <- c(0.16,8.69,0.04,0.13,0.13,0.02)

mossSynDF <- data.frame(mossDivTime, mossMeanSynBlockSize, mossSynBlockSum, mossPercSynGenes)
mossSynDF <- cbind(mossSynDF, "Mosses")
colnames(mossSynDF) <- c("DivTime", "SynBlockSize", "SynBlockSum", "PercSynGenes", "Phylum")


###### Combined DFs
ang_hwt_synDF <- rbind(angSynDF, hwtSynDF)
ang_hwt_moss_synDF <- rbind(ang_hwt_synDF, mossSynDF)
combined_synDF <- rbind(ang_hwt_moss_synDF, data.frame("DivTime" = 303, "SynBlockSize" = 8.351111, "SynBlockSum" = 1879, "PercSynGenes" = 9.96, "Phylum" = "Liverworts"))


blocksize <- ggplot(data = ang_hwt_moss_synDF, aes(x = DivTime, y = SynBlockSize, color = Phylum))+
  geom_point() +
  geom_smooth(method='lm', formula = y~I(1/x)) +
  xlab("Divergence Time (Ma)") +
  ylab("Avg. Synentic Block Size") +
  theme_bw()

blocksum <- ggplot(data = ang_hwt_moss_synDF, aes(x = DivTime, y = SynBlockSum, color = Phylum))+
  geom_point() +
  geom_smooth(method='lm', formula = y~log(x)) +
  xlab("Divergence Time (Ma)") +
  ylab("Total Sytentic Genes") +  
  theme_bw()

ggplot(data = combined_synDF, aes(x = DivTime, y = PercSynGenes, color = Phylum))+
  geom_point(show.legend = F) +
  geom_smooth(method='lm', formula = y~log(x), show.legend = F) +
  xlab("Divergence Time (Ma)") +
  ylab("% genes in syntenic blocks") +  
  theme_bw() +
  scale_color_manual(values = c("orange", "green4", "red", "blue")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggarrange(blocksize, blocksum, blockperc, nrow = 3)

### PPRs

speciesPPR <- c("AnagrBONN", "AnagrOXF","Anang","Anfus","Anpun","Ledus","Mefla","Noorb","Papea","Phcar","Phsp","Phchi","Phphy")
countsPPR <- as.numeric(c(1887, 397,918,998,346,216,988,1500,1188,676,1935,540,1275))
PPRdf <- data.frame(speciesPPR, countsPPR)
colnames(PPRdf) <- c("Species","Counts")
ggplot(data = PPRdf, aes(x = Species, y = Counts)) +
  geom_bar(stat="identity") +
  theme_bw() +
  xlab("") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())





