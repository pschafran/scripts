#### Isa's plot
utg83 <- read.delim(file.choose(), sep = "\t", header = F)
allbut83 <- read.delim(file.choose(), sep = "\t", header = F)

# Subset just 2nd and 3rd columns of utg83
utg83_2 <- utg83[, 2:3]
# Add "Contig 83" in new column
utg83_2 <- cbind(utg83_2, "Contig 83")

# Add "Other contigs" in new column
allbut83_2 <- cbind(allbut83,"Other contigs")

# Add column names to each 
colnames(allbut83_2) <- c("Position","Depth","Source")
colnames(utg83_2) <- c("Position","Depth","Source")

# Combine dataframes to long format
longformat <- rbind(utg83_2,allbut83_2)

# Make plot
depthPlot <- ggplot(data = longformat, aes(x=Position,y=Depth,color=Source))+
  geom_point(size=0.01, alpha=0.002)+
  geom_smooth(level = 0.99, se = TRUE ) +
  xlim(0,2.5e+06)+
  ylim(0,40)

depthPlot + theme_bw() + theme(text = element_text(size=20))