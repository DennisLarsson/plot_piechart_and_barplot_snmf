#source("/home/biogeoanalysis/RAD/spicatumGroup/06populations_50miss_mac2/PieChartMap-1.2/piechartMap.R")
library(rworldmap)
library(plotrix)

# version 1.21

ck=require("LEA")
if (ck==FALSE) {
  install.packages("devtools")
  devtools::install_github("bcm-uga/LEA")
}
library(parallel)

#give path to csv file
# version 1.1

setwd("/home/biogeoanalysis/RAD/spicatumGroup/06populations_50miss_mac2/snmf/")    #give the path to the work directory where you files are and/or where you would like to output files to go
inputfile.name = "spicgrp_K4_snmf_popAverage_offset.csv"
outputfile.name="spicgrp_K4_snmf"

print.label = F
pieChartSize = 0.6

inputfile <- read.csv(inputfile.name,header=FALSE, as.is=TRUE)
chart.data <- inputfile[4:ncol(inputfile)]
plot.location <- inputfile[2:3]
labels <- inputfile[1]

actual.location = read.csv("spicgrp_K4_snmf_popAverage.csv",header=FALSE, as.is=TRUE)
actual.location <- actual.location[2:3]

blwt.color <- colorRampPalette(c("grey10", "white"))

plot.piechart <- function(actual.loc, plot.loc, chart.data, labels, lab.con, pCS) {
  #plots map and axes
  #worldmap <- getMap(resolution = "low")
  #plot(worldmap, xlim = c(min(plot.loc[,1]), max(plot.loc[,1])), ylim = c(min(plot.loc[,2]), max(plot.loc[,2])))
  #box(which="plot")
  #axis(1)
  #axis(2)
  elev_map <- raster("/home/biogeoanalysis/iDDC/wc2.1_30s_elev.tif")
  range = as(extent(min(plot.loc[,1])-1, max(plot.loc[,1])+1, min(plot.loc[,2])-1, max(plot.loc[,2])+1),'SpatialPolygons') # W E S N
  crs(range) <- "+proj=longlat +datum=WGS84 +no_defs"
  elev_eur <- crop(elev_map, range)
  image(elev_eur, col = blwt.color(100))
  bordCol = "black"
  
  #plots the piecharts, admixData$AverK3_1 etc can be changed and extended if needed, example: admixData$AverK4_1[x],admixData$AverK4_2[x],admixData$AverK4_3[x],admixData$AverK4_4[x]
  #Just make sure the column names in your input csv matches those in the for loop!
  col.list <- c("yellow", "blue", "blue", "blue", "red", "red", "red", "red", "green", "green", "green", "green", "green", "green", "green", 
                "orange", "orange", "orange", "orange", "orange", "orange", "purple", "purple")
  
  for (x in 1:nrow(chart.data)) {
    arrows(plot.loc[x,1], plot.loc[x,2], actual.loc[x,1], actual.loc[x,2], length = 0.07, col = "white", code=0, lwd = 1.5)
    floating.pie(plot.loc[x,1], plot.loc[x,2], unlist(chart.data[x,]), 
                 radius=pCS, col=c("red", "blue", "orange", "green","purple", "brown","darkgrey", "yellow", "darkgreen", "cyan"), border = NA)
    
    draw.circle(plot.loc[x,1], plot.loc[x,2],radius = pCS, border=col.list[x], lwd = 3)
    
  }
  if (lab.con) {
    #function to put labels in map
    for (x in 1:nrow(actual.loc)) {
      text(plot.loc[x,1], plot.loc[x,2], labels=labels[x,1], cex= 0.7, pos = 2, offset = 1)
    }
  }
}

# version 1.21

setwd("/home/biogeoanalysis/RAD/spicatumGroup/06populations_50miss_mac2/snmf")  #give the path to the work directory where you files are and/or where you would like to output files to go
filename="/home/biogeoanalysis/RAD/spicatumGroup/06populations_50miss_mac2/populations.snps.vcf"           #give the name of the file if it is in the work directory or the full path if it is somewhere else
outputname="spicgrp"             #give the prefix for all of the output files
popmap="/home/biogeoanalysis/RAD/spicatumGroup/popmap_spicGroup_sub_Grp"      #give the name of the popmap if it is in the work directory or the full path if it is somewhere else

NrCores=NULL # NULL can be replaced with a number of cores if you don't want to use all.
if (is.null(NrCores)) {NrCores=detectCores()} #detects how many cores are available

pop <- read.delim(popmap, header = FALSE)
pop_sorted<-pop[order(pop[,2]),]

#vcf2geno(filename, output.file = paste(outputname,".geno",sep=""), force = TRUE) #converts the vcf to geno format 

#obj.snmf <- snmf(paste(outputname,".geno",sep=""), K = 1:10, repetitions = 100, project = "new", CPU = NrCores, entropy = T, iterations = 2000) #runs the actual analysis

#if a projects has already been run, you can load it with the below command. Just uncomment it (remove the starting '#') and change the path to the file
obj.snmf <- load.snmfProject("/home/biogeoanalysis/RAD/spicatumGroup/06populations_50miss_mac2/snmf/spicgrp.snmfProject") 

#this set up a list of cross entropy for each run and K.
CE_AllK <- list()
for (K in 1:10) {
  CE <- vector()
  x <- 1
  for (i in seq(K,length(obj.snmf@runs),10)){
    CE[x] <- obj.snmf@runs[[i]]@crossEntropy
    x=x+1
  }
  CE_AllK[[K]] <- CE
}

#this calculates the mean cross entropy
mean.ce = vector()
for (i in 1:10) {
  mean.ce[i] <- mean(unlist(CE_AllK[i]))
}
diff.mean.ce <- vector()
i=1
while (i <= length(mean.ce)-1){
  diff.mean.ce[i] <- mean.ce[i]-mean.ce[i+1]
  i = i +1
}

plot(obj.snmf, pch = 19, col = "blue",main = "Cross-entropy", xaxt="n")
axis(side=1, at=1:10, labels = seq(1,10))

#sets up a loop to plot the admixture ratios for each K
K=4
bestRun <- match(min(unlist(CE_AllK[K])),unlist(CE_AllK[K])) #find the run number of the run with the lowest cross entropy

qmatrix = Q(obj.snmf, K = K, run=bestRun) #extracts the best (lowest cross entropy) qmatrix (matrix of admixture ratios) from the full dataset.
colorsPlot = c("red", "blue", "orange", "green","purple", "brown","darkgrey", "yellow", "darkgreen", "cyan")
barplot(t(qmatrix), border = NA, space = 0, ylab = "Ancestry coefficients", col = colorsPlot, main = paste("Ancestry coefficients for K=",K,sep = ""))

#plot lines and names of populations into the plot
axis(1, tapply(1:nrow(pop), pop[,2], mean), unique(pop_sorted[,2]), las=2, cex.axis=1, tick = F, line = -0.8)
abline(v=tapply(1:nrow(pop), pop[,2],max), lty=2, lwd=0.5)

#write the qmatrix to a csv file that can be used in other applications.
write.csv(qmatrix,paste(outputname,"_K",K,"_snmf.csv",sep=""), row.names = pop$V1)


png(file="snmf_combined.png", height = 700, width = 700, pointsize = 20)
mat = c(1,2,3,3,3,3)
layout(matrix(mat, nrow = 3, ncol = 2, byrow = TRUE))
#plot 1
par(mar = c(2, 3, 2, 1))
plot(obj.snmf, pch = 19, col = "blue",main = "Cross-entropy", xaxt="n", las =2)
axis(side=1, at=1:10, labels = seq(1,10))

#plot 2
par(mar = c(2, 2.5, 2, 0))
barplot(t(qmatrix), border = NA, space = 0, ylab = "Ancestry coefficients", col = colorsPlot, main = paste("Ancestry coefficients for K=",K,sep = ""), las =2)
#axis(1, tapply(1:nrow(pop), pop[,2], mean), unique(pop_sorted[,2]), las=2, cex.axis=1, tick = F, line = -0.8)
text(x = tapply(1:nrow(pop), pop[,2], mean), y = -0.04, labels = unique(pop_sorted[,2]), xpd = NA, srt = 55, adj = 1, cex = 1)
abline(v=tapply(1:nrow(pop), pop[,2],max), lty=1, lwd=1)

#plot 3
par(mar = c(2, 2.5, 1.5, 0.5))
plot.piechart(actual.location, plot.location, chart.data, labels, print.label, pieChartSize)
dev.off()



pdf(file="snmf_combined.pdf", height = 10, width = 10, pointsize = 16)
mat = c(1,2,3,3,3,3)
layout(matrix(mat, nrow = 3, ncol = 2, byrow = TRUE))
#plot 1
par(mar = c(2, 3, 2, 1))
plot(obj.snmf, pch = 19, col = "blue",main = "Cross-entropy", xaxt="n", las =2)
axis(side=1, at=1:10, labels = seq(1,10))

#plot 2
par(mar = c(2, 2.5, 2, 0))
barplot(t(qmatrix), border = NA, space = 0, ylab = "Ancestry coefficients", col = colorsPlot, main = paste("Ancestry coefficients for K=",K,sep = ""), las =2)
#axis(1, tapply(1:nrow(pop), pop[,2], mean), unique(pop_sorted[,2]), las=2, cex.axis=1, tick = F, line = -0.8)
text(x = tapply(1:nrow(pop), pop[,2], mean), y = -0.04, labels = unique(pop_sorted[,2]), xpd = NA, srt = 55, adj = 1, cex = 1)
abline(v=tapply(1:nrow(pop), pop[,2],max), lty=1, lwd=1)

#plot 3
par(mar = c(2, 2, 1.5, 0.5))
plot.piechart(actual.location, plot.location, chart.data, labels, print.label, pieChartSize)
dev.off()
