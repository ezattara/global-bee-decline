####HEADER####
# Global Bee Decline 
# An analysis of multidecadal trends in bee species richness
# from GBIF data
# Eduardo E. Zattara & Marcelo A. Aizen
# Grupo de Ecología de la Polinización, INIBIOMA, 
# Universidad Nacional del Comahue-CONICET, 
# Quintral 1250, Bariloche (8400), Argentina.
# Correspondence to ezattara@comahue-conicet.gob.ar

#Original data: https://doi.org/10.15468/dl.h52qyh 
#License: CC BY-NC 4.0 
#File: 415 MB CSV (Zipped)
#Involved datasets: 1026
#Filters: Scientific Name = Hymenoptera AND Basis of Record = Preserved Specimen

####Library installation and loading####
#install.packages("viridis", "iNEXT","data.table","dplyr","vegan", "tigerstats","RcolorBrewer")  # Uncomment to install packages

#Load libraries
library(iNEXT)
library("data.table")
library(dplyr)
library("viridis")
library("vegan")
library("tigerstats")
library(RColorBrewer)

####Read and format data####
#Select fields for GBIF data
fields<- c("gbifID","speciesKey","family","genus","species","countryCode","decimalLatitude","decimalLongitude","year","institutionCode","datasetKey","taxonRank","basisOfRecord")

#Select Hymenopteran families for the analysis
outgroup.fams <- c("Formicidae","Sphecidae","Crabronidae")
anthophila.fams <-c("Melittidae","Andrenidae","Halictidae","Colletidae","Megachilidae","Apidae")
all.fams <- c(outgroup.fams, anthophila.fams)
all.famcols=c("gray","blue","blue","violet","red","red","red","darkorange","darkorange")

#Select working directory where the following data files are found:
# "0001373-190621201848488.csv" # GBIF Download DOI: 10.15468/dl.h52qyh 
# "country-and-continent-codes-list.csv" # Country and continent codes, from https://datahub.io/JohnSnowLabs/country-and-continent-codes-list/r/0.html
# "0012833-190320150433242.csv" # GBIF Download DOI: 10.15468/dl.o73fzx

setwd(choose.dir())

#Create a country code to country/continent table
countries <- read.csv("country-and-continent-codes-list.csv")
country2continent <- countries[,c(1,4)]
colnames(country2continent)<-c("Continent","countryCode")

#Read data from GBIF data file (WARNING: change nThread parameter to an appropriate value)
bees.raw <- fread("0001373-190621201848488.csv", select=fields, nThread=6) %>% filter(family %in% all.fams, year!="")


####Count records, generate counts per species, year and family####
#Calculate a decade column for the records
bees.raw$decade <- floor(bees.raw$year/10)*10

#Remove records without ID at the species level
bees <- bees.raw %>% filter(species!="")

#From data, generate a table counting records per species and year
bee.species.yr <- count(bees,family,species, year)

#Add a species presence/absence field
bee.species.yr$r <- ifelse(bee.species.yr$n>0,1,0)

#Add total records per year by family
counts.by.yr <- as.data.frame(tapply(bee.species.yr$n, list(bee.species.yr$year, bee.species.yr$family), sum, na.rm=TRUE), keep.rownames = TRUE)

#Create a column holding the year as a numeric value, as they are wiped out from the row names by the mutate operation 
row.names(counts.by.yr) -> counts.by.yr$year
as.numeric(as.character(counts.by.yr$year))-> counts.by.yr$year

#Replace NA values for 0s
counts.by.yr <- counts.by.yr %>% mutate_if(is.integer, ~replace(., is.na(.), 0))

#Count total species recorded per year by family
species.by.yr <- as.data.frame(tapply(bee.species.yr$r, list(bee.species.yr$year, bee.species.yr$family), sum, na.rm=TRUE), keep.rownames = TRUE)

#Create a column holding the year as a numeric value, as they are wiped out from the row names by the mutate operation 
row.names(species.by.yr) -> species.by.yr$year
as.numeric(as.character(species.by.yr$year))-> species.by.yr$year

#Replace NA values for 0s
species.by.yr <- species.by.yr %>% mutate_if(is.numeric, ~replace(., is.na(.), 0))

#### Anthophila global trends####

anthophila.sp.yr <- bee.species.yr %>% filter(family %in% anthophila.fams) 

#Add total records per year
anto.counts.by.yr <- as.data.frame(tapply(anthophila.sp.yr$n, anthophila.sp.yr$year, sum, na.rm=TRUE), keep.rownames = TRUE)

#Create a column holding the year as a numeric value, as they are wiped out from the row names by the mutate operation 
row.names(anto.counts.by.yr) -> anto.counts.by.yr$year
as.numeric(as.character(anto.counts.by.yr$year))-> anto.counts.by.yr$year

#Replace NA values for 0s
anto.counts.by.yr <- anto.counts.by.yr %>% mutate_if(is.integer, ~replace(., is.na(.), 0))
names(anto.counts.by.yr)<-c("Records","year")

#Count total species recorded per year
anto.species.by.yr <- as.data.frame(tapply(anthophila.sp.yr$r, anthophila.sp.yr$year, sum, na.rm=TRUE), keep.rownames = TRUE)

#Create a column holding the year as a numeric value, as they are wiped out from the row names by the mutate operation 
row.names(anto.species.by.yr) -> anto.species.by.yr$year
as.numeric(as.character(anto.species.by.yr$year))-> anto.species.by.yr$year

#Replace NA values for 0s
anto.species.by.yr <- anto.species.by.yr %>% mutate_if(is.numeric, ~replace(., is.na(.), 0))
names(anto.species.by.yr)<-c("Richness","year")


####Figure 1A-ANTHOPHILA-Yearly counts of species and records####
plot(anto.species.by.yr$year, anto.species.by.yr$Richness, xlim=c(1920,2020), col="red", pch=20, cex=1, ylab = "Species richness", xlab = "Year")
lines(predict(loess(anto.species.by.yr$Richness~anto.species.by.yr$year, span=0.2))~anto.species.by.yr$year,col="red",  lwd=2)
par(new = TRUE)
plot(anto.counts.by.yr$year, anto.counts.by.yr$Records, xlim=c(1920,2018), pch=20, col=adjustcolor("red", alpha.f = 0.7) ,cex=0.5,xaxt = "n", yaxt = "n", ylab = "", xlab = "")
lines(predict(loess(anto.counts.by.yr$Records~anto.counts.by.yr$year, span=0.3))~anto.counts.by.yr$year, col=adjustcolor("red", alpha.f = 0.7))
axis(side = 4)

####Calculate Chao global stats for anthophila####
#Count occurrences per species and decade starting from 1950
ant.species.by.dec <- count(bees %>% filter(year>1949, family %in% anthophila.fams),decade,family,species)
anthophila.by.dec.mat <-as.data.frame(tapply(ant.species.by.dec$n, list(ant.species.by.dec$species,ant.species.by.dec$decade),sum, na.rm=T))
anthophila.by.dec.mat[is.na(anthophila.by.dec.mat)] <- 0
bees.global.chao.decades <- iNEXT(anthophila.by.dec.mat, q=0, datatype = "abundance")

####Figure 1B-ANTHOPHILA-iNEXT decadal curves####
plot(bees.global.chao.decades, show.main = F, col=viridis(length(levels(bees.global.chao.decades$AsyEst$Site)),direction=-1))

####Figure 1C-ANTHOPHILA-Asymptotic estimators per decade####
par(mfcol=c(1,1))
asyest.global.bees <- bees.global.chao.decades$AsyEst[seq(1,21,3),]
asyest.global.bees$decade <- as.numeric(levels(asyest.global.bees$Site))[asyest.global.bees$Site] #convert factor into numeric variable
plot(asyest.global.bees$decade, asyest.global.bees$Estimator, type="p", ylim=c(min(asyest.global.bees$LCL)*0.9,max(asyest.global.bees$UCL)*1.1), pch=16, ylab="Estimated richness", xlab="")
lines(asyest.global.bees$decade, asyest.global.bees$Estimator)
arrows(asyest.global.bees$decade, asyest.global.bees$LCL, asyest.global.bees$decade, asyest.global.bees$UCL, length=0.05, angle=90, code=3)

# To calculate the fraction of diversity, use the following line
asyest.global.bees[,4]/max(asyest.global.bees[,4])

####Family level analyses, including outgroups####
####Figure 2A - FAMILIES - Plotting records and species by year ####
# Plot number of species (primary, left axis) and number of records (secondary, right axis) by year (dots) for each family
par(mfcol=c(9,1))
par(mar=c(1,2,0,2))
for (i in 1:(length(all.fams)-1)) 
{
plot(species.by.yr$year, species.by.yr[,all.fams[i]], xlim=c(1920,2020), pch=20, col=all.famcols[i],cex=1, xaxt = "n", ylab = "", xlab = "")
lines(predict(loess(species.by.yr[,all.fams[i]]~species.by.yr$year, span=0.2))~species.by.yr$year, col=all.famcols[i], lwd=2)
par(new = TRUE)
plot(counts.by.yr$year, counts.by.yr[,all.fams[i]], xlim=c(1920,2020), pch=20, col=adjustcolor(all.famcols[i], alpha.f = 0.3) ,cex=0.5,xaxt = "n", yaxt = "n", ylab = "", xlab = "")
lines(predict(loess(counts.by.yr[,all.fams[i]]~counts.by.yr$year, span=0.2))~counts.by.yr$year, col=adjustcolor(all.famcols[i], alpha.f = 0.3))
axis(side = 4)
text(1930, 0.9*max(counts.by.yr[,all.fams[i]]), labels = paste(all.fams[i]))
}
par(mar=c(2,2,0,2))
i <- length(all.fams)
plot(species.by.yr$year, species.by.yr[,all.fams[i]], xlim=c(1920,2020), pch=20, col=all.famcols[i],cex=1, ylab = "", xlab = "")
lines(predict(loess(species.by.yr[,all.fams[i]]~species.by.yr$year, span=0.2))~species.by.yr$year, col=all.famcols[i], lwd=2)
par(new = TRUE)
plot(counts.by.yr$year, counts.by.yr[,all.fams[i]], xlim=c(1920,2020), pch=20, col=adjustcolor(all.famcols[i], alpha.f = 0.3),cex=0.5,xaxt = "n", yaxt = "n", ylab = "", xlab = "")
lines(predict(loess(counts.by.yr[,all.fams[i]]~counts.by.yr$year, span=0.2))~counts.by.yr$year, col=adjustcolor(all.famcols[i], alpha.f = 0.3))
axis(side = 4)
text(1930, 0.9*max(counts.by.yr[,all.fams[i]]), labels = paste(all.fams[i]))


####Calculate global Chao statistics for each family####
bee.species.fam <- count(bees,family,species)
global.chao <- list()
for (i in all.fams)
{
  famdata <- bee.species.fam %>% filter(family==i)
  famdata <- as.matrix(famdata[,c(-1,-2)])
  global.chao[[i]] <- iNEXT(famdata, q=0, datatype = "abundance")
}
#If need to extract Species richness indicators, use
global.chao[["Sphecidae"]]$AsyEst[1,]

#### Generate decadal counts and estimate Chao numbers with iNEXT ####
#Count occurrences per species and decade starting from 1950
species.by.dec <- count(bees %>% filter(year>1949),decade,family,species)
species.by.dec.mat <-as.data.frame(tapply(species.by.dec$n, list(species.by.dec$species,species.by.dec$decade),sum, na.rm=T))

#Add species names column and then merge this table with the species to family mapping table
bee.species <- count(bees,family,species)[,-3]
species.by.dec.mat$species<-row.names(species.by.dec.mat)
species.by.dec.mat <- merge(species.by.dec.mat, bee.species)

#Then create an empty list to hold the iNEXT output for each family and run iNEXT in each family
global.chao.decades <- list()
for (i in all.fams)
{
  famdata <- species.by.dec.mat %>% filter(family==i)
  famdata <- famdata[,c(-1,-9,-10)]
  famdata[is.na(famdata)] <- 0
  global.chao.decades[[i]] <- iNEXT(famdata, q=0, datatype = "abundance")
}
### WARNING! Running iNEXT for q values higher than 0 can result in much longer runtimes for the analysis.

#If need to extract Species richness indicators, use
global.chao.decades[["Apidae"]]$AsyEst[seq(1,21,3),]


####Figure 2B-All families-iNEXT decadal curves by family####
#Plot the interpolation-extrapolation curves by decade for each family
for (i in 1:(length(all.fams)-1)) 
{
plot(global.chao.decades[[i]], show.legend = F, show.main = F, col=viridis(length(levels(global.chao.decades[[i]]$AsyEst$Site)),direction=-1))
text(0,max(global.chao.decades[[i]]$AsyEst$Estimator)*0.95, labels=all.fams[i],adj=0)
}
plot(global.chao.decades[[9]], show.legend = T, show.main = F, col=viridis(length(levels(global.chao.decades[[9]]$AsyEst$Site)),direction=-1))
text(0,max(global.chao.decades[[9]]$AsyEst$Estimator)*0.95, labels=all.fams[9],adj=0)

####Figure 2C-All families-asymptotic estimators per decade and family####
#par(mfcol=c(1,1))
for (i in all.fams){
asyest <- global.chao.decades[[i]]$AsyEst[seq(1,21,3),]
asyest$decade <- as.numeric(levels(asyest$Site))[asyest$Site] #convert factor into numeric variable
plot(asyest$decade, asyest$Estimator, type="p", ylim=c(min(asyest$LCL)*0.9,max(asyest$UCL)*1.1), pch=16, ylab="Estimated richness", xlab="", main=i)
lines(asyest$decade, asyest$Estimator)
arrows(asyest$decade, asyest$LCL, asyest$decade, asyest$UCL, length=0.05, angle=90, code=3)
}

## Fraction of diversity per family
perc.decline<-as.numeric()
j<-1
for(i in all.fams)
  { perc.decline[j]<-global.chao.decades[[i]]$AsyEst[19,4]/max(global.chao.decades[[i]]$AsyEst[seq(1,21,3),4])
  j<-j+1}
names(perc.decline)<-all.fams

1-max(perc.decline[4:9])
1-min(perc.decline[4:9])


#### Figure 2 COMPLETE - All families, year and decadal plots by family ####
par(mfcol=c(9,3))
par(mar=c(1,2,1,2))
for (i in 1:(length(all.fams)-1)) 
{
  plot(species.by.yr$year, species.by.yr[,all.fams[i]], xlim=c(1920,2020), pch=20, col=all.famcols[i],cex=1, xaxt = "n", ylab = "", xlab = "")
  lines(predict(loess(species.by.yr[,all.fams[i]]~species.by.yr$year, span=0.2))~species.by.yr$year, col=all.famcols[i], lwd=2)
  par(new = TRUE)
  plot(counts.by.yr$year, counts.by.yr[,all.fams[i]], xlim=c(1920,2020), pch=20, col=adjustcolor(all.famcols[i], alpha.f = 0.3) ,cex=0.5,xaxt = "n", yaxt = "n", ylab = "", xlab = "")
  lines(predict(loess(counts.by.yr[,all.fams[i]]~counts.by.yr$year, span=0.2))~counts.by.yr$year, col=adjustcolor(all.famcols[i], alpha.f = 0.3))
  axis(side = 4)
}
par(mar=c(2,2,1,2))
i <- length(all.fams)
plot(species.by.yr$year, species.by.yr[,all.fams[i]], xlim=c(1920,2020), pch=20, col=all.famcols[i],cex=1, ylab = "", xlab = "")
lines(predict(loess(species.by.yr[,all.fams[i]]~species.by.yr$year, span=0.2))~species.by.yr$year, col=all.famcols[i], lwd=2)
par(new = TRUE)
plot(counts.by.yr$year, counts.by.yr[,all.fams[i]], xlim=c(1920,2020), pch=20, col=adjustcolor(all.famcols[i], alpha.f = 0.3),cex=0.5,xaxt = "n", yaxt = "n", ylab = "", xlab = "")
lines(predict(loess(counts.by.yr[,all.fams[i]]~counts.by.yr$year, span=0.2))~counts.by.yr$year, col=adjustcolor(all.famcols[i], alpha.f = 0.3))
axis(side = 4)

par(mar=c(1,2,1,2))
for (i in 1:(length(all.fams)-1)) 
{
  plot(global.chao.decades[[i]], show.legend = F, show.main = F, col=viridis(length(levels(global.chao.decades[[i]]$AsyEst$Site)),direction=-1), yaxt = "n")
  axis(side = 4)
}
par(mar=c(2,2,1,2))
i <- length(all.fams)
plot(global.chao.decades[[i]], show.legend = F, show.main = F, col=viridis(length(levels(global.chao.decades[[i]]$AsyEst$Site)),direction=-1), yaxt = "n")
axis(side = 4)

for (i in all.fams){
  asyest <- global.chao.decades[[i]]$AsyEst[seq(1,21,3),]
  asyest$decade <- as.numeric(levels(asyest$Site))[asyest$Site] #convert factor into numeric variable
  plot(asyest$decade, asyest$Estimator, type="p", ylim=c(min(asyest$LCL)*0.9,max(asyest$UCL)*1.1), pch=16, ylab="Estimated richness", xlab="")
  lines(asyest$decade, asyest$Estimator)
  arrows(asyest$decade, asyest$LCL, asyest$decade, asyest$UCL, length=0.05, angle=90, code=3)
  mtext(i, side=4)
}

####Figure 3 - Pielou's eveness index over time####
#ESTIMATE GLOBAL EQUITATIVITY USING PIELOU'S J EVENNESS
anto.sp.yr.table <- xtabs(n~year+species, data=anthophila.sp.yr %>% filter(year>1899, year<2019))
H <- diversity(anto.sp.yr.table)
J <- H/log(specnumber(anto.sp.yr.table))

par(mfrow=c(1,1))
plot(as.numeric(names(J)),J, pch=16, ylab="Pielou's index",xlab="Year")
lines(predict(loess(J~as.numeric(names(J)), span=0.3))~as.numeric(names(J)),col="red",  lwd=2)



####Re-run analyses including HUMAN_OBSERVATION records####
#Original data: https://doi.org/10.15468/dl.o73fzx
#License: CC BY-NC 4.0 
#File: 535 MB CSV (Zipped)
#Involved datasets: 1977
#fields<- c("gbifID","speciesKey","family","genus","species","countryCode","decimalLatitude","decimalLongitude","year","institutionCode","datasetKey", "basisOfRecord")
morebees <- fread("0012833-190320150433242.csv", select=fields, nThread=6) %>% filter(family %in% anthophila.fams, species!="", year!="")
#Calculate a decade column for the records
morebees$decade <- floor(morebees$year/10)*10

####Figure S1####
basis.per.decade <- table(morebees$basisOfRecord,morebees$decade)
par(mfrow=c(1,1))
par(mar=c(4,5,1,1))
barplot(basis.per.decade[c(2,1),27:33]/1000, xlab="Decade", col=c("darkorange","red"), 
        legend = rownames(basis.per.decade[c(2,1),27:33]), args.legend = list(x="topleft", bty="n"),
        ylab="Records (thousands)")

#From data, generate a table counting records per species and year
anthophila.sp.yr <- count(morebees,family,species, year)
#Add a species presence/absence field
anthophila.sp.yr$r <- ifelse(anthophila.sp.yr$n>0,1,0)
anthophila.sp.yr <- anthophila.sp.yr %>% filter(year<2019)

#Add total records per year
anto.counts.by.yr <- as.data.frame(tapply(anthophila.sp.yr$n, anthophila.sp.yr$year, sum, na.rm=TRUE), keep.rownames = TRUE)
#Create a column holding the year as a numeric value, as they are wiped out from the row names by the mutate operation 
row.names(anto.counts.by.yr) -> anto.counts.by.yr$year
as.numeric(as.character(anto.counts.by.yr$year))-> anto.counts.by.yr$year
#Replace NA values for 0s
anto.counts.by.yr <- anto.counts.by.yr %>% mutate_if(is.integer, ~replace(., is.na(.), 0))
names(anto.counts.by.yr)<-c("Records","year")

#Count total species recorded per year
anto.species.by.yr <- as.data.frame(tapply(anthophila.sp.yr$r, anthophila.sp.yr$year, sum, na.rm=TRUE), keep.rownames = TRUE)
#Create a column holding the year as a numeric value, as they are wiped out from the row names by the mutate operation 
row.names(anto.species.by.yr) -> anto.species.by.yr$year
as.numeric(as.character(anto.species.by.yr$year))-> anto.species.by.yr$year
#Replace NA values for 0s
anto.species.by.yr <- anto.species.by.yr %>% mutate_if(is.numeric, ~replace(., is.na(.), 0))
names(anto.species.by.yr)<-c("Richness","year")

#Count occurrences per species and decade starting from 1950
species.by.dec <- count(morebees %>% filter(year>1949),decade,family,species)
species.by.dec.mat <-as.data.frame(tapply(species.by.dec$n, list(species.by.dec$species,species.by.dec$decade),sum, na.rm=T))

#Add species names column and then merge this table with the species to family mapping table
bee.species <- count(morebees,family,species)[,-3]
species.by.dec.mat$species<-row.names(species.by.dec.mat)
species.by.dec.mat <- merge(species.by.dec.mat, bee.species)

#Now calculate Chao global stats for anthophila  
anthophila.by.dec.mat <- species.by.dec.mat[,c(-1,-9,-10)]
anthophila.by.dec.mat[is.na(anthophila.by.dec.mat)] <- 0
allbees.global.chao.decades <- iNEXT(anthophila.by.dec.mat, q=0, datatype = "abundance")


####Figure S2####
par(mfrow=c(1,3))
#Plot yearly counts of species and records
plot(anto.species.by.yr$year, anto.species.by.yr$Richness, xlim=c(1920,2020), col="red", pch=20, cex=1, ylab = "Species Richness", xlab = "Year")
lines(predict(loess(anto.species.by.yr$Richness~anto.species.by.yr$year, span=0.2))~anto.species.by.yr$year,col="red",  lwd=2)
par(new = TRUE)
plot(anto.counts.by.yr$year, anto.counts.by.yr$Records, xlim=c(1920,2018), pch=20, col="orange" ,cex=0.5,xaxt = "n", yaxt = "n", ylab = "", xlab = "")
lines(predict(loess(anto.counts.by.yr$Records~anto.counts.by.yr$year, span=0.3))~anto.counts.by.yr$year, col="orange")
axis(side = 4)

#Plot iNEXT curves
plot(allbees.global.chao.decades, show.legend = T, show.main = F, col=viridis(length(levels(allbees.global.chao.decades$AsyEst$Site)),direction=-1))

#Plot asymptotic estimators per decade
asyest <- allbees.global.chao.decades$AsyEst[seq(1,21,3),]
asyest$decade <- as.numeric(levels(asyest$Site))[asyest$Site] #convert factor into numeric variable
plot(asyest$decade, asyest$Estimator, type="p", ylim=c(min(asyest$LCL)*0.9,max(asyest$UCL)*1.1), pch=16, ylab="Estimated richness", xlab="")
lines(asyest$decade, asyest$Estimator)
arrows(asyest$decade, asyest$LCL, asyest$decade, asyest$UCL, length=0.05, angle=90, code=3)

#ESTIMATE GLOBAL EQUITATIVITY USING PIELOU'S J EVENNESS
anto.sp.yr.table <- xtabs(n~year+species, data=anthophila.sp.yr %>% filter(year>1899, year<2019))
H <- diversity(anto.sp.yr.table)
J <- H/log(specnumber(anto.sp.yr.table))
par(mfrow=c(1,1))
plot(as.numeric(names(J)),J, pch=16, ylab="Pielou's evenness",xlab="")
lines(predict(loess(J~as.numeric(names(J)), span=0.3))~as.numeric(names(J)),col="red",  lwd=2)

####Figure S3 - ITIS species counts vs asymptotic estimators####

#Species counts obtained from the Integrated Taxonomic Information System (ITIS, www.itis.gov).
itis.sp.cts <-c(10213,784,8871,206,3010,4186,2493,5062,5677)
names(itis.sp.cts)<-all.fams

chao.sp.cts <-as.numeric()
for(i in all.fams) {chao.sp.cts[i]<- global.chao[[i]]$AsyEst[1,4]}
chao.sp.cts <-unlist(chao.sp.cts)

par(mfcol=c(1,1))
plot(itis.sp.cts/1000, chao.sp.cts/1000, pch=16, col=all.famcols,
     xlab = "ITIS Species Count (thousands)",
     ylab = "Asymptotic Estimator of Richness (thousands)")
text(itis.sp.cts/1000, chao.sp.cts/1000, label=names(itis.sp.cts))
segments(0,0,10,10, lty=3)
legend("topleft",legend=names(itis.sp.cts), col=all.famcols, pch=16, bty="n")

#Fraction of known species estimated using Chao for all families in the analyses
sum(chao.sp.cts)/sum(itis.sp.cts)

#Fraction of known species estimated using Chao for anthophila families in the analyses
sum(chao.sp.cts[4:9])/sum(itis.sp.cts[4:9])

plot(chao.sp.cts/itis.sp.cts, pch=16, col=all.famcols)


####Figure S4 - Counting fraction of species with no ID####
#Filter anthophila families from unfiltered original GBIF data
antho.raw <- bees.raw %>% filter(family %in% anthophila.fams)
#Label and count yearly number of records with no ID to species
antho.raw$sp.id <- ifelse(antho.raw$species!="",1,0)
specific.n.yr <- count(antho.raw, sp.id, year)
spid.yr<- xtabs(n~sp.id+year, data=specific.n.yr%>%filter(year>1899, year<2019))
#Estimate and plot yearly fraction of records missing species ID
spid.yr.frq <- colPerc(spid.yr)/100
freq.misid <- data.frame(as.numeric(names(spid.yr.frq[1,])),spid.yr.frq[1,])
names(freq.misid)<-c("year","freq.misid")
plot(freq.misid, pch=16, xlab="Year",ylab="Fraction of records missing species ID")
lines(predict(loess(freq.misid$freq.misid~freq.misid$year, span=0.3))~freq.misid$year,col="red",  lwd=2)


#### Generate trends for Anthophila per continent####

#Count occurrences per decade and country
bees.cont <- merge(bees %>% filter(countryCode!="", year>1899, year<2019, family %in% anthophila.fams),country2continent) 

#From data, generate a table counting records per species and year
bee.species.cont.yr <- count(bees.cont,Continent,species, year)

#Add a species presence/absence field
bee.species.cont.yr$r <- ifelse(bee.species.cont.yr$n>0,1,0)

#Add total records per year by family
cont.counts.by.yr <- as.data.frame(tapply(bee.species.cont.yr$n, list(bee.species.cont.yr$year, bee.species.cont.yr$Continent), sum, na.rm=TRUE), keep.rownames = TRUE)
#Create a column holding the year as a numeric value, as they are wiped out from the row names by the mutate operation 
row.names(cont.counts.by.yr) -> cont.counts.by.yr$year
as.numeric(as.character(cont.counts.by.yr$year))-> cont.counts.by.yr$year
#Replace NA values for 0s
cont.counts.by.yr <- cont.counts.by.yr %>% mutate_if(is.integer, ~replace(., is.na(.), 0))
#Remove Antarctica
cont.counts.by.yr<-cont.counts.by.yr[,-2]

#Count total species recorded per year by family
cont.species.by.yr <- as.data.frame(tapply(bee.species.cont.yr$r, list(bee.species.cont.yr$year, bee.species.cont.yr$Continent), sum, na.rm=TRUE), keep.rownames = TRUE)
#Create a column holding the year as a numeric value, as they are wiped out from the row names by the mutate operation 
row.names(cont.species.by.yr) -> cont.species.by.yr$year
as.numeric(as.character(cont.species.by.yr$year))-> cont.species.by.yr$year
#Replace NA values for 0s
cont.species.by.yr <- cont.species.by.yr %>% mutate_if(is.numeric, ~replace(., is.na(.), 0))
#Remove Antarctica
cont.species.by.yr<-cont.species.by.yr[,-2]

#### Generate decadal counts and estimate Chao numbers with iNEXT ####

#Count occurrences per species and decade starting from 1950
species.by.dec.cont <- count(bees.cont %>% filter(year>1949),decade,Continent,species)

#Create an empty list to hold the iNEXT output for each continent and run iNEXT 
continental.anto.chao.decades <- list()
for (i in unique(species.by.dec.cont$Continent))
{
  contsp <- species.by.dec.cont %>% filter(Continent==i)
  contsp.mat <-as.data.frame(tapply(contsp$n, list(contsp$species,contsp$decade), sum, na.rm=T))
  contsp.mat[is.na(contsp.mat)] <- 0
  continental.anto.chao.decades[[i]] <- iNEXT(contsp.mat, q=0, datatype = "abundance")
}

####Figure S5 - records by decade and continent####
continent.per.decade <- table(bees.cont$Continent,bees.cont$decade)

#Remove Antarctica and decades before the 1950's
continent.per.decade <-continent.per.decade[-2,6:12]
continent.per.decade <-continent.per.decade[6:1,]
par(mfrow=c(1,1))
par(mar=c(4,5,1,1))
par(mfcol=c(1,2))
barplot(continent.per.decade/1000, xlab="Decade", legend = rownames(continent.per.decade), args.legend = list(x="topleft", bty="n"),ylab="Records (thousands)", col = brewer.pal(n = 6, name = "RdBu"))
barplot(colPerc(continent.per.decade)[-7,], xlab="Decade", ylab="% of records", col = brewer.pal(n = 6, name = "RdBu"))

####Figure S6 Alt2 - Continental patterns ####
par(mfcol=c(6,3))
par(mar=c(1,2,1,2))
for (i in 1:6) 
{
  plot(cont.species.by.yr$year, cont.species.by.yr[,i], xlim=c(1920,2020), pch=20, cex=1, ylab = "", xlab = "")
  lines(predict(loess(cont.species.by.yr[,i]~cont.species.by.yr$year, span=0.2))~cont.species.by.yr$year, lwd=2)
  par(new = TRUE)
  plot(cont.counts.by.yr$year, cont.counts.by.yr[,i], xlim=c(1920,2020), pch=20, col="gray",cex=0.5,xaxt = "n", yaxt = "n", ylab = "", xlab = "")
  lines(predict(loess(cont.counts.by.yr[,i]~cont.counts.by.yr$year, span=0.2))~cont.counts.by.yr$year, col="gray")
  axis(side = 4)
}

#Plot iNEXT curves by continent
for (i in unique(species.by.dec.cont$Continent))
{ plot(continental.anto.chao.decades[[i]], show.legend = F, show.main = F, col=viridis(length(levels(continental.anto.chao.decades[[i]]$AsyEst$Site)),direction=-1))
}

#Plot asymptotic estimators per decade and continent
for (i in unique(species.by.dec.cont$Continent)){
  asyest <- continental.anto.chao.decades[[i]]$AsyEst[seq(1,21,3),]
  asyest$decade <- as.numeric(levels(asyest$Site))[asyest$Site] #convert factor into numeric variable
  plot(asyest$decade, asyest$Estimator, type="p", ylim=c(min(asyest$LCL)*0.9,max(asyest$UCL)*1.1), pch=16, xlab="", ylab="")
  lines(asyest$decade, asyest$Estimator)
  arrows(asyest$decade, asyest$LCL, asyest$decade, asyest$UCL, length=0.05, angle=90, code=3)
}

## Fraction of diversity per continent
perc.decline.cont<-as.numeric()
j<-1
for(i in unique(species.by.dec.cont$Continent))
{ perc.decline.cont[j]<-continental.anto.chao.decades[[i]]$AsyEst[19,4]/max(continental.anto.chao.decades[[i]]$AsyEst[seq(1,21,3),4])
j<-j+1}
names(perc.decline.cont)<-unique(species.by.dec.cont$Continent)

1-max(perc.decline[4:9])
1-min(perc.decline[4:9])

#####Estimate annual fraction of records due to a given species######
#From data, generate a table counting records per species and year
bee.species.yr <- count(bees,family,species, year)
bee.species.yr.table <- tapply(bee.species.yr$n, list(bee.species.yr$species, bee.species.yr$year), sum, na.rm=TRUE)
bee.species.yr.table[is.na(bee.species.yr.table)] <- 0
bee.species.yr.table.Perc <- colPerc(bee.species.yr.table)
data.years <- as.numeric(names(bee.species.yr.table.Perc[1,]))

####Figure S7 - Fraction of records due to Apis mellifera####
apis.exp.fit <- lm(log(bee.species.yr.table.Perc["Apis mellifera",114:233]+0.00001)~data.years[114:233])
apis.exp.pred <- exp(predict(apis.exp.fit,list(Time=data.years[114:233])))
par(mfrow=c(1,1), mar=c(4,4,2,2))
plot(data.years[114:233],bee.species.yr.table.Perc["Apis mellifera",114:233], pch=16, ylab="% records due to Apis mellifera", xlab="Year", ylim=c(0,8))
lines(apis.exp.pred~data.years[114:233], col="red", lwd=2)

####END OF SCRIPT####
#PARTIAL OR TOTAL USE THIS CODE OR A DERIVATIVE FOR A DIFFERENT ANALYSIS IS ALLOWED, 
#BUT PLEASE DO NOT FORGET TO GIVE PROPER ATTRIBUTION AND TO CITE THE COMPANION PAPER TO THIS CODE
#THANKS!
#####