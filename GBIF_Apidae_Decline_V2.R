####HEADER####
# Global decline of species richness in records of bee occurrence
#
# An analysis of multidecadal trends in bee species richness from GBIF data
#
# Eduardo E. Zattara & Marcelo A. Aizen
#
# Grupo de Ecología de la Polinización, INIBIOMA, 
# Universidad Nacional del Comahue-CONICET, 
# Quintral 1250, Bariloche (8400), Argentina.
# Correspondence to ezattara@comahue-conicet.gob.ar

#Data source: https://doi.org/10.15468/dl.ysjm4x
#License: CC BY-NC 4.0 
#File: 715 MB CSV (Zipped)
#Involved datasets: 2374
#Filters: 
#    Scientific Name == "Hymenoptera" 
#    Basis of Record == c("PRESERVED_SPECIMEN", "HUMAN_OBSERVATION")
#    Year < 2020

#### Setting up data ---------------------------

#Load necessary libraries
library(iNEXT)
library("data.table")
library(dplyr)
library(tidyr)
library(ggplot2)
library(egg)
library(viridisLite)
library(vegan)
library(nlme)

#Select working directory where the following data files are found:
# "0001373-190621201848488.csv" # GBIF Download from https://doi.org/10.15468/dl.ysjm4x
# "country-and-continent-codes-list.csv" # Country and continent codes, from https://datahub.io/JohnSnowLabs/country-and-continent-codes-list/r/0.html

#Enter here current GBIF download to analyze -- Change if using different file
GBIFdata <- "0058082-200221144449610.csv" # All Hymenoptera PRESERVED_SPECIMEN + HUMAN_OBSERVATION until 2020

#Create a country code to country/continent table
countries <- read.csv("country-and-continent-codes-list.csv")
country2continent <- countries[,c(1,4)]
colnames(country2continent)<-c("Continent","countryCode")

#Select fields for GBIF data
fields<- c("gbifID","speciesKey","family","genus","species","countryCode","decimalLatitude","decimalLongitude","year","institutionCode","datasetKey","taxonRank","basisOfRecord","collectionCode")

#Select Hymenopteran families for the analysis
anthophila_fams <-c("Melittidae","Andrenidae","Halictidae","Colletidae","Megachilidae","Apidae")
outgroup_fams <- c("Formicidae","Sphecidae","Crabronidae")
all_fams <- c(outgroup_fams,anthophila_fams)

#Read data from GBIF data file (WARNING: change nThread parameter to an appropriate value) for Anthophila, filtering out records without a species id
bees_all <-fread(GBIFdata, select=fields, nThread=6) %>% filter(family %in% anthophila_fams, taxonRank=="SPECIES", year!="", species!="")

#Calculate the "idecade" columns for the records
bees_all$idecade <- floor((bees_all$year + 4)/10)*10


#### Count numbers records and species reported per year ----------------
#Count number of records by species and year (for full and specimens-only datasets)
bees_all_recs_by_sp_year <- bees_all %>% group_by(species, family, year, idecade) %>% count()
bees_prsv_recs_by_sp_year <- bees_all %>% filter(basisOfRecord=="PRESERVED_SPECIMEN")%>% group_by(species, family, year, idecade) %>% count()

#Histograms of abundance of records
ggplot(bees_all_recs_by_sp_year %>% filter(idecade>1940), aes(log(n))) +
  geom_histogram() +
  facet_grid(.~idecade) + scale_y_log10()

ggplot(bees_prsv_recs_by_sp_year %>% filter(idecade>1940), aes(log(n))) +
  geom_histogram() +
  facet_grid(.~idecade) + scale_y_log10()

#Count number of species per year (for complete and preserved specimens only datasets)
bees_all_sp_year <- bees_all_recs_by_sp_year %>% select(-n)  %>% group_by(year) %>% count() %>% rename(sp = n)
bees_prsv_sp_year <- bees_prsv_recs_by_sp_year %>% select(-n) %>% group_by(year) %>% count() %>% rename(sp = n)

#Count number of records per year (for complete and preserved specimens only datasets)
bees_all_recs_year <- bees_all %>% group_by(year) %>% count() %>% rename(records = n)
bees_prsv_recs_year <- bees_all %>% filter(basisOfRecord=="PRESERVED_SPECIMEN")%>% group_by(year) %>% count()%>% rename(records = n)

#Count number of datasets per year (for complete and preserved specimens only datasets)
bees_all_recs_by_dataset_year <- bees_all %>% group_by(year, datasetKey) %>% count() 
bees_all_datasets_year <- bees_all_recs_by_dataset_year %>% select(-n) %>% group_by(year) %>% count() %>% rename(datasets = n)

bees_prsv_recs_by_dataset_year <- bees_all %>% filter(basisOfRecord=="PRESERVED_SPECIMEN")%>% group_by(year, datasetKey) %>% count()
bees_prsv_datasets_year <- bees_prsv_recs_by_dataset_year %>% select(-n) %>% group_by(year) %>% count() %>% rename(datasets = n)

#Count number of collections per year (for complete and preserved specimens only datasets)
bees_all_recs_by_collection_year <- bees_all %>% group_by(year, collectionCode) %>% count()
bees_all_collections_year <- bees_all_recs_by_collection_year %>% select(-n)  %>% group_by(year) %>% count() %>% rename(collections = n)

bees_prsv_recs_by_collection_year <- bees_all %>% filter(basisOfRecord=="PRESERVED_SPECIMEN")%>% group_by(year, collectionCode) %>% count()
bees_prsv_collection_year <- bees_prsv_recs_by_collection_year %>% select(-n) %>% group_by(year) %>% count() %>% rename(collections = n)

#Count number of institutions per year (for complete and preserved specimens only datasets)
bees_all_recs_by_institution_year <- bees_all %>% group_by(year, institutionCode) %>% count()
bees_all_institution_year <- bees_all_recs_by_institution_year %>% select(-n) %>% group_by(year) %>% count() %>% rename(institutions = n)

bees_prsv_recs_by_institution_year <- bees_all %>% filter(basisOfRecord=="PRESERVED_SPECIMEN")%>% group_by(year, institutionCode) %>% count()
bees_prsv_institution_year <- bees_prsv_recs_by_institution_year %>% select(-n) %>% group_by(year) %>% count() %>% rename(institutions = n)

#Merge all yearly counts
bees_all_gbif_year <- full_join(bees_all_sp_year, bees_all_recs_year, by = "year") %>% full_join(bees_all_datasets_year, by = "year") %>% full_join(bees_all_collections_year, by = "year") %>% full_join(bees_all_institution_year, by = "year")
bees_prsv_gbif_year <- full_join(bees_prsv_sp_year, bees_prsv_recs_year, by = "year") %>% full_join(bees_prsv_datasets_year, by = "year") %>% full_join(bees_prsv_collection_year, by = "year") %>% full_join(bees_prsv_institution_year, by = "year")

#### Testing for negative trend on species counts after 1985 ----

# Linear model on full dataset
bees_all_lm <- lm(sp ~  year+ records + collections + institutions + datasets, data = bees_all_gbif_year %>% filter(year>1985, year<2016))
plot(bees_all_lm)
summary(bees_all_lm)

# Generalized least squares on full dataset
bees_all_gls <- gls(sp ~ year + records + collections + institutions + datasets, data = bees_all_gbif_year %>% filter(year>1985, year<2016))
summary(bees_all_gls)

# Generalized least squares with an autocorrelation-moving average correlation structure of order 1,0 on full dataset
bees_all_glsau <- gls(data=bees_all_gbif_year %>% filter(year>1985, year<2016), 
                      sp ~  year+ records + collections + institutions + datasets,
                      correlation = corARMA(form = ~year, p=1, q=0))
summary(bees_all_glsau)

#Comparing gls models with and without autocorrelation
anova(bees_all_gls, bees_all_glsau)

# Linear model on specimens-only dataset
bees_prsv_lm <- lm(sp ~  year+ records + collections + institutions + datasets, data = bees_prsv_gbif_year %>% filter(year>1985, year<2016))
plot(bees_prsv_lm)
summary(bees_prsv_lm)

# Generalized least squares on specimens-only dataset
bees_prsv_gls <- gls(sp ~  year + records + collections + institutions + datasets, data = bees_prsv_gbif_year %>% filter(year>1985, year<2016))
summary(bees_prsv_gls)

# Generalized least squares with an autocorrelation-moving average correlation structure of order 1,0 on specimens-only dataset
bees_prsv_glsau <- gls(data=bees_prsv_gbif_year %>% filter(year>1985, year<2016), 
                      sp ~  year+ records + collections + institutions + datasets,
                      correlation = corARMA(form = ~year, p=1, q=0))
summary(bees_prsv_glsau)

#Comparing gls models with and without autocorrelation
anova(bees_prsv_gls, bees_prsv_glsau)


####Plot raw number of species and records per year for specimens-only and full datasets--------------------------
#Plot number of records by year
plot_recs <-
ggplot() +
  geom_point(data = bees_all_recs_year %>% filter(year<2020, year>1900), aes(x = year, y = records), color = "#00AFBB") +
  geom_point(data = bees_all_recs_year %>% filter(year<2020, year>2015), aes(x = year, y = records+500), color = "black", shape = 8) +
  geom_smooth(data = bees_all_recs_year %>% filter(year<2016, year>1900), aes(x = year, y = records), color = "#00AFBB", fill = "lightblue") +
  geom_point(data = bees_prsv_recs_year %>% filter(year<2020, year>1900), aes(x = year, y = records), color = "red") +
  geom_point(data = bees_prsv_recs_year %>% filter(year<2020, year>2015), aes(x = year, y = records+500), color = "black", shape = 8) +
  geom_smooth(data = bees_prsv_recs_year %>% filter(year<2016, year>1900), aes(x = year, y = records), color = "red", fill = "pink") +
  ylab("No. records") + xlab("Year")


#Plot number of species by year
plot_sp <-
ggplot() +
  geom_point(data = bees_all_sp_year %>% filter(year<2020, year>1900), aes(x = year, y = sp), color = "#00AFBB") +
  geom_point(data = bees_all_sp_year %>% filter(year<2020, year>2015), aes(x = year, y = sp), color = "blue", shape = 19) +
  geom_point(data = bees_all_sp_year %>% filter(year<2020, year>2015), aes(x = year, y = sp+20), color = "black", shape = 8) +
  geom_smooth(data = bees_all_sp_year %>% filter(year<2016, year>1900), aes(x = year, y = sp), color = "#00AFBB", fill = "lightblue") +
  geom_point(data = bees_prsv_sp_year %>% filter(year<2020, year>1900), aes(x = year, y = sp), color = "red") +
  geom_point(data = bees_prsv_sp_year %>% filter(year<2020, year>2015), aes(x = year, y = sp), color = "darkred", shape = 19) +
  geom_point(data = bees_prsv_sp_year %>% filter(year<2020, year>2015), aes(x = year, y = sp+20), color = "black", shape = 8) +
  geom_smooth(data = bees_prsv_sp_year %>% filter(year<2016, year>1900), aes(x = year, y = sp), color = "red", fill = "pink") +
  ylab("species richness") + xlab("Year")


#Plot number of datasets by year
ggplot() +
  geom_point(data = bees_all_datasets_year %>% filter(year<2020, year>1920), aes(x = year, y = datasets), color = "#00AFBB") +
  geom_smooth(data = bees_all_datasets_year %>% filter(year<2016, year>1920), aes(x = year, y = datasets), color = "#00AFBB", fill = "lightblue") +
  geom_point(data = bees_prsv_datasets_year %>% filter(year<2020, year>1920), aes(x = year, y = datasets), color = "red") +
  geom_smooth(data = bees_prsv_datasets_year %>% filter(year<2016, year>1920), aes(x = year, y = datasets), color = "red", fill = "pink") +
  ylim(0, max(bees_all_datasets_year$datasets)*1.05)

#Plot number of collections by year
ggplot() +
  geom_point(data = bees_all_collections_year %>% filter(year<2020, year>1920), aes(x = year, y = collections), color = "#00AFBB") +
  geom_smooth(data = bees_all_collections_year %>% filter(year<2016, year>1920), aes(x = year, y = collections), color = "#00AFBB", fill = "lightblue") +
  geom_point(data = bees_prsv_collection_year %>% filter(year<2020, year>1920), aes(x = year, y = collections), color = "red") +
  geom_smooth(data = bees_prsv_collection_year %>% filter(year<2016, year>1920), aes(x = year, y = collections), color = "red", fill = "pink") 

#Plot number of institutions by year
ggplot() +
  geom_point(data = bees_all_institution_year %>% filter(year<2020, year>1920), aes(x = year, y = institutions), color = "#00AFBB") +
  geom_smooth(data = bees_all_institution_year %>% filter(year<2016, year>1920), aes(x = year, y = institutions), color = "#00AFBB", fill = "lightblue") +
  geom_point(data = bees_prsv_institution_year %>% filter(year<2020, year>1920), aes(x = year, y = institutions), color = "red") +
  geom_smooth(data = bees_prsv_institution_year %>% filter(year<2016, year>1920), aes(x = year, y = institutions), color = "red", fill = "pink") 

#### Use iNEXT to generate rarefaction curves and estimate species diversity for each idecade ----------------

#Count number of records by species and idecade (for complete and preserved specimens only datasets)
bees_all_recs_by_sp_idec <- bees_all %>% group_by(species, family, idecade) %>% count()
bees_prsv_recs_by_sp_idec <- bees_all %>% filter(basisOfRecord=="PRESERVED_SPECIMEN")%>% group_by(species, family, idecade) %>% count()

#Generate idecadal matrices and calculate Chao stats
bees_all_recs_by_sp_idec.mat <-with(bees_all_recs_by_sp_idec %>% filter(idecade>1940, idecade<2020), as.data.frame(tapply(n, list(species,idecade),sum, na.rm=T)))
bees_all_recs_by_sp_idec.mat[is.na(bees_all_recs_by_sp_idec.mat)] <- 0
bees_all_global_chao_idec <- iNEXT(bees_all_recs_by_sp_idec.mat, q=0, datatype = "abundance")

bees_prsv_recs_by_sp_idec.mat <-with(bees_prsv_recs_by_sp_idec %>% filter(idecade>1940, idecade<2020), as.data.frame(tapply(n, list(species,idecade),sum, na.rm=T)))
bees_prsv_recs_by_sp_idec.mat[is.na(bees_prsv_recs_by_sp_idec.mat)] <- 0
bees_prsv_global_chao_idec <- iNEXT(bees_prsv_recs_by_sp_idec.mat, q=0, datatype = "abundance")


#### Plot iNEXT curves ------------------
plot_chao_prsv <- 
  ggiNEXT(bees_prsv_global_chao_idec, type=1) +
  scale_shape_manual(values=c(19,19,19,19,19,19,19)) + xlim(0,1200000) +
  scale_colour_manual(values=viridis(length(levels(bees_prsv_global_chao_idec$AsyEst$Site)),direction=-1)) +
  theme(legend.position="none", axis.text=element_text(size=10),axis.title=element_text(size=12) ) + 
  xlab("No. records (full dataset)") + ylab("Species richness")


plot_chao_all <- 
  ggiNEXT(bees_all_global_chao_idec, type=1, ) +
  scale_shape_manual(values=c(19,19,19,19,19,19,19))  + xlim(0,1200000) +
  scale_colour_manual(values=viridis(length(levels(bees_all_global_chao_idec$AsyEst$Site)),direction=-1)) +
  theme(legend.position="none", axis.text=element_text(size=10),axis.title=element_text(size=12)) + 
  xlab("No. records (specimens-only dataset)") + ylab("Species richness")

grid.arrange(plot_chao_all, plot_chao_prsv, ncol = 1)


#### Extract and plot estimators as barplots ----

#Specimens-only dataset
chaoasyest_global_bees_prsv_idec <- bees_prsv_global_chao_idec$AsyEst[seq(1,21,3),]

plot_est_prsv <- 
ggplot(chaoasyest_global_bees_prsv_idec, aes(Site, Estimator)) +
  geom_bar(stat="Identity") +
  coord_cartesian(ylim = c(4000, 7000)) + 
  geom_errorbar(aes(ymin=LCL, ymax=UCL), width=.2,
                position=position_dodge(.9)) +
  xlab("iDecade")  + ylab("Species richness")


#Calculate loss in specimen dataset during last two idecades relative to the average of previous idecades
1-chaoasyest_global_bees_prsv_idec[6,4]/mean(chaoasyest_global_bees_prsv_idec[1:5,4]) # Percent loss during the 2000s
1-chaoasyest_global_bees_prsv_idec[7,4]/mean(chaoasyest_global_bees_prsv_idec[1:5,4]) # Percent loss during the 2010s

#Full dataset 
chaoasyest_global_bees_all_idec <- bees_all_global_chao_idec$AsyEst[seq(1,21,3),]

plot_est_all <- 
  ggplot(chaoasyest_global_bees_all_idec, aes(Site, Estimator)) +
  geom_bar(stat="Identity") +
  coord_cartesian(ylim = c(4000, 7000)) + 
  geom_errorbar(aes(ymin=LCL, ymax=UCL), width=.2,
                position=position_dodge(.9)) +
  xlab("iDecade") + ylab("Species richness")

#Calculate loss in full dataset during last two idecades relative to the average of previous idecades
1-chaoasyest_global_bees_all_idec[6,4]/mean(chaoasyest_global_bees_all_idec[1:5,4]) # Percent loss during the 2000s
1-chaoasyest_global_bees_all_idec[7,4]/mean(chaoasyest_global_bees_all_idec[1:5,4]) # Percent loss during the 2010s

#### Figure 1 --------
grid.arrange(
  plot_recs,
  plot_sp,
  plot_chao_prsv,plot_est_prsv,
  plot_chao_all,plot_est_all,  
  widths = c(1, 1, 1),
  layout_matrix = rbind(c(1, 1, 1),
                        c(2, 2, 2),
                        c(3, 3, 4),
                        c(5, 5, 6)))


#### Calculate stats for each family --------

#Re-read data from GBIF data file (WARNING: change nThread parameter to an appropriate value) for Anthophila+Outgroups, filtering out records without a species id
bwa_all <-fread(GBIFdata, select=fields, nThread=6) %>% filter(family %in% all_fams, taxonRank=="SPECIES", year!="", species!="")
bwa_all$idecade <- floor((bwa_all$year + 4)/10)*10
bwa_all$family <- factor(bwa_all$family,levels=all_fams)
bwa_prsv <- bwa_all %>% filter(basisOfRecord=="PRESERVED_SPECIMEN")

#Total Records per Family
records_per_fam <- bwa_all %>% group_by(family) %>% count()

#Count number of records by species and year (for complete and preserved specimens only datasets)
bwa_all_recs_by_sp_year <- bwa_all %>% group_by(species, family, year, idecade) %>% count()
bwa_prsv_recs_by_sp_year <- bwa_prsv %>% group_by(species, family, year, idecade) %>% count()

bwa_all_recs_year <- bwa_all_recs_by_sp_year %>% group_by(family, year) %>% count() %>% rename(records = n)
bwa_prsv_recs_year <- bwa_prsv_recs_by_sp_year %>% group_by(family, year) %>% count() %>% rename(records = n)

#Plot number of records per year, by family, for full (blue) and specimen-only (red) datasets
ggplot() +
  geom_point(data = bwa_all_recs_year %>% filter(year<2016, year>1945), aes(x = year, y = records), color = "#00AFBB") +
  geom_smooth(data = bwa_all_recs_year %>% filter(year<2016, year>1945), aes(x = year, y = records), color = "#00AFBB", fill = "lightblue") +
  geom_point(data = bwa_prsv_recs_year %>% filter(year<2016, year>1945), aes(x = year, y = records), color = "red") +
  geom_smooth(data = bwa_prsv_recs_year %>% filter(year<2016, year>1945), aes(x = year, y = records), color = "red", fill = "pink") +
  facet_grid(family~., scales = "free") 


#Count number of species by year, grouped by family (for complete and preserved specimens only datasets)
bwa_all_sp_year <- bwa_all_recs_by_sp_year  %>% select(-n)  %>% group_by(family, year) %>% count() %>% rename(sp = n)
bwa_prsv_sp_year <- bwa_prsv_recs_by_sp_year  %>% select(-n)  %>% group_by(family, year) %>% count() %>% rename(sp = n)


all_famcols=c("gray","blue","blue","violet","red","red","red","darkorange","darkorange")
names(all_famcols)<- all_fams

#Plot species per year, by family, specimen-only dataset
ggplot(bwa_prsv_sp_year %>% filter(year<2016, year>1945), aes(x = year, y = sp, color = family)) +
  geom_point(show.legend = F)+
  geom_smooth(show.legend = F)+
  facet_grid(family~., scales = "free_y")+
  scale_color_manual(values = all_famcols)

#Plot species per year, by family, full dataset
ggplot(bwa_all_sp_year %>% filter(year<2016, year>1945), aes(x = year, y = sp, color = family)) +
  geom_point(show.legend = F)+
  geom_smooth(show.legend = F)+
  facet_grid(family~., scales = "free_y")+
  scale_color_manual(values = all_famcols)

####Calculate Chao statistics for each FAMILY####

####Count number of records by species and idecade (for complete and preserved specimens only datasets)----------------
bwa_all_recs_by_sp_idec <- bwa_all %>% group_by(species, family, idecade) %>% count()
bwa_prsv_recs_by_sp_idec <- bwa_prsv %>% group_by(species, family, idecade) %>% count()

#Then create an empty list to hold the iNEXT output for each family and run iNEXT in each family
bwa_all_idec_chao <- list()
for (i in all_fams)
{
  famdata <- with(bwa_all_recs_by_sp_idec %>% filter(idecade>1940, idecade<2020, family == i), as.data.frame(tapply(n, list(species,idecade),sum, na.rm=T)))
  famdata[is.na(famdata)] <- 0
  bwa_all_idec_chao[[i]] <- iNEXT(famdata, q=0, datatype = "abundance")
}

bwa_prsv_idec_chao <- list()
for (i in all_fams)
{
  famdata <- with(bwa_prsv_recs_by_sp_idec %>% filter(idecade>1940, idecade<2020, family == i), as.data.frame(tapply(n, list(species,idecade),sum, na.rm=T)))
  famdata[is.na(famdata)] <- 0
  bwa_prsv_idec_chao[[i]] <- iNEXT(famdata, q=0, datatype = "abundance")
}

for (i in all_fams){ 
  p <- 
  ggiNEXT(bwa_prsv_idec_chao[[i]], type=1) +
  scale_shape_manual(values=c(16,16,16,16,16,16,16)) +
  scale_size_manual(values = rep(2,9)) +
  scale_colour_manual(values=viridis(length(levels(bwa_prsv_idec_chao[[i]]$AsyEst$Site)),direction=-1)) +
  theme(legend.position="none", axis.title.x = element_blank(),axis.title.y = element_blank(), text = element_text(size = 10))

 
  assign(paste0("chao_prsv_plot_", i), p)
}

grid.arrange(chao_prsv_plot_Formicidae,
             chao_prsv_plot_Sphecidae,
             chao_prsv_plot_Crabronidae,
             chao_prsv_plot_Melittidae,
             chao_prsv_plot_Andrenidae,
             chao_prsv_plot_Halictidae,
             chao_prsv_plot_Colletidae,
             chao_prsv_plot_Megachilidae,
             chao_prsv_plot_Apidae,
             ncol = 1)


for (i in all_fams){
  asyest <- bwa_prsv_idec_chao[[i]]$AsyEst[seq(1,21,3),]
  asyest$decade <- as.numeric(levels(asyest$Site))[asyest$Site] #convert factor into numeric variable
  plot(asyest$decade, asyest$Estimator, type="p", ylim=c(min(asyest$LCL)*0.9,max(asyest$UCL)*1.1), pch=16, ylab="Estimated richness", xlab="", main=i)
  lines(asyest$decade, asyest$Estimator)
  arrows(asyest$decade, asyest$LCL, asyest$decade, asyest$UCL, length=0.05, angle=90, code=3)
}

for (i in all_fams){ 
  asyest <- bwa_prsv_idec_chao[[i]]$AsyEst[seq(1,21,3),]
  p <- 
  ggplot(asyest, aes(Site, Estimator)) +
  geom_bar(stat="Identity") +
  coord_cartesian(ylim = c(0.9*min(asyest$LCL), max(asyest$UCL))) + 
  geom_errorbar(aes(ymin=LCL, ymax=UCL), width=.2, position=position_dodge(.9)) +
  theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(), text = element_text(size = 8))
  
  assign(paste0("chao_prsv_est_", i), p)
}

grid.arrange(chao_prsv_est_Formicidae,
             chao_prsv_est_Sphecidae,
             chao_prsv_est_Crabronidae,
             chao_prsv_est_Melittidae,
             chao_prsv_est_Andrenidae,
             chao_prsv_est_Halictidae,
             chao_prsv_est_Colletidae,
             chao_prsv_est_Megachilidae,
             chao_prsv_est_Apidae,
             ncol = 1)



#### Plot all family based plots -----

#Plot species per year, by family, specimen-only dataset
famplot <- ggplot(bwa_prsv_sp_year %>% filter(year<2016, year>1945), aes(x = year, y = sp, color = family)) +
  geom_point(show.legend = F)+
  geom_smooth(show.legend = F)+
  facet_grid(family~., scales = "free_y")+
  scale_color_manual(values = all_famcols) +
  theme(legend.position="none", axis.title.x = element_blank(),axis.title.y = element_blank(), text = element_text(size = 10))

#### Figure 2 ------
grid.arrange(famplot,
             chao_prsv_plot_Formicidae, #2
             chao_prsv_plot_Sphecidae,
             chao_prsv_plot_Crabronidae,
             chao_prsv_plot_Melittidae,
             chao_prsv_plot_Andrenidae,
             chao_prsv_plot_Halictidae,
             chao_prsv_plot_Colletidae,
             chao_prsv_plot_Megachilidae,
             chao_prsv_plot_Apidae,
             chao_prsv_est_Formicidae, #11
             chao_prsv_est_Sphecidae,
             chao_prsv_est_Crabronidae,
             chao_prsv_est_Melittidae,
             chao_prsv_est_Andrenidae,
             chao_prsv_est_Halictidae,
             chao_prsv_est_Colletidae,
             chao_prsv_est_Megachilidae,
             chao_prsv_est_Apidae,
             widths = c(2, 2, 1),
             layout_matrix = rbind(c(1, 2, 11),
                                   c(1, 3, 12),
                                   c(1, 4, 13),
                                   c(1, 5, 14),
                                   c(1, 6, 15),
                                   c(1, 7, 16),
                                   c(1, 8, 17),
                                   c(1, 9, 18),
                                   c(1, 10,19)))

#Calculate loss in specimen dataset during last two idecades relative to the average of previous idecades

for(i in all_fams){
  print(
    paste(i, "loss during the 2000s:",
          (1- round( (bwa_prsv_idec_chao[[i]]$AsyEst %>% filter(Diversity=="Species richness"))[6,"Estimator"] /
                       mean((bwa_prsv_idec_chao[[i]]$AsyEst %>% filter(Diversity=="Species richness"))[1:5,"Estimator"]),2))*100,
          "%"))
  print(
    paste(i, "loss during the 2010s:",
          (1- round( (bwa_prsv_idec_chao[[i]]$AsyEst %>% filter(Diversity=="Species richness"))[7,"Estimator"] /
                       mean((bwa_prsv_idec_chao[[i]]$AsyEst %>% filter(Diversity=="Species richness"))[1:5,"Estimator"]),2))*100,
          "%"))
}



####ITIS species counts vs asymptotic estimators####

#Species counts obtained from the Integrated Taxonomic Information System (ITIS, www.itis.gov).
itis_sp_cts <-c(10213,784,8871,206,3010,4186,2493,5062,5677)
names(itis_sp_cts)<-all_fams

bwa_all_species_family <- bwa_all %>% group_by(species, family) %>% count()

global_chao <- list()
for (i in all_fams)
{
  famdata <- bwa_all_species_family %>% filter(family==i)
  famdata <- as.matrix(famdata[,c(-1,-2)])
  global_chao[[i]] <- iNEXT(famdata, q=0, datatype = "abundance")
}

chao_sp_cts <-as.numeric()
for(i in all_fams) {chao_sp_cts[i]<- global_chao[[i]]$AsyEst[1,4]}
chao_sp_cts <-unlist(chao_sp_cts)

itis_chao <- bind_cols(as_tibble(cbind(itis_sp_cts, chao_sp_cts)), all_famcols)
names(itis_chao) <- c("itis","chao","Family")

#### Figure S1 -----------
ggplot(data= itis_chao, aes(x= itis, y=chao)) +
  geom_point() + xlab ("ITIS Species Count") +
  ylab("Asymptotic Estimator of Richness") +
  geom_text(aes(label = all_fams), nudge_y= 200) +
  geom_abline(slope = 1)

####Counting fraction of species with no ID####
#Filter anthophila families from unfiltered original GBIF data
bees_raw <-fread(GBIFdata, select=fields, nThread=6) %>% filter(family %in% anthophila_fams, year!="")

total_recs_year <- bees_raw %>% group_by(year) %>% count() %>% rename(all = n)
missid_year <- bees_raw %>% filter(taxonRank %in% c("GENUS","FAMILY", "UNRANKED")) %>% group_by(year) %>% count() %>% rename(misid = n)

missid_year <- right_join(missid_year, total_recs_year, by="year")
missid_year$fraction <- missid_year$misid/missid_year$all

#### Figure S2 -----
ggplot(missid_year %>% filter(year>1899, year<2016), aes(x = year, y = fraction)) +
  geom_point() +
  geom_smooth() +
  ylab("Fraction of record missing species-rank ID") + xlab("Year")

#### Estimate trends for Anthophila per continent####
continents <- c("Africa","Asia","Europe","North America","Oceania","South America")

#Count records per decade and continent
bees_all_cont <-left_join(bees_all %>% filter(countryCode!=""), country2continent, by="countryCode") 
bees_all_cont_recs_by_sp_year <- bees_all_cont %>% group_by(species, Continent, year, idecade) %>% summarize(records = length(gbifID))
bees_prsv_cont_recs_by_sp_year <- bees_all_cont %>% filter(basisOfRecord=="PRESERVED_SPECIMEN")  %>% group_by(species, Continent, year, idecade) %>% summarize(records = length(gbifID))

#Count records per continent
bees_all_recs_by_cont <- bees_all_cont %>% filter(Continent %in% continents) %>% group_by(Continent, idecade) %>% summarize(records = length(gbifID))

#Plot Records By continent
cont_stack <- ggplot(data=bees_all_recs_by_cont %>% filter(idecade>1940, idecade<2020), aes(fill= Continent, y=records, x = idecade))+
  geom_bar(position="stack", stat="identity", show.legend = F) +
  xlab("iDecade")
cont_perc <- ggplot(data=bees_all_recs_by_cont %>% filter(idecade>1940, idecade<2020), aes(fill= Continent, y=records, x = idecade))+
  geom_bar(position="fill", stat="identity") +
  xlab("iDecade") + ylab("Percent of decadal records")

#### Figure S3 -----
grid.arrange(cont_stack, cont_perc, nrow = 1)


#Count number of species by year, grouped by family (for complete and preserved specimens only datasets)
bees_all_cont_sp_year <- bees_all_cont_recs_by_sp_year  %>% select(-records)  %>% group_by(Continent, year) %>% count() %>% rename(sp = n)
bees_all_cont_recs_year <- bees_all_cont_recs_by_sp_year %>% group_by(Continent, year) %>% summarize(records = sum(records))

bees_prsv_cont_sp_year <- bees_prsv_cont_recs_by_sp_year  %>% select(-records)  %>% group_by(Continent, year) %>% count() %>% rename(sp = n)
bees_prsv_cont_recs_year <- bees_prsv_cont_recs_by_sp_year %>% group_by(Continent, year) %>% summarize(records = sum(records))

cont_recs <- 
ggplot() +
  geom_point(data = bees_all_cont_recs_year %>% filter(year<2016, year>1945, !is.na(Continent)), aes(x = year, y = records), color = "#00AFBB") +
  geom_smooth(data = bees_all_cont_recs_year %>% filter(year<2016, year>1945, !is.na(Continent)), aes(x = year, y = records), color = "#00AFBB", fill = "lightblue") +
  geom_point(data = bees_prsv_cont_recs_year %>% filter(year<2016, year>1945, !is.na(Continent)), aes(x = year, y = records), color = "red") +
  geom_smooth(data = bees_prsv_cont_recs_year %>% filter(year<2016, year>1945, !is.na(Continent)), aes(x = year, y = records), color = "red", fill = "pink") +
  facet_grid(Continent~., scales = "free_y") +
  theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank(), text = element_text(size = 8), strip.background = element_blank(),strip.text.y = element_blank())
 

cont_sp <- ggplot() +
  geom_point(data = bees_all_cont_sp_year %>% filter(year<2016, year>1945, !is.na(Continent)), aes(x = year, y = sp), color = "#00AFBB") +
  geom_smooth(data = bees_all_cont_sp_year %>% filter(year<2016, year>1945, !is.na(Continent)), aes(x = year, y = sp), color = "#00AFBB", fill = "lightblue") +
  geom_point(data = bees_prsv_cont_sp_year %>% filter(year<2016, year>1945, !is.na(Continent)), aes(x = year, y = sp), color = "red") +
  geom_smooth(data = bees_prsv_cont_sp_year %>% filter(year<2016, year>1945, !is.na(Continent)), aes(x = year, y = sp), color = "red", fill = "pink") +
  facet_grid(Continent~., scales = "free_y") +
  theme(legend.position="none", axis.title.x = element_blank(),axis.title.y = element_blank(), text = element_text(size = 8), strip.background = element_blank(),strip.text.y = element_blank())


####Calculate Chao statistics for each CONTINENT ####


## FULL dataset 
#Create an empty list to hold the iNEXT output for each family and run iNEXT in each family
cont_all_idec_chao <- list()
for (i in continents)
{
  contdata <- with(bees_all_cont_recs_by_sp_year %>% filter(idecade>1940, idecade<2020, Continent == i), as.data.frame(tapply(records, list(species,idecade),sum, na.rm=T)))
  contdata[is.na(contdata)] <- 0
  cont_all_idec_chao[[i]] <- iNEXT(contdata, q=0, datatype = "abundance")
}

# Generate ggiNEXT plots for each continent
for (i in continents){ 
  p <- 
    ggiNEXT(cont_all_idec_chao[[i]], type=1) +
    scale_shape_manual(values=c(16,16,16,16,16,16,16)) +
    scale_size_manual(values = rep(2,9)) +
    scale_colour_manual(values=viridis(length(levels(cont_prsv_idec_chao[[i]]$AsyEst$Site)),direction=-1)) +
    theme(legend.position="none", axis.title.x = element_blank(),axis.title.y = element_blank(), text = element_text(size = 10))
  assign(paste0("chao_all_plot_", sub(' ', '', i)), p)
}

# Arrange and display ggiNEXT plots for each continent
grid.arrange(chao_all_plot_Africa,
             chao_all_plot_Asia,
             chao_all_plot_Europe,
             chao_all_plot_NorthAmerica,
             chao_all_plot_Oceania,
             chao_all_plot_SouthAmerica,
             ncol = 1)

# Extract idecadal richness estimator statistics for each continent and generate barplots
for (i in continents){ 
  asyest <- cont_prsv_idec_chao[[i]]$AsyEst[seq(1,21,3),]
  p <- 
    ggplot(asyest, aes(Site, Estimator)) +
    geom_bar(stat="Identity") +
    coord_cartesian(ylim = c(0.9*min(asyest$LCL), max(asyest$UCL))) + 
    geom_errorbar(aes(ymin=LCL, ymax=UCL), width=.2, position=position_dodge(.9)) +
    theme(legend.position="none", axis.title.x = element_blank(),axis.title.y = element_blank(), text = element_text(size = 6))
  assign(paste0("chao_prsv_est_",  sub(' ', '', i)), p)
}

# Arrange and display barplots for each continent
grid.arrange(chao_prsv_est_Africa,
             chao_prsv_est_Asia,
             chao_prsv_est_Europe,
             chao_prsv_est_NorthAmerica,
             chao_prsv_est_Oceania,
             chao_prsv_est_SouthAmerica,
             ncol = 1)

## SPECIMEN dataset 
#Create an empty list to hold the iNEXT output for each family and run iNEXT in each family
cont_prsv_idec_chao <- list()
for (i in continents)
{
  contdata <- with(bees_prsv_cont_recs_by_sp_year %>% filter(idecade>1940, idecade<2020, Continent == i), as.data.frame(tapply(records, list(species,idecade),sum, na.rm=T)))
  contdata[is.na(contdata)] <- 0
  cont_prsv_idec_chao[[i]] <- iNEXT(contdata, q=0, datatype = "abundance")
}

# Generate ggiNEXT plots for each continent
for (i in continents){ 
  p <- 
    ggiNEXT(cont_prsv_idec_chao[[i]], type=1) +
    scale_shape_manual(values=c(16,16,16,16,16,16,16)) +
    scale_size_manual(values = rep(2,9)) +
    scale_colour_manual(values=viridis(length(levels(cont_prsv_idec_chao[[i]]$AsyEst$Site)),direction=-1)) +
    theme(legend.position="none", axis.title.x = element_blank(),axis.title.y = element_blank(), text = element_text(size = 10))
  assign(paste0("chao_prsv_plot_", sub(' ', '', i)), p)
}

# Arrange and display ggiNEXT plots for each continent
grid.arrange(chao_prsv_plot_Africa,
             chao_prsv_plot_Asia,
             chao_prsv_plot_Europe,
             chao_prsv_plot_NorthAmerica,
             chao_prsv_plot_Oceania,
             chao_prsv_plot_SouthAmerica,
             ncol = 1)

# Extract idecadal richness estimator statistics for each continent and generate barplots
for (i in continents){ 
  asyest <- cont_prsv_idec_chao[[i]]$AsyEst[seq(1,21,3),]
  p <- 
    ggplot(asyest, aes(Site, Estimator)) +
    geom_bar(stat="Identity") +
    coord_cartesian(ylim = c(0.9*min(asyest$LCL), max(asyest$UCL))) + 
    geom_errorbar(aes(ymin=LCL, ymax=UCL), width=.2, position=position_dodge(.9)) +
    theme(legend.position="none", axis.title.x = element_blank(),axis.title.y = element_blank(), text = element_text(size = 9))
  assign(paste0("chao_prsv_est_",  sub(' ', '', i)), p)
}

# Arrange and display barplots for each continent
grid.arrange(chao_prsv_est_Africa,
             chao_prsv_est_Asia,
             chao_prsv_est_Europe,
             chao_prsv_est_NorthAmerica,
             chao_prsv_est_Oceania,
             chao_prsv_est_SouthAmerica,
             ncol = 1)

#### Figure S4 ------
## Multipanel figure for continent-level analyses
grid.arrange(cont_recs, cont_sp,
             chao_prsv_plot_Africa, #3
             chao_prsv_plot_Asia,
             chao_prsv_plot_Europe,
             chao_prsv_plot_NorthAmerica,
             chao_prsv_plot_Oceania,
             chao_prsv_plot_SouthAmerica,
             chao_prsv_est_Africa, #9
             chao_prsv_est_Asia,
             chao_prsv_est_Europe,
             chao_prsv_est_NorthAmerica,
             chao_prsv_est_Oceania,
             chao_prsv_est_SouthAmerica,
             widths = c(2, 2, 2, 1),
             layout_matrix = rbind(c(1, 2, 3, 9),
                                   c(1, 2, 4, 10),
                                   c(1, 2, 5, 11),
                                   c(1, 2, 6, 12),
                                   c(1, 2, 7, 13),
                                   c(1, 2, 8, 14)))


####Pielou's eveness index over time####
#ESTIMATE GLOBAL EQUITATIVITY USING PIELOU'S J EVENNESS

#FULL dataset
bees_all_species_year_table <- xtabs(~year + species, data = bees_all %>% filter(year>1899, year<2016))
year <- bees_all %>% filter(year>1899, year<2016) %>% select(year) %>% unique()
H_all <- diversity(bees_all_species_year_table)
J_all <- H_all/log(specnumber(bees_all_species_year_table))

years <- names(J_all)
J_all_tbl <- as_tibble(cbind(as.numeric(years), as.numeric(J_all)))
colnames(J_all_tbl)<-c("year", "J")

#SPECIMEN dataset
bees_prsv_species_year_table <- xtabs(~year + species, data = bees_all %>% filter(year>1899, year<2016, basisOfRecord=="PRESERVED_SPECIMEN"))
year_prsv <- bees_all %>% filter(year>1899, year<2016, basisOfRecord=="PRESERVED_SPECIMEN") %>% select(year) %>% unique()
H_prsv <- diversity(bees_prsv_species_year_table)
J_prsv <- H_prsv/log(specnumber(bees_prsv_species_year_table))

years <- names(J_prsv)
J_prsv_tbl <- as_tibble(cbind(as.numeric(years), as.numeric(J_prsv)))
colnames(J_prsv_tbl)<-c("year", "J")

####Figure 3 -----
#Plot J over time for both datasets
ggplot() +
  geom_point(data = J_all_tbl, aes(x = year, y = J), color = "#00AFBB") +
  geom_smooth(data = J_all_tbl, aes(x = year, y = J), color = "#00AFBB", fill = "lightblue") +
  geom_point(data = J_prsv_tbl, aes(x = year, y = J), color = "red") +
  geom_smooth(data = J_prsv_tbl, aes(x = year, y = J), color = "red", fill = "pink") +
  xlab("Year")
  

####Calculate fraction of Apis mellifera records over time for both datasets####
apis_all_recs_by_year <- right_join(bees_all_recs_year, bees_all_recs_by_sp_year %>% filter(species=="Apis mellifera"), by="year")
apis_all_recs_by_year <- apis_all_recs_by_year %>% mutate(fraction = n / records)

apis_prsv_recs_by_year <- right_join(bees_prsv_recs_year, bees_prsv_recs_by_sp_year %>% filter(species=="Apis mellifera"), by="year")
apis_prsv_recs_by_year <- apis_prsv_recs_by_year %>% mutate(fraction = n / records)

####Figure S5####
ggplot() +
  geom_point(data = apis_all_recs_by_year %>% filter(year>1899, year<2016), aes(x = year, y = fraction), color = "#00AFBB") +
  geom_smooth(data = apis_all_recs_by_year %>% filter(year>1899, year<2016), aes(x = year, y = fraction), color = "#00AFBB", fill = "lightblue") +
  geom_point(data = apis_prsv_recs_by_year %>% filter(year>1899, year<2016), aes(x = year, y = fraction), color = "red") +
  geom_smooth(data = apis_prsv_recs_by_year %>% filter(year>1899, year<2016), aes(x = year, y = fraction), color = "red", fill = "pink")+
  xlab("Year") + ylab("Fraction of total records reporting Apis mellifera")

#### Country specific observations ####

####United States of America ----
bees_us_recs_by_sp_year <- bees_all %>% filter(countryCode=="US", basisOfRecord=="PRESERVED_SPECIMEN") %>% group_by(species, family, year, idecade) %>% count()
bees_us_recs_year <- bees_all  %>% filter(countryCode=="US", basisOfRecord=="PRESERVED_SPECIMEN") %>% group_by(year) %>% count() %>% rename(records = n)
bees_us_sp_year <- bees_us_recs_by_sp_year %>% select(-n)  %>% group_by(year) %>% count() %>% rename(sp = n)

bees_us_recs_sp_by_year <- inner_join(bees_us_recs_year, bees_us_sp_year, by="year")

plot_us_recs <- 
  ggplot(data = bees_us_recs_sp_by_year %>% filter(year>1899, year<2016), aes(x = year, y = records)) +
  geom_point() +
  geom_smooth() +
  xlab("Year") + ylab("No. records")
plot_us_sp <- 
ggplot(data = bees_us_recs_sp_by_year %>% filter(year>1899, year<2016), aes(x = year, y = sp)) +
  geom_point() +
  geom_smooth() +
  xlab("Year") + ylab("Species richness")

sum(bees_us_recs_sp_by_year$records) # Total number of records

#Count number of records by species and idecade
bees_us_recs_by_sp_idec <- bees_all %>% filter(countryCode=="US", basisOfRecord=="PRESERVED_SPECIMEN") %>% group_by(species, family, idecade) %>% count()

#Generate idecadal matrices and calculate Chao stats
bees_us_recs_by_sp_idec.mat <-with(bees_us_recs_by_sp_idec %>% filter(idecade>1940, idecade<2020), as.data.frame(tapply(n, list(species,idecade),sum, na.rm=T)))
bees_us_recs_by_sp_idec.mat[is.na(bees_us_recs_by_sp_idec.mat)] <- 0
bees_us_chao_idec <- iNEXT(bees_us_recs_by_sp_idec.mat, q=0, datatype = "abundance")

#iNEXT curves
plot_chao_us <- 
  ggiNEXT(bees_us_chao_idec, type=1) +
  scale_shape_manual(values=c(19,19,19,19,19,19,19)) + xlim(0,4e+05) +
  scale_colour_manual(values=viridis(length(levels(bees_us_chao_idec$AsyEst$Site)),direction=-1)) +
  theme(legend.position="none", axis.text=element_text(size=10),axis.title=element_text(size=12) ) + 
  xlab("No. records") + ylab("Species richness")

#Extract and plot estimators as barplots
chaoasyest_us_bees_idec <- bees_us_chao_idec$AsyEst[seq(1,21,3),]

plot_est_us <- 
  ggplot(chaoasyest_us_bees_idec, aes(Site, Estimator)) +
  geom_bar(stat="Identity") +
  coord_cartesian(ylim = c(min(chaoasyest_us_bees_idec$LCL*0.5), max(chaoasyest_us_bees_idec$UCL))) + 
    geom_errorbar(aes(ymin=LCL, ymax=UCL), width=.2,
                position=position_dodge(.9)) +
  xlab("idecade")  + ylab("Species richness")



####Brazil----
bees_br_recs_by_sp_year <- bees_all %>% filter(countryCode=="BR", basisOfRecord=="PRESERVED_SPECIMEN") %>% group_by(species, family, year, idecade) %>% count()
bees_br_recs_year <- bees_all  %>% filter(countryCode=="BR", basisOfRecord=="PRESERVED_SPECIMEN") %>% group_by(year) %>% count() %>% rename(records = n)
bees_br_sp_year <- bees_br_recs_by_sp_year %>% select(-n)  %>% group_by(year) %>% count() %>% rename(sp = n)

bees_br_recs_sp_by_year <- inner_join(bees_br_recs_year, bees_br_sp_year, by="year")

plot_br_recs <- 
  ggplot(data = bees_br_recs_sp_by_year %>% filter(year>1899, year<2016), aes(x = year, y = records)) +
  geom_point() +
  geom_smooth() +
  xlab("Year") + ylab("No. records")
plot_br_sp <- 
  ggplot(data = bees_br_recs_sp_by_year %>% filter(year>1899, year<2016), aes(x = year, y = sp)) +
  geom_point() +
  geom_smooth() +
  xlab("Year") + ylab("Species richness")

sum(bees_br_recs_sp_by_year$records) # Total number of records

#Count number of records by species and idecade
bees_br_recs_by_sp_idec <- bees_all %>% filter(countryCode=="BR", basisOfRecord=="PRESERVED_SPECIMEN") %>% group_by(species, family, idecade) %>% count()

#Generate idecadal matrices and calculate Chao stats
bees_br_recs_by_sp_idec.mat <-with(bees_br_recs_by_sp_idec %>% filter(idecade>1940, idecade<2020), as.data.frame(tapply(n, list(species,idecade),sum, na.rm=T)))
bees_br_recs_by_sp_idec.mat[is.na(bees_br_recs_by_sp_idec.mat)] <- 0
bees_br_chao_idec <- iNEXT(bees_br_recs_by_sp_idec.mat, q=0, datatype = "abundance")

#iNEXT curves
plot_chao_br <- 
  ggiNEXT(bees_br_chao_idec, type=1) +
  scale_shape_manual(values=c(19,19,19,19,19,19,19)) + 
  scale_colour_manual(values=viridis(length(levels(bees_br_chao_idec$AsyEst$Site)),direction=-1)) +
  theme(legend.position="none", axis.text=element_text(size=10),axis.title=element_text(size=12) ) + 
  xlab("No. records") + ylab("Species richness")

#Extract and plot estimators as barplots
chaoasyest_br_bees_idec <- bees_br_chao_idec$AsyEst[seq(1,21,3),]

plot_est_br <- 
  ggplot(chaoasyest_br_bees_idec, aes(Site, Estimator)) +
  geom_bar(stat="Identity") +
  coord_cartesian(ylim = c(min(chaoasyest_br_bees_idec$LCL*0.5), max(chaoasyest_br_bees_idec$UCL))) + 
  geom_errorbar(aes(ymin=LCL, ymax=UCL), width=.2,
                position=position_dodge(.9)) +
  xlab("idecade")  + ylab("Species richness")

####Great Britain ----
bees_gb_recs_by_sp_year <- bees_all %>% filter(countryCode=="GB", basisOfRecord=="PRESERVED_SPECIMEN") %>% group_by(species, family, year, idecade) %>% count()
bees_gb_recs_year <- bees_all  %>% filter(countryCode=="GB", basisOfRecord=="PRESERVED_SPECIMEN") %>% group_by(year) %>% count() %>% rename(records = n)
bees_gb_sp_year <- bees_gb_recs_by_sp_year %>% select(-n)  %>% group_by(year) %>% count() %>% rename(sp = n)

bees_gb_recs_sp_by_year <- inner_join(bees_gb_recs_year, bees_gb_sp_year, by="year")

plot_gb_recs <- 
  ggplot(data = bees_gb_recs_sp_by_year %>% filter(year>1899, year<2016), aes(x = year, y = records)) +
  geom_point() +
  geom_smooth() +
  xlab("Year") + ylab("No. records")
plot_gb_sp <- 
  ggplot(data = bees_gb_recs_sp_by_year %>% filter(year>1899, year<2016), aes(x = year, y = sp)) +
  geom_point() +
  geom_smooth() +
  xlab("Year") + ylab("Species richness")

sum(bees_gb_recs_sp_by_year$records) # Total number of records

#Count number of records by species and idecade
bees_gb_recs_by_sp_idec <- bees_all %>% filter(countryCode=="GB", basisOfRecord=="PRESERVED_SPECIMEN") %>% group_by(species, family, idecade) %>% count()

#Generate idecadal matrices and calculate Chao stats
bees_gb_recs_by_sp_idec.mat <-with(bees_gb_recs_by_sp_idec %>% filter(idecade>1940, idecade<2020), as.data.frame(tapply(n, list(species,idecade),sum, na.rm=T)))
bees_gb_recs_by_sp_idec.mat[is.na(bees_gb_recs_by_sp_idec.mat)] <- 0
bees_gb_chao_idec <- iNEXT(bees_gb_recs_by_sp_idec.mat, q=0, datatype = "abundance")

#iNEXT curves
plot_chao_gb <- 
  ggiNEXT(bees_gb_chao_idec, type=1) +
  scale_shape_manual(values=c(19,19,19,19,19,19,19)) + 
  scale_colour_manual(values=viridis(length(levels(bees_gb_chao_idec$AsyEst$Site)),direction=-1)) +
  theme(legend.position="none", axis.text=element_text(size=10),axis.title=element_text(size=12) ) + 
  xlab("No. records") + ylab("Species richness")

#Extract and plot estimators as barplots
chaoasyest_gb_bees_idec <- bees_gb_chao_idec$AsyEst[seq(1,21,3),]

plot_est_gb <- 
  ggplot(chaoasyest_gb_bees_idec, aes(Site, Estimator)) +
  geom_bar(stat="Identity") +
  coord_cartesian(ylim = c(min(chaoasyest_gb_bees_idec$LCL*0.5), max(chaoasyest_gb_bees_idec$UCL))) + 
  geom_errorbar(aes(ymin=LCL, ymax=UCL), width=.2,
                position=position_dodge(.9)) +
  xlab("idecade")  + ylab("Species richness")


####Panama  ----
bees_pa_recs_by_sp_year <- bees_all %>% filter(countryCode=="PA", basisOfRecord=="PRESERVED_SPECIMEN") %>% group_by(species, family, year, idecade) %>% count()
bees_pa_recs_year <- bees_all  %>% filter(countryCode=="PA", basisOfRecord=="PRESERVED_SPECIMEN") %>% group_by(year) %>% count() %>% rename(records = n)
bees_pa_sp_year <- bees_pa_recs_by_sp_year %>% select(-n)  %>% group_by(year) %>% count() %>% rename(sp = n)

bees_pa_recs_sp_by_year <- inner_join(bees_pa_recs_year, bees_pa_sp_year, by="year")

plot_pa_recs <- 
  ggplot(data = bees_pa_recs_sp_by_year %>% filter(year>1899, year<2016), aes(x = year, y = records)) +
  geom_point() +
  geom_smooth() +
  xlab("Year") + ylab("No. records")
plot_pa_sp <- 
  ggplot(data = bees_pa_recs_sp_by_year %>% filter(year>1899, year<2016), aes(x = year, y = sp)) +
  geom_point() +
  geom_smooth() +
  xlab("Year") + ylab("Species richness")

sum(bees_pa_recs_sp_by_year$records) # Total number of records

#Count number of records by species and idecade
bees_pa_recs_by_sp_idec <- bees_all %>% filter(countryCode=="PA", basisOfRecord=="PRESERVED_SPECIMEN") %>% group_by(species, family, idecade) %>% count()

#Generate idecadal matrices and calculate Chao stats
bees_pa_recs_by_sp_idec.mat <-with(bees_pa_recs_by_sp_idec %>% filter(idecade>1940, idecade<2020), as.data.frame(tapply(n, list(species,idecade),sum, na.rm=T)))
bees_pa_recs_by_sp_idec.mat[is.na(bees_pa_recs_by_sp_idec.mat)] <- 0
bees_pa_chao_idec <- iNEXT(bees_pa_recs_by_sp_idec.mat, q=0, datatype = "abundance")

#iNEXT curves
plot_chao_pa <- 
  ggiNEXT(bees_pa_chao_idec, type=1) +
  scale_shape_manual(values=c(19,19,19,19,19,19,19)) + 
  scale_colour_manual(values=viridis(length(levels(bees_pa_chao_idec$AsyEst$Site)),direction=-1)) +
  theme(legend.position="none", axis.text=element_text(size=10),axis.title=element_text(size=12) ) + 
  xlab("No. records") + ylab("Species richness")

#Extract and plot estimators as barplots
chaoasyest_pa_bees_idec <- bees_pa_chao_idec$AsyEst[seq(1,21,3),]

plot_est_pa <- 
  ggplot(chaoasyest_pa_bees_idec, aes(Site, Estimator)) +
  geom_bar(stat="Identity") +
  coord_cartesian(ylim = c(min(chaoasyest_pa_bees_idec$LCL*0.5), max(chaoasyest_pa_bees_idec$UCL))) + 
  geom_errorbar(aes(ymin=LCL, ymax=UCL), width=.2,
                position=position_dodge(.9)) +
  xlab("idecade")  + ylab("Species richness")

#### Figure S6 ####
## Multipanel figure for continent-level analyses
grid.arrange(plot_us_recs, plot_us_sp,plot_chao_us,plot_est_us, 
             plot_br_recs, plot_br_sp,plot_chao_br,plot_est_br,
             plot_gb_recs, plot_gb_sp,plot_chao_gb,plot_est_gb,
             plot_pa_recs, plot_pa_sp,plot_chao_pa,plot_est_pa,
             widths = c(2, 2, 2, 2),
             layout_matrix = rbind(c(1, 2, 3, 4),
                                   c(5, 6, 7, 8),
                                   c(9, 10, 11, 12),
                                   c(13, 14, 15, 16)))

######################### END OF SCRIPT ---------