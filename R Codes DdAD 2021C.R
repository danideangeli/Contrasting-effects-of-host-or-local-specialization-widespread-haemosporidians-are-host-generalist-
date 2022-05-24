#### Code to run the analyses of the study: "Contrasting effects of host or local specialization: widespread haemosporidians are host generalist whereas local specialists are locally abundant", DeAngeli et al. 2021 - Global Ecology and Biogeography

#### Packages ####

library(ape)
library(dplyr)
library(fossil)
library(hillR)
library(data.table)
library(GeoRange)
library(tibble)
library(reshape2)
library(vegan)
library(ape)
library(ips)
library(MCMCtreeR)
library(treeio)
library(phytools)
library(caper)
library(ggplot2)

##### Creating Dataframes for analyses ######

DF <- read.csv2("SupTable1.csv") #importing data filtered for duplicated sequences

# To run analyses for only Plasmodium/Haemoproteus genera or South America and Europe only it is necessary to filter the data here

DF <- filter(DF, parasiteGenus == "Plasmodium")
DF <- filter(DF, parasiteGenus == "Haemoproteus")
DF <- filter(DF, continent == "Europe")
DF <- filter(DF, continent == "South_America")

#ps: run the script five times using a distinct dataframe each time

# Creating a dataframe with the number of observations of lineages per host species

Int_list = DF %>%
  group_by(Lineage) %>%
  group_by(species, add = TRUE) %>%
  summarise(abundance = n())
Int_list <- data.frame(Int_list)

# Creating a dataframe with the number of localities the host species were observed, and with the total number of occurrences per host species
Occ_list = DF %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(n_sites = n_distinct(site),
                   n_bird_individuals = n())

Occ_list <- as.data.frame(Occ_list)

##### Importing and filtering Bird's Phylogeny ######

# AllBirdsHackett1.tre must be downloaded from https://birdtree.org/
allTrees <- readTree("AllBirdsHackett1.tre")
random_trees <-sample(allTrees@treelst, size = 100)
random_trees1 <- as(random_trees, 'TreeMen')
random_trees2 <- as(random_trees1, 'multiPhylo')

tree <- as.phylo(random_trees2[[1]])

fullbirds <- as.data.frame(tree$tip.label)
mybirds <- as.data.frame(Int_list$species)
names(mybirds) <- c("species")

# Removing from the phylogeny the species which are not present in the MalAvi
todrop<-tree$tip.label[which(tree$tip.label%in%as.character(mybirds[,1])==FALSE)]
mybirds_tree<-drop.tip(tree,todrop)

# Filtering The DF2 for species present in the Phylogeny
Occ_list1<-
  Occ_list %>%
  filter(species %in% mybirds_tree$tip.label)
Occ_list1 <- data.frame(Occ_list1)

# Filtering The DF1 for species present in the Phylogeny

Int_list1 <- 
  Int_list%>%
  filter(species %in% Occ_list1$species)
Int_list1 <- data.frame(Int_list1)

##### Calculating Host's Phylogenetic Range #####

Int_mat <- create.matrix(Int_list1, tax.name ="Lineage", locality ="species", abund = TRUE, abund.col = "abundance")

Phy_Hill <- hill_phylo(Int_mat, mybirds_tree)
Phy_Hill <- as.data.frame(Phy_Hill)

Phy_Hill1 <- Phy_Hill[rownames(Phy_Hill),]
Phy_Hill1 <- as.data.frame(Phy_Hill1)

Phy_Hill1 <- rownames_to_column(Phy_Hill, var = "Lineage") %>% as_tibble()

#### Calculating Geographical Range ####

## creating occurence matrix ##

Occ_mat <- create.matrix(DF, tax.name ="Lineage", locality ="IDLoc", abund = FALSE)
Occ_mat1 <- t(Occ_mat)
Occ_mat1 <- as.data.frame(Occ_mat1) 
setDT(Occ_mat1, keep.rownames = "IDLoc")
Loc_coords <- DF[ , c("IDLoc", "Latitude", "Longitude")] 
Loc_coords1 <- distinct(Loc_coords)
Loc_coords1$IDLoc <- as.character(Loc_coords1$IDLoc)
Occ_mat2 <- inner_join(Occ_mat1, Loc_coords1, by = c("IDLoc" = "IDLoc"))
setcolorder(Occ_mat2, c("Longitude", "Latitude"))
Occ_mat2[Occ_mat2 == 0] <- NA
Occ_mat2$IDLoc <- NULL

Occ_mat2$Longitude <- as.numeric(paste(Occ_mat2$Longitude))
Occ_mat2$Latitude <- as.numeric(paste(Occ_mat2$Latitude))

Geo_Range <- GeoRange_MultiTaxa(OccMatrix= Occ_mat2, TaxaStart=3)
Geo_Range <- rownames_to_column(Geo_Range, var = "Lineage") %>% as_tibble()

Geo_Range1 <- Geo_Range[ , c("Lineage", "MST")] 

## fixing bug from GeoRange function ##

Geo_Range2 <- Geo_Range1 %>%
  mutate(MST=lag(MST, n = 2L))
Geo_Range3 <- Geo_Range1[1601:1602,]
Geo_Range2[1:2,]$MST = Geo_Range3$MST

##### Calculating Parasites' Environmental Range #####

# Crating a df with selected columns to analyzes (Long, Lat, Continent and Country)
coords <- cbind.data.frame(x = DF$Longitude, y = DF$Latitude, Continent = DF$continent, Country = DF$country, ID = DF$IDLoc)
coords <- coords[!duplicated(coords[,c(1,2,5)]),]
coords <- coords[complete.cases(coords),]
coords <- subset(coords, abs(coords$x) <= 180 & abs(coords$y) <= 180)

### Importing worldclim data ###

#dir.create(path = "BioClim") # creating a directory for saving BioClim data
bioclim.data <- raster::getData(name = "worldclim",
                                var = c("bio"),
                                res = 10) # Loading BioClim Data from WorldClim

coords.range <- c(min(coords$x),max(coords$x),min(coords$y),max(coords$y)) # Creating a vector with coords range
points <- SpatialPoints(coords[,c(1,2)], proj4string = bioclim.data@crs) # Creating a spatial object from coords
values <- raster::extract(bioclim.data,points) # Extracting the BioClim values to our Sites coords

# Creating a Dataframe with all variables
BioVars <- cbind.data.frame(coordinates(points),coords[,-c(1,2)],values)
BioVars <- BioVars[complete.cases(BioVars),] # removing sites with NA values for BioClim
rownames(BioVars) <- BioVars$ID

BioVars.st <- cbind.data.frame(BioVars[,c(1:5)],scale(BioVars[,-c(1:5)]))
colnames(BioVars.st)[c(1,2)] <- c("long","lat")

#### Cluster of Localities ####

Dist.Locs <- vegdist(BioVars.st[,-c(1:5)], method = "euclidian")
Clust.Locs <- as.phylo(hclust(Dist.Locs))

Occ_mat3 <- Occ_mat[,Clust.Locs$tip.label]
Occ_mat3 <- Occ_mat3[rowSums(Occ_mat3) != 0,]

## Calculating environmental hill ##

Env_Hill <- hill_phylo(comm = Occ_mat3, tree = Clust.Locs)
Env_Hill <- as.data.frame(Env_Hill)

Env_Hill1 <- rownames_to_column(Env_Hill, var = "Lineage") %>% as_tibble()

#### Joining all data #####

DF$LATLOG <- paste(DF$Latitude, DF$Longitude, sep = "_")

NObsdata <- DF %>%
  group_by(Lineage) %>%
  summarise(NObs = n(),
            NLocs = n_distinct(LATLOG))


DF1 <- left_join(Geo_Range2, Phy_Hill1, by = c("Lineage" = "Lineage"))
DF1 <- left_join(DF1, NObsdata, by = "Lineage")
DF1$Nmediolocal <- DF1$NObs/DF1$NLocs
DF2 <- filter(DF1, NLocs >1)
DF2 <- filter(DF2, Phy_Hill > 0)
DF2 <- data.frame(DF2)
DF2$NObs <- as.numeric(paste(DF2$NObs))
DF2$MST <- as.numeric(paste(DF2$MST))
DF2$Phy_Hill <- as.numeric(paste(DF2$Phy_Hill))
DF2 <- left_join(DF2, Env_Hill1, by =  "Lineage")
DF2 <- filter(DF2, Env_Hill > 0)


#### Importing and filtering Parasite Phylogeny #####


phylo <- read.mrbayes("HaemSeq5.tre")
phylo <- as.phylo(phylo)

# rooting tree
phylo_rooted <- reroot(phylo, node.number = 1)

#filtering parasite phylogeny
matches2 <-match(DF2$Lineage, phylo_rooted$tip.label)
matches2 <-na.omit(matches2)
phylotree <-drop.tip(phylo_rooted, phylo_rooted$tip.label[-matches2])

sptest=as.factor(phylotree$tip.label)
sp=as.data.frame(sptest)

DF3=DF2 %>%
  filter(Lineage %in% sp$sptest)

#### Running PGLS models at Global Scale ####

DF4 <- scale(DF3[,2:7])
DF4 <- as.data.frame(DF4)
DF4$Lineage <- DF3$Lineage 
setcolorder(DF4, c("Lineage"))

DF5 <- comparative.data(phylotree, DF4, Lineage,vcv=TRUE, vcv.dim=3)

mod1 <- pgls(MST ~ Phy_Hill + Nmediolocal + NObs, DF5)
summary.pgls(mod1)

mod2 <- pgls(Env_Hill ~ Phy_Hill + Nmediolocal + NObs, DF5)
summary.pgls(mod2)

#### Testing correlation between Geographical and Environmental Range ####

teste <- cor(DF3$Phy_Hill, DF3$Env_Hill)
teste

#### Plotting results #####

p1 <- ggplot(data = DF3, aes(Phy_Hill,MST)) + labs(x="Phylogenetic Host Range", y="Geographical Range (km)") + 
  geom_point() + geom_point(data = DF3, aes(y = MST), colour = 'red', size = 1.2) +
  geom_smooth(method = "lm", colour = "red")

p1

p1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black"))


p4 <- ggplot(data = DF3, aes(Phy_Hill,Env_Hill)) + labs(x="Phylogenetic Host Range", y="Parasite's Environmental Diversity") + 
  geom_point() + geom_point(data = DF3, aes(y = Env_Hill), colour = 'purple', size = 1.2) +
  geom_smooth(method = "lm", colour = "purple")

p4

p4 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black"))

