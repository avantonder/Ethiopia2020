###############################
### Load required libraries ###
###############################

library(tidyverse)
library(lubridate)
library(ggmap)
library(ggrepel)
library(geosphere)
library(patchwork)
library(contoureR)
library(ade4)
library(tidygraph)
library(ggraph)
library(igraph)
library(ggtree)
library(TipDatingBeast)

####################
### Read in data ###
####################

######################
## Read in metadata ##
######################

EthiopiaMetadata <- read_csv("Ethiopia_all_metadata.csv") 

#############################################
## Convert original dates to decimal dates ##
#############################################

EthiopiaMetadata <- EthiopiaMetadata %>% 
  mutate(decimalDate = decimal_date(dmy(Date)))

##################################################
## Save decimal dates to add to BEAST alignment ##
##################################################

EthiopiaAllDates <- EthiopiaMetadata %>%
  mutate(new_id = paste(lane_id, 
                        decimalDate, 
                        sep = "_")) %>% 
  select(lane_id, 
         new_id)

write_delim(EthiopiaAllDates, 
            "Ethiopia_new_dates.txt", 
            delim = "\t", 
            col_names = FALSE)

#########################################
## Define study lat/long for stamenmap ##
#########################################

EthiopiaMap <- c(left = 33, 
                 bottom = 3, 
                 right = 48, 
                 top = 15)

##############################################
## Download stamenmap for above coordinates ##
##############################################

EthiopiaStamenMap <- get_stamenmap(EthiopiaMap, 
                                   zoom = 7, 
                                   maptype = "terrain", 
                                   color = "bw", 
                                   force = TRUE)

####################################################
## Plot all points and colour by host (Figure 1A) ##
####################################################

EthiopiaStamenMap_allHost <- ggmap(EthiopiaStamenMap) + 
  geom_point(data = EthiopiaMetadata, 
             aes(x = Longitude, 
                 y = Latitude, 
                 colour = Host, 
                 shape = Host), 
             alpha = 1, 
             fill = NA, 
             size = 3) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  labs(colour = 'Host')

#########################################################
## Extract single example of herd locations for labels ##
#########################################################

EthiopiaHerdCoords <- EthiopiaMetadata %>% 
  group_by(Herd_code) %>% 
  filter(row_number()==1)

####################################################
## Create a map with herd locations plotted on it ##
## (Figure 3A)                                    ##
####################################################

EthiopiaStamenMapHerds <- ggmap(EthiopiaStamenMap) +
  geom_point(data = EthiopiaHerdCoords, 
             aes(x = Longitude, 
                 y = Latitude), 
             alpha = 1, 
             fill = NA, 
             size = 2) +
  geom_label_repel(data = EthiopiaHerdCoords, 
                   aes(x = Longitude, 
                       y = Latitude, 
                       label = Herd_code)) + 
  xlab("Longitude") + 
  ylab("Latitude")

####################################################
## Download stamen map with zoomed in coordinates ##
####################################################

EthiopiaMapZoom <- c(left = 38.4, 
                     bottom = 8.7, 
                     right = 39.2, 
                     top = 9.3)

EthiopiaStamenMapZoom <- get_stamenmap(EthiopiaMapZoom, 
                                       zoom = 10, 
                                       maptype = "terrain", 
                                       color = "bw", 
                                       force = TRUE)

###########################################################
## Create a zoomed map with herd locations plotted on it ##
## (Figure 3A insert)                                    ##
###########################################################

EthiopiaStamenMapHerdsZoom <- ggmap(EthiopiaStamenMapZoom) +
  geom_point(data = EthiopiaHerdCoords, 
             aes(x = Longitude, 
                 y = Latitude), 
             alpha = 1, 
             fill = NA, 
             size = 2) +
  geom_label_repel(data = EthiopiaHerdCoords, 
                   aes(x = Longitude, 
                       y = Latitude, 
                       label = Herd_code)) + 
  xlab("Longitude") + 
  ylab("Latitude")

############################################################
## Calculate and plot convex hull for each clonal complex ##
############################################################

Ethiopiasubset <- NULL

Ethiopiahull <- NULL

Ethiopiapolygon <- NULL

for(cc in unique(EthiopiaMetadata$CC)){
  
  # Get the subset of the data for the current cc
  Ethiopiasubset[[cc]] <- EthiopiaMetadata[EthiopiaMetadata$CC == cc, ]
  
  # Calculate the convex hull
  Ethiopiahull[[cc]] <- getConvexHull(x=Ethiopiasubset[[cc]]$Longitude, 
                                      y=Ethiopiasubset[[cc]]$Latitude)
  
  # Calculate polygon coordinates
  Ethiopiapolygon[[cc]] <- cbind(Ethiopiasubset[[cc]][Ethiopiahull[[cc]], "Longitude"], 
                                 Ethiopiasubset[[cc]][Ethiopiahull[[cc]], "Latitude"])
  
}

#################################################
## Collapse polygon list into single dataframe ##
#################################################

EthiopiaPolygonAll <- bind_rows(Ethiopiapolygon, 
                                .id = "CC")

#####################################################
## Overlay polygons on coordinate plot (Figure 2A) ##
#####################################################

EthiopiaCoordPlotPolygon <- ggplot() +
  geom_point(data = EthiopiaMetadata, 
             aes(x = Longitude, 
                 y = Latitude, 
                 color = CC)) +
  geom_polygon(data = EthiopiaPolygonAll, 
               aes(x = Longitude, 
                   y = Latitude, 
                   group = CC, 
                   fill = CC), 
               alpha = 0.3) +
  theme_bw() +
  labs(color = "Clonal complex", 
       fill = "Clonal complex") + 
  theme(legend.position = "none")

#################################
## Read in pairsnp output file ##
#################################

EthiopiaPWsnps <- read_csv("all_Ethiopia_300320_masked_noref_snps.csv")

###########################################
## Convert pairsnp output into dataframe ##
###########################################

EthiopiaPWsnpsDecon <- data.frame( t(combn(names(EthiopiaPWsnps),2)), 
                                   dist=t(EthiopiaPWsnps)[lower.tri(EthiopiaPWsnps)] )
colnames(EthiopiaPWsnpsDecon) <- c("Taxon1", 
                                   "Taxon2", 
                                   "dist")
EthiopiaPWsnpsDecon$Taxon1 <- gsub('# ', 
                                   "",
                                   EthiopiaPWsnpsDecon$Taxon1)

##########################################
### Geographical localization analyses ###
##########################################

################################################################
## Create a copy of deconvoluted pairwise SNP distance matrix ##
################################################################

EthiopiaSNPsGeoDists <- EthiopiaPWsnpsDecon

###########################################################################
## Create a new dataframe with lane_id, Spoligotype, Longitude, Latitude ##
###########################################################################

EthiopiaMetadataGeoDists <- EthiopiaMetadata %>% 
  select(lane_id, 
         Region, 
         Spoligotype, 
         CC, 
         Host, 
         Longitude, 
         Latitude, 
         Animal_ID, 
         Herd_code)

#################################
## Add metadata for each taxon ##
#################################

EthiopiaSNPsGeoDists <- EthiopiaSNPsGeoDists %>% 
  left_join(EthiopiaMetadataGeoDists, 
            by = c("Taxon1" = "lane_id")) %>% 
  left_join(EthiopiaMetadataGeoDists, 
            by = c("Taxon2" = "lane_id"))

######################################################
## Calculate pw geographical distance for each pair ##
######################################################

EthiopiaSNPsGeoDists <- EthiopiaSNPsGeoDists %>% 
  mutate(geo_dist = distHaversine(cbind(Longitude.x, 
                                        Latitude.x), 
                                  cbind(Longitude.y, 
                                        Latitude.y))/1000)

#################################################################
## Create new dataframe for geographical localization analyses ##
#################################################################

EthiopiaSNPsGeoDistsFiltered <- EthiopiaSNPsGeoDists %>%
  select(Taxon1, 
         Taxon2, 
         CC.x, 
         CC.y, 
         Host.x, 
         Host.y, 
         Animal_ID.x, 
         Animal_ID.y, 
         Herd_code.x, 
         Herd_code.y, 
         dist, 
         geo_dist) %>% 
  mutate(CC_match = ifelse(CC.x == CC.y, 
                           paste("Within clonal complex"), 
                           paste("Between clonal complex"))) %>% 
  mutate(host_match = ifelse(Host.x == Host.y, 
                             paste("Same host"), 
                             paste("Different host"))) %>% 
  mutate(animal_match = ifelse(Animal_ID.x == Animal_ID.y, 
                               paste("Same animal"), 
                               paste("Different animal")))

#############################################################
## Plot pairwise SNP distance against geographic distance. ##
## Facet by clonal complex (Figure 2B)                     ##
#############################################################

EthiopiaDistancePlotSameCC <- ggplot(data = EthiopiaSNPsGeoDistsFilteredSameCC) + 
  geom_point(aes(x = dist, 
                 y = geo_dist)) +
  xlab("Pairwise SNP distance") +
  ylab("Geographic distance (kilometers)") +
  theme_bw() +
  facet_wrap(~CC.x, 
             ncol = 3)

##########################################################
## Plot histogram of pairwise SNP distances (Figure 2C) ##
##########################################################

EthiopiaSNPHistogram <- ggplot(data = EthiopiaPWsnpsDecon, 
                               aes(x = dist)) + 
  geom_histogram(bins = 200, 
                 alpha = 0.8, 
                 position = "identity") + 
  theme_bw()  + 
  xlab("Pairwise SNP distance") + 
  ylab("Frequency") 

##############################################################
## Create box plot of pairwise distances per CC (Figure 2D) ##
##############################################################

EthiopiaSNPsGeoDistsFilteredSameCC <- EthiopiaSNPsGeoDistsFiltered %>% 
  filter(CC_match == "Within clonal complex")

EthiopiaSNPsCCBox <- ggplot(data = EthiopiaSNPsGeoDistsFilteredSameCC, 
                              aes(x = CC.x, 
                                  y = dist, 
                                  fill = CC.x)) + 
  geom_jitter(color="black", 
              size=0.4, 
              alpha=0.6) +
  geom_boxplot(alpha = 0.6) +
  theme_bw() + 
  xlab("Clonal complex") + 
  ylab("Pairwise SNP distance") + 
  theme(legend.position = "none")

##############################################################################
## Create box plot of pairwise distances divided by within and between host ##
## (Figure 2E)                                                              ##
##############################################################################

EthiopiaSNPsAnimalBox <- ggplot(data = EthiopiaSNPsGeoDistsFiltered, 
                                aes(x = animal_match, 
                                    y = dist, 
                                    fill = animal_match)) + 
  geom_jitter(color="black", size=0.4, alpha=0.6) +
  geom_boxplot(alpha = 0.6) +
  theme_bw() + 
  xlab("") + 
  ylab("Pairwise SNP distance") + 
  theme(legend.position = "none")

###########################################################
## Create new dataframe to examine within animal pw SNPs ##
###########################################################

EthiopiaSNPsGeoDistsFilteredAnimal <- EthiopiaSNPsGeoDistsFiltered %>%
  filter(animal_match != "Different animal") %>%
  filter(Host.x != "Dromedary" & Host.y != "Dromedary")

#############################################################################
## Plot histogram of pairwise SNP distances for within animals (Figure 2F) ##
#############################################################################

EthiopiaSNPHistogramSameAnimal <- ggplot(data = EthiopiaSNPsGeoDistsFilteredAnimal, 
                                         aes(x = dist)) + 
  geom_histogram(bins = 200, 
                 alpha = 0.8, 
                 position = "identity") + 
  theme_bw()  + 
  xlab("Pairwise SNP distance") + 
  ylab("Frequency")

########################################################################
## Zoom in on above histogram to look at peaks between 0 and 100 SNPs ##
########################################################################

EthiopiaSNPHistogramSameAnimalZoom <- EthiopiaSNPHistogramSameAnimal + 
  coord_cartesian(xlim = c(0,
                           100))

############################################
## Create patchwork of SNP/geo dist plots ##
############################################

EthiopiaFig2 <- (EthiopiaCoordPlotPolygon|EthiopiaDistancePlotSameCC)/
  (EthiopiaSNPHistogram|EthiopiaSNPsCCBox)/
  (EthiopiaSNPsAnimalBox|EthiopiaSNPHistogramSameAnimal)

##########################################
## Identify and remove abbatoir samples ##
##########################################

EthiopiaSNPsGeoDistsFilteredSameCCnoUnknown8 <- EthiopiaSNPsGeoDistsFilteredSameCC %>% 
  filter(CC.x != "Unknown8") %>% 
  filter(!Herd_code.x %in% c("AA","BA","FA", "GA", "NA", "SA")) %>% 
  filter(!Herd_code.y %in% c("AA","BA","FA", "GA", "NA", "SA"))

#########################################
## Perform Mantel tests on Af2 and Eu3 ##
#########################################

EthiopiaCCsubset <- NULL

EthiopiaCCSNPdists <- NULL

EthiopiaCCGeoDists <- NULL

EthiopiaCCMantel <- NULL

for(cluster in unique(EthiopiaSNPsGeoDistsFilteredSameCCnoUnknown8$CC.x)){
  
  #################################################
  # Get the subset of the data for the current CC #
  #################################################
  
  EthiopiaCCsubset[[cluster]] <- EthiopiaSNPsGeoDistsFilteredSameCC[EthiopiaSNPsGeoDistsFilteredSameCCnoUnknown8$CC.x == cluster, ]
  
  ####################################################
  # Create a matrix of the SNP distances for each CC #
  ####################################################
  
  EthiopiaCCSNPdists[[cluster]] <- dist(EthiopiaCCsubset[[cluster]]$dist)
  
  ############################################################
  # Create a matrix of the geographic distances for each CC ##
  ############################################################
  
  EthiopiaCCGeoDists[[cluster]] <- dist(EthiopiaCCsubset[[cluster]]$geo_dist)
  
  ##################################
  # Conduct Mantel test on each CC #
  ##################################
  
  EthiopiaCCMantel[[cluster]] <- mantel.rtest(EthiopiaCCSNPdists[[cluster]], 
                                              EthiopiaCCGeoDists[[cluster]], 
                                              nrepet = 999)
}

################################################
## Create new df containing nodes for network ##
################################################

EthiopiaNetworkMetadata <- EthiopiaMetadata

EthiopiaNetworkMetadata_nodes <- EthiopiaNetworkMetadata %>%
  mutate(id = row_number()) %>% 
  select(id, 
         lane_id, 
         Longitude, 
         Latitude, 
         Host, 
         Animal_ID, 
         Herd_code)

##############################################################
## Extract pairwise SNP comparisons < 16 (edges of network) ##
##############################################################

EthiopiaPWsnpsDecon_edges_15_snp <- EthiopiaPWsnpsDecon %>%
  filter(dist < 16) %>% 
  select(Taxon1, 
         Taxon2, 
         dist)

EthiopiaPWsnpsDecon_edges_15_snp <- EthiopiaPWsnpsDecon_edges_15_snp %>%
  left_join(EthiopiaNetworkMetadata_nodes[,c(1,2)], 
            by = c("Taxon1" = "lane_id")) %>% 
  rename(from.id = id)

EthiopiaPWsnpsDecon_edges_15_snp <- EthiopiaPWsnpsDecon_edges_15_snp %>% 
  left_join(EthiopiaNetworkMetadata_nodes[,c(1,2)], 
            by = c("Taxon2" = "lane_id")) %>% 
  select(from.id, 
         to.id = id, 
         dist)

###########################
## Create network object ##
###########################

EthiopiaNetworkMetadata_nodes_routes <- tbl_graph(nodes = EthiopiaNetworkMetadata_nodes, 
                                                  edges = EthiopiaPWsnpsDecon_edges_15_snp, 
                                                  directed = TRUE)

EthiopiaNetworkMetadata_nodes_routes_components <- components(EthiopiaNetworkMetadata_nodes_routes)

EthiopiaNetworkMetadata_nodes_networked <- EthiopiaNetworkMetadata_nodes %>% 
  mutate(network_id = EthiopiaNetworkMetadata_nodes_routes_components$membership)

###################################
## Filter out network singletons ##
###################################

EthiopiaNetworkMetadata_nodes_networkedFiltered <- EthiopiaNetworkMetadata_nodes_networked %>% 
  group_by(network_id) %>%
  filter(n() > 1) %>%
  ungroup() %>% 
  mutate(id = row_number()) %>% 
  select(id, 
         lane_id, 
         Longitude, 
         Latitude, 
         Host, 
         Animal_ID, 
         Herd_code)

##############################################################
## Extract pairwise SNP comparisons < 16 (edges of network) ##
##############################################################

EthiopiaPWsnpsDecon_edges_15_snpFiltered <- EthiopiaPWsnpsDecon %>%
  filter(dist < 16) %>% 
  select(Taxon1, 
         Taxon2, 
         dist)

EthiopiaPWsnpsDecon_edges_15_snpFiltered <- EthiopiaPWsnpsDecon_edges_15_snpFiltered %>%
  left_join(EthiopiaNetworkMetadata_nodes_networkedFiltered[,c(1,2,6)], 
            by = c("Taxon1" = "lane_id")) %>% 
  rename(from.id = id)

EthiopiaPWsnpsDecon_edges_15_snpFiltered <- EthiopiaPWsnpsDecon_edges_15_snpFiltered %>% 
  left_join(EthiopiaNetworkMetadata_nodes_networkedFiltered[,c(1,2,6)], 
            by = c("Taxon2" = "lane_id")) %>%
  mutate(animal_match = ifelse(Animal_ID.x == Animal_ID.y, 
                               paste("Same animal"), 
                               paste("Different animal"))) %>% 
  select(from.id, 
         to.id = id, 
         animal_match, 
         dist)

###############################
## Create new network object ##
###############################

EthiopiaNetworkMetadata_nodes_routesFiltered <- tbl_graph(nodes = EthiopiaNetworkMetadata_nodes_networkedFiltered, 
                                                          edges = EthiopiaPWsnpsDecon_edges_15_snpFiltered, 
                                                          directed = TRUE)

###############################################################
## Plot network without isolates not in networks (Figure 3B) ##
###############################################################

Ethiopia_15_snp_networkHostFiltered <- ggraph(EthiopiaNetworkMetadata_nodes_routesFiltered, 
                                              layout = "nicely") + 
  geom_edge_link(aes(edge_colour = animal_match, 
                     label = dist), 
                 edge_width = 1) +
  geom_edge_link(aes(filter = dist < 6), 
                 edge_colour = 'black', 
                 edge_width = 1, 
                 edge_linetype = "dashed") + 
  geom_node_point(aes(colour = Animal_ID, 
                      shape = Host), 
                  size = 6) +
  geom_node_text(aes(label = Herd_code), 
                 size = 4) +
  theme_graph() + 
  labs(shape = "Host")+ 
  theme(legend.position = "none")

#####################################################################
### Select random samples for EU1 and EU2 for global context tree ###
#####################################################################

##########################################
## Read in list of EU1 and EU2 isolates ##
##########################################

EthEU1EU2 <- read_csv("Eth_EU1_EU2.csv")

################################
## Create new dfs for each CC ##
################################

EthEU1 <- EthEU1EU2 %>% 
  filter(Loiseau_CC == "EU1")

EthEU2 <- EthEU1EU2 %>% 
  filter(Loiseau_CC == "EU2")

##############################################
## Randomly sample 100 isolates for each CC ##
##############################################

EthEU1Random <- sample_n(EthEU1,
                         100)

EthEU2Random <- sample_n(EthEU2,
                         100)

################################
## Write out sampled isolates ##
################################

write_csv(EthEU1Random, "EthEU1Random.csv")

write_csv(EthEU2Random, "EthEU2Random.csv")

####################
### Eu3 analyses ###
####################

############################
## Read in BEAST MCC tree ##
############################

EU3BEASTtree <- read.beast("EU3_no_outgroup_030420_masked_snps_new_dated_strict_constant_MCC.tree")

###########################################
## Read in Eu3 metadata for all isolates ##
###########################################

EU3Metadata <- read_csv("EU3_metadata.csv")

#####################################################
## Remove isolate dates from BEAST tree tip labels ##
#####################################################

EU3BEASTtree@phylo$tip.label <- EU3Metadata$Isolate_ID[match(EU3BEASTtree@phylo$tip.label,
                                                             EU3Metadata$tip_label)]

########################################################
## Plot BEAST tree with vertical lines for date range ##
########################################################

EU3BEASTggtree <- ggtree(EU3BEASTtree, 
                         right=TRUE, 
                         mrsd="2018-04-23") + 
  geom_tiplab(size = 2, 
              align = TRUE, 
              offset = 5) +
  geom_range("height_0.95_HPD", 
             color='red', 
             size=1, 
             alpha=.5, 
             branch.length="height") + 
  theme_tree2() +
  scale_x_continuous(breaks=seq(1570, 
                                2020, 
                                50), 
                     minor_breaks=seq(1970, 
                                      2020, 
                                      10)) + 
  theme(panel.grid.major = element_line(color="black", 
                                        size=.2),
        panel.grid.minor   = element_line(color="grey", 
                                          size=.2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) + 
  xlim_tree(2020)

############################################
## Add metadata to BEAST tree (Figure 4B) ##
############################################

EU3MetadataDF <- as.data.frame(EU3Metadata[,3:5])

rownames(EU3MetadataDF) = EU3Metadata$Isolate_ID

EU3BEASTggtreeMeta2 <- gheatmap(EU3BEASTggtree, 
                                EU3MetadataDF, 
                                colnames_position = "top", 
                                font.size=3, 
                                offset=25, 
                                width=0.15, 
                                colnames_angle = 90, 
                                hjust = 0) 

##############################
## Convert BEAST tree to df ##
##############################

EU3BEASTtreeTibble <- as_tibble(EU3BEASTtree)

##########################
### BEAST DTR analysis ###
##########################

#####################################################################
## Read in BEAST xml file and create 20 replicates with randomized ##
## tip dates                                                       ##
#####################################################################

RandomDates(name = "EU3_no_outgroup_030420_masked_snps_new_dated_strict_constant_1", 
            reps = 20)

#######################################################
## Calculate and plot DRT results for each replicate ##
#######################################################

PlotDRT(name = "EU3_no_outgroup_030420_masked_snps_new_dated_strict_constant_1", 
        reps = 20, 
        burnin = 0.1)

#########################
## Read in clock rates ##
#########################

DTRfile <- read_csv("EU3_clock.rate.stats.csv")

#################################################
## Prepare colours for observed and replicates ##
#################################################

DTRdfColour <- DTRfile %>%
  mutate(colour = ifelse(calibr == 0, paste("red"), 
                         paste("black")))

DTRdfLines <- DTRdfColour[DTRdfColour$calibr == 0,]

#############################################
## Plot BEAST DTR results (Supp. Figure 1) ##
#############################################

DTRplot <- ggplot(data = DTRdfColour) + 
  geom_point(aes(x = factor(calibr), 
                 y = median, 
                 color = colour)) + 
  geom_errorbar(aes(x = factor(calibr), 
                    ymin = lowerHPD, 
                    ymax = HigherHPD, 
                    color = colour), 
                width = 0.2, 
                position = position_dodge(0.9)) + 
  theme_bw() + 
  ylab("Clock rate") + 
  xlab("Replicate") +
  geom_hline(data = DTRdfLines, 
             aes(yintercept = lowerHPD), 
             colour = "red", 
             linetype = "dashed") +
  geom_hline(data = DTRdfLines, 
             aes(yintercept = HigherHPD), 
             colour = "red", 
             linetype = "dashed") +
  scale_color_manual(values = c("black", 
                                "red")) + 
  theme(legend.position = "none")

