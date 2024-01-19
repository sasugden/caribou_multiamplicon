###### PREPARE WORKSPACE FOR STATISTICAL ANALYSIS ####
## [1] Prepare sample data ####
sample_data <- tibble::add_column(sample_data, .before=1, SampleID = rownames(sample_data))
sample_data$Group <- factor(sample_data$Group, levels=c("deep-snow", "shallow-snow",
                                                        "Revelstoke pen", "LARS"))

# Create a separate metadata sheet for each amplicon or amplicon subset.
# This will create space to store data on e.g. species richness, Shannon diversity.
sample_data.mb <- list(sample_data, sample_data, sample_data, sample_data, 
                       sample_data, sample_data, sample_data)
names(sample_data.mb) <- c("X16S", "X18S", "FITS", "PITS", "protists", "algae", "lichens")

sample_data.mb[[1]]$Amplicon <- rep("X16S")
sample_data.mb[[2]]$Amplicon <- rep("X18S")
sample_data.mb[[3]]$Amplicon <- rep("FITS")
sample_data.mb[[4]]$Amplicon <- rep("PITS")
sample_data.mb[["protists"]]$Amplicon <- "protists"
sample_data.mb[["algae"]]$Amplicon <- "algae"
sample_data.mb[["lichens"]]$Amplicon <- "lichens"

# Define theme for plots.
plot_theme <- theme(plot.tag = element_text(color="black", face="bold"),
                    strip.text = element_text(size=9),
                    panel.grid = element_blank(),
                    panel.border = element_rect(color="black", fill=NA, size=1),
                    axis.text = element_text(color= "black", size=8),
                    axis.title = element_text(color="black", size=9),
                    legend.title = element_text(color="black", size=9),
                    plot.title=element_text(color="black", size=10),
                    legend.text = element_text(color="black", size=8))

# Name the different phyloseq objects.
names(pseq.clr) <- c("X16S", "X18S", "FITS", "PITS")
names(pseq.rarefied) <- c("X16S", "X18S", "FITS", "PITS")

df <- as.data.frame(cbind(sample_data(pseq.clr[[4]])))
df <- tibble::add_column(df, .before=1, SampleID = rownames(df))
sample_data(pseq.clr[[4]]) <- df
sample_data(pseq.rarefied[[4]]) <- df
rm(df)

# Remove the designations "uncultured" from the taxonomy table and fix the protist taxonomy
temp <- as.data.frame(cbind(tax_table(pseq.clr[[2]])))
temp[temp=="uncultured"] <- NA
tax_table(pseq.clr[[2]]) <- as.matrix(temp)

temp <- as.data.frame(cbind(tax_table(pseq.rarefied[[2]])))
temp[temp=="uncultured"] <- NA
tax_table(pseq.rarefied[[2]]) <- as.matrix(temp)

rm(temp)

## [2] Replace 'NAs' in the taxonomy table with higher-level taxon ranks. ####
# For example, if an ASV is classified as 'Lachnospiraceae' (family) but lacks a genus classification,
# this fills in the genus with 'Uncl_Lachnospiraceae.'
# This avoids later errors with the 'tax_glom' function in phyloseq.
pseq.rarefied.sub <- pseq.rarefied

for(p in c(1:4)){
  if(p %in% c(1,3,4)){k==6}
  if(p==2){k==7}
  
  for(i in 1:nrow(tax_table(pseq.rarefied.sub[[p]]))){
    for(j in 2:k){
      if(is.na(tax_table(pseq.rarefied.sub[[p]])[i,j])==TRUE){
        if(substr(tax_table(pseq.rarefied.sub[[p]])[i,j-1], 1, 4)=="Uncl"){
          tax_table(pseq.rarefied.sub[[p]])[i,j] <- tax_table(pseq.rarefied.sub[[p]])[i,j-1]}
        else {
          tax_table(pseq.rarefied.sub[[p]])[i,j] <- paste0("Uncl_", tax_table(pseq.rarefied.sub[[p]])[i,j-1])}}
    }}
}

pseq.clr.sub <- pseq.clr

for(p in c(1:4)){
  if(p %in% c(1,3,4)){k==6}
  if(p==2){k==7}
  
  for(i in 1:nrow(tax_table(pseq.clr.sub[[p]]))){
    for(j in 2:k){
      if(is.na(tax_table(pseq.clr.sub[[p]])[i,j])==TRUE){
        if(substr(tax_table(pseq.clr.sub[[p]])[i,j-1], 1, 4)=="Uncl"){
          tax_table(pseq.clr.sub[[p]])[i,j] <- tax_table(pseq.clr.sub[[p]])[i,j-1]}
        else {
          tax_table(pseq.clr.sub[[p]])[i,j] <- paste0("Uncl_", tax_table(pseq.clr.sub[[p]])[i,j-1])}}
    }}
}

## [3] Create separate phyloseq objects for protists, algae, and lichens. ####
# For protists and algae, use the built-in classifications from the taxonomy table.
pseq.clr$protists <- subset_taxa(pseq.clr[[2]], Kingdom=="Protista")
pseq.clr$algae <-  subset_taxa(pseq.clr[[2]], Kingdom=="Algae")

pseq.rarefied$protists <-  subset_taxa(pseq.rarefied[[2]], Kingdom=="Protista")
pseq.rarefied$algae <- subset_taxa(pseq.rarefied[[2]], Kingdom=="Algae")

pseq.clr.sub$protists <- subset_taxa(pseq.clr.sub[[2]], Kingdom=="Protista")
pseq.clr.sub$algae <- subset_taxa(pseq.clr.sub[[2]], Kingdom=="Algae")

pseq.rarefied.sub$protists <- subset_taxa(pseq.rarefied.sub[[2]], Kingdom=="Protista")
pseq.rarefied.sub$algae <- subset_taxa(pseq.rarefied.sub[[2]], Kingdom=="Algae")

# For lichens, create lists of lichenized taxa.
# This is the list of lichenized fungi curated by Lucking et al. (2016)
lichenized.fungi <- read.csv("~/caribou/github_publish/lichenized_fungi.csv")

# Subset FUNGuild identifications to only lichenized fungi.
lichenized.guild <- subset(fungal.guilds, guild=="Lichenized")

pseq.clr$lichens <- subset_taxa(pseq.clr[[3]], Genus %in% lichenized.fungi$Genus | 
                                  taxa_names(pseq.clr[[3]]) %in% rownames(lichenized.guild))
pseq.rarefied$lichens <- subset_taxa(pseq.rarefied[[3]], Genus %in% lichenized.fungi$Genus | 
                                       taxa_names(pseq.rarefied[[3]]) %in% rownames(lichenized.guild))
pseq.clr.sub$lichens <- subset_taxa(pseq.clr.sub[[3]], Genus %in% lichenized.fungi$Genus | 
                                      taxa_names(pseq.clr.sub[[3]]) %in% rownames(lichenized.guild))
pseq.rarefied.sub$lichens <- subset_taxa(pseq.rarefied.sub[[3]], Genus %in% lichenized.fungi$Genus | 
                                           taxa_names(pseq.rarefied.sub[[3]]) %in% rownames(lichenized.guild))

# Remove Sarea, which is not a true lichen.
pseq.clr$lichens <- subset_taxa(pseq.clr$lichens, Genus !="Sarea")
pseq.rarefied$lichens <- subset_taxa(pseq.rarefied$lichens, Genus !="Sarea")
pseq.clr.sub$lichens <- subset_taxa(pseq.clr.sub$lichens, Genus !="Sarea")
pseq.rarefied.sub$lichens <- subset_taxa(pseq.rarefied.sub$lichens, Genus !="Sarea")

# Remove samples with no reads from these taxon subsets.
pseq.clr$protists <- subset_samples(pseq.clr$protists, sample_sums(pseq.clr$protists) > 0)
pseq.clr.sub$protists <- subset_samples(pseq.clr.sub$protists, sample_sums(pseq.clr.sub$protists) > 0)
pseq.rarefied$protists <- subset_samples(pseq.rarefied$protists, sample_sums(pseq.rarefied$protists) > 0)
pseq.rarefied.sub$protists <- subset_samples(pseq.rarefied.sub$protists, sample_sums(pseq.rarefied.sub$protists) > 0)

pseq.clr$algae <- subset_samples(pseq.clr$algae, sample_sums(pseq.clr$algae) > 0)
pseq.clr.sub$algae <- subset_samples(pseq.clr.sub$algae, sample_sums(pseq.clr.sub$algae) > 0)
pseq.rarefied$algae <- subset_samples(pseq.rarefied$algae, sample_sums(pseq.rarefied$algae) > 0)
pseq.rarefied.sub$algae <- subset_samples(pseq.rarefied.sub$algae, sample_sums(pseq.rarefied.sub$algae) > 0)

pseq.clr$lichens <- subset_samples(pseq.clr$lichens, sample_sums(pseq.clr$lichens) > 0)
pseq.clr.sub$lichens <- subset_samples(pseq.clr.sub$lichens, sample_sums(pseq.clr.sub$lichens) > 0)
pseq.rarefied$lichens <- subset_samples(pseq.rarefied$lichens, sample_sums(pseq.rarefied$lichens) > 0)
pseq.rarefied.sub$lichens <- subset_samples(pseq.rarefied.sub$lichens, sample_sums(pseq.rarefied.sub$lichens) > 0)

## [4] Create a summary of all taxon abundances and prevalences, for subsetting purposes later ####
prev.data <- list()

for(p in c(1:7)){
  prev.data[[p]] <- list()
  
  if(p %in% c(1,3,7)){taxranks <- c("Phylum", "Class", "Order", "Family", "Genus")}
  if(p == 4){taxranks <- c("Class", "Order", "Family", "Genus")}
  if(p %in% c(2,5,6)){taxranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")}
  
  for(i in c(1:(length(taxranks)+1))){
    temp.physeq <- pseq.clr.sub[[p]]
    
    # Agglomerate to the taxonomic rank of interest.
    if(i < (length(taxranks)+1)){
      temp.physeq <- tax_glom(temp.physeq, taxrank=taxranks[[i]])
      taxa_names(temp.physeq) <- as.data.frame(cbind(tax_table(temp.physeq)))[,taxranks[[i]]]
    }
    
    # Calculate prevalence (number of samples in which each taxon occurs)
    prev <- as.data.frame(microbiome::prevalence(temp.physeq, detection=0, count=TRUE))
    colnames(prev) <- "Prevalence"
    prev$Taxon <- rownames(prev)
    
    # Convert taxon abundances to relative abundances.
    temp.physeq <- transform_sample_counts(temp.physeq, function(x) 100*x/sum(x))
    temp.otu <- as.data.frame(otu_table(temp.physeq))
    
    # Save the mean taxon relative abundance across all samples.
    prev$Abundance <- colMeans(temp.otu)
    
    # Add population and herd information to the OTU table.
    temp.otu <- transform(merge(sample_data(pseq.clr.sub[[p]])[,c("Group", "Herd")], temp.otu, 
                    by=0, all.x=TRUE, all.y=FALSE),
              row.names=Row.names, Row.names=NULL)
    
    # Calculate mean relative abundance per population.
    temp.grp <- temp.otu %>% 
      dplyr::group_by(Group) %>% 
      dplyr::summarise_if(is.numeric, mean) %>%
      tibble::column_to_rownames(., "Group") %>%
      t() %>%
      as.data.frame()
    
    # Calculate mean relative abundance per herd.
    temp.hrd <- temp.otu %>%  
      dplyr::group_by(Herd) %>%
      dplyr::summarise_if(is.numeric, mean) %>%
      tibble::column_to_rownames(., "Herd") %>%
      t() %>%
      as.data.frame() %>%
      dplyr::select(., -`Revelstoke pen`, -`LARS`)
    
    # Save final data.
    temp.abnd <- as.data.frame(dplyr::bind_cols(temp.grp, temp.hrd))
    prev <- cbind(prev, temp.abnd)
    prev.data[[p]][[i]] <- prev
    
    # Clean workspace.
    rm(prev, temp.physeq, temp.otu, temp.grp, temp.hrd, temp.abnd)
  }
  names(prev.data[[p]]) <- c(taxranks, "ASV")
}

names(prev.data) <- c("X16S", "X18S", "FITS", "PITS", "protists", "algae", "lichens")

####### INITIAL CALCULATIONS ####
## [1] Calculate rarefaction curves. ####
# Define rarefaction curve calculation function (obtained from https://github.com/joey711/phyloseq/issues/143)
calculate_rarefaction_curves <- function(psdata, measures, depths) {
  require('plyr') # ldply
  require('reshape2') # melt
  
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}

# Calculate rarefaction curves.
# (Note that this is only done for the four sequenced)
rarefaction_curves <- list()

for(p in c(1:4)){
  rarefaction_curve_data <- calculate_rarefaction_curves(pseq.clr[[p]],
                                                         c('Observed', 'Shannon'),
                                                         rep(c(1, 10, 100, seq(1000, 75000, by=1000)), each = 10))
  
  # Summarize alpha diversity and add sample data.
  rarefaction_curve_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'),
                                     summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))
  rarefaction_curve_summary$Sample <- gsub("X1", "1", rarefaction_curve_summary$Sample)
  rarefaction_curve_summary$Sample <- gsub("X2", "2", rarefaction_curve_summary$Sample)
  rarefaction_curve_summary$Sample <- gsub(".", "-", rarefaction_curve_summary$Sample, fixed=TRUE)
  rarefaction_curve_summary_verbose <- merge(rarefaction_curve_summary, 
                                             data.frame(sample_data(pseq.clr[[p]])),
                                             by.x = 'Sample', by.y = 'row.names')
  
  rarefaction_curves[[p]] <- rarefaction_curve_summary_verbose
  
  rm(rarefaction_curve_summary_verbose, rarefaction_curve_summary, rarefaction_curve_data)
}

ggplot(subset(rarefaction_curves[[3]], Measure=="Observed"), aes(x=Depth, y=Alpha_diversity_mean)) + 
  geom_line(aes(group=Sample)) + 
  facet_wrap(~Group)

## [2] Calculate richness and diversity using rarefaction and extrapolation from unrarefied data ####
# Run iNEXT on each phyloseq object.
iNext.outputs.3 <- list(
  X16S = iNEXT::iNEXT(as.data.frame(t(otu_table(pseq.clr.sub[[1]]))), q=0, datatype="abundance"),
  X18S = iNEXT::iNEXT(as.data.frame(t(otu_table(pseq.clr.sub[[2]]))), q=0, datatype="abundance"),
  FITS = iNEXT::iNEXT(as.data.frame(t(otu_table(pseq.clr.sub[[3]]))), q=0, datatype="abundance"),
  PITS = iNEXT::iNEXT(as.data.frame(t(otu_table(pseq.clr.sub[[4]]))), q=0, datatype="abundance"),
  protists = iNEXT::iNEXT(as.data.frame(t(otu_table(pseq.clr.sub[[5]]))), q=0, datatype="abundance"),
  algae = iNEXT::iNEXT(as.data.frame(t(otu_table(pseq.clr.sub[[6]]))), q=0, datatype="abundance"),
  lichens = iNEXT::iNEXT(as.data.frame(t(otu_table(pseq.clr.sub[[7]]))), q=0, datatype="abundance")
)

# Extract richness and diversity estimates and save them as part of the sample data for each data set.
for(p in c(1:7)){
  temp <- iNext.outputs.3[[p]]$AsyEst
  
  temp.richness <- subset(temp, Diversity=="Species richness")
  temp.richness$SampleID <- temp.richness$Assemblage
  temp.richness$Observed_Extrap <- temp.richness$Estimator
  temp.richness$Observed_Extrap <- round(temp.richness$Observed_Extrap, digits = 0)
  
  temp.diversity <- subset(temp, Diversity=="Shannon diversity")
  temp.diversity$SampleID <- temp.diversity$Assemblage
  temp.diversity$Shannon_Extrap <- log(temp.diversity$Estimator)
  
  pseq.est <- estimate_richness(pseq.rarefied[[p]], measures=c("Observed", "Chao1", "Shannon"))
  pseq.est$SampleID <- sample_data(pseq.rarefied[[p]])$SampleID
  
  sample_data.mb[[p]] <- merge(sample_data.mb[[p]], pseq.est,
                               by="SampleID", all=TRUE)
  sample_data.mb[[p]] <- merge(sample_data.mb[[p]], temp.richness[,c("SampleID", "Observed_Extrap")],
                               by="SampleID", all=TRUE)
  sample_data.mb[[p]] <- merge(sample_data.mb[[p]], temp.diversity[,c("SampleID", "Shannon_Extrap")], 
                               by="SampleID", all=TRUE)
  
  rm(temp.diversity, temp.richness, temp, pseq.est)
}
 
# Calculate Faith's phylogenetic diversity for each sample and save it.
for(p in c(1:7)){ 
  temp.phylo.class = phy_tree(pseq.rarefied[[p]])
  temp.vegan.otu.table = as.data.frame(otu_table(pseq.rarefied[[p]]))
  picante.pd.result <- pd(temp.vegan.otu.table, temp.phylo.class, include.root=FALSE)
  
  picante.results <- data.frame(SampleID=rownames(picante.pd.result),
                                PD = picante.pd.result$PD)
  
  sample_data.mb[[p]] <- merge(sample_data.mb[[p]], picante.results, by="SampleID", all=TRUE)
  
  rm(picante.results, temp.phylo.class, temp.vegan.otu.table, picante.pd.result)
  
}


## [3] Calculate distance matrices ####
# Because there are so few lichen reads present in the LARS population, remove LARS from the lichen analysis.
lichens.at.LARS <- list(
  clr = subset_samples(pseq.clr$lichens, Group=="LARS"),
  rarefied = subset_samples(pseq.rarefied$lichens, Group=="LARS")
)

lichens.at.LARS[[1]] <- prune_taxa(taxa_sums(lichens.at.LARS[[1]]) > 0, lichens.at.LARS[[1]])
lichens.at.LARS[[2]] <- prune_taxa(taxa_sums(lichens.at.LARS[[2]]) > 0, lichens.at.LARS[[2]])

pseq.clr$lichens <- subset_samples(pseq.clr$lichens, Group !="LARS")
pseq.clr.sub$lichens <- subset_samples(pseq.clr.sub$lichens, Group !="LARS")
pseq.rarefied$lichens <- subset_samples(pseq.rarefied$lichens, Group !="LARS")
pseq.rarefied.sub$lichens <- subset_samples(pseq.rarefied.sub$lichens, Group !="LARS")

# Distance matrices for study populations.
dist.cb.list <- list()

for(p in c(1:4)){
  dist.cb.list[[p]] <- list()
  dist.cb.list[[p]][["Bray-Curtis"]] <- vegdist(as.data.frame(otu_table(pseq.rarefied[[p]])), method="bray")
  dist.cb.list[[p]][["Jaccard"]] <- vegdist(as.data.frame(otu_table(pseq.rarefied[[p]])), method="jaccard")
  dist.cb.list[[p]][["Aitchison"]] <- vegdist(as.data.frame(otu_table(microbiome::transform(pseq.clr[[p]], "clr"))), method="euclidean")
  dist.cb.list[[p]][["wUF"]] <- UniFrac(pseq.rarefied[[p]], weighted=TRUE, normalized=TRUE)
  dist.cb.list[[p]][["uwUF"]] <- UniFrac(pseq.rarefied[[p]], weighted=FALSE, normalized=TRUE)  
}

# For the protists, algae, and lichens, calculate distance matrices based on relative abundances.
for(p in c(5:8)){
  temp <- transform_sample_counts(pseq.clr.sub[[p]], function(x) 100*x/sum(x))
  
  dist.cb.list[[p]] <- list()
  dist.cb.list[[p]][["Bray-Curtis"]] <- vegdist(as.data.frame(otu_table(temp)), method="bray")
  dist.cb.list[[p]][["Jaccard"]] <- vegdist(as.data.frame(otu_table(temp)), method="jaccard")
  dist.cb.list[[p]][["Aitchison"]] <- vegdist(as.data.frame(otu_table(microbiome::transform(pseq.clr[[p]], "clr"))), method="euclidean")
  dist.cb.list[[p]][["wUF"]] <- UniFrac(temp, weighted=TRUE, normalized=TRUE)
  dist.cb.list[[p]][["uwUF"]] <- UniFrac(temp, weighted=FALSE, normalized=TRUE) 
  
  rm(temp)
}

names(dist.cb.list) <- c("X16S", "X18S", "FITS", "PITS", "protists", "algae", "lichens")

# Repeat for comparisons among deep-snow caribou herds.
dist.deepsnow.list <- list()

for(p in c(1:4)){
  temp.rarefied <- subset_samples(pseq.rarefied[[p]], Group=="deep-snow")
  temp.rarefied <- prune_taxa(taxa_sums(temp.rarefied) > 0, temp.rarefied)
  temp.clr <- subset_samples(pseq.clr[[p]], Group=="deep-snow")
  temp.clr <- prune_taxa(taxa_sums(temp.clr) > 0, temp.clr)
  
  dist.deepsnow.list[[p]] <- list()
  dist.deepsnow.list[[p]][["Bray-Curtis"]] <- vegdist(as.data.frame(otu_table(temp.rarefied)), method="bray")
  dist.deepsnow.list[[p]][["Jaccard"]] <- vegdist(as.data.frame(otu_table(temp.rarefied)), method="jaccard")
  dist.deepsnow.list[[p]][["Aitchison"]] <- vegdist(as.data.frame(otu_table(microbiome::transform(temp.clr, "clr"))), method="euclidean")
  dist.deepsnow.list[[p]][["wUF"]] <- UniFrac(temp.rarefied, weighted=TRUE, normalized=TRUE)
  dist.deepsnow.list[[p]][["uwUF"]] <- UniFrac(temp.rarefied, weighted=FALSE, normalized=TRUE)
  
  rm(temp.rarefied, temp.clr)
}

for(p in c(5:8)){
  temp <- transform_sample_counts(pseq.clr.sub[[p]], function(x) 100*x/sum(x))
  temp <- subset_samples(temp, Group=="deep-snow")
  temp <- prune_taxa(taxa_sums(temp) > 0, temp)
  
  temp2 <- subset_samples(pseq.clr.sub[[p]], Group=="deep-snow")
  temp2 <- prune_taxa(taxa_sums(temp2) > 0, temp2)
  
  dist.deepsnow.list[[p]] <- list()
  dist.deepsnow.list[[p]][["Bray-Curtis"]] <- vegdist(as.data.frame(otu_table(temp)), method="bray")
  dist.deepsnow.list[[p]][["Jaccard"]] <- vegdist(as.data.frame(otu_table(temp)), method="jaccard")
  dist.deepsnow.list[[p]][["Aitchison"]] <- vegdist(as.data.frame(otu_table(microbiome::transform(temp2, "clr"))), method="euclidean")
  dist.deepsnow.list[[p]][["wUF"]] <- UniFrac(temp, weighted=TRUE, normalized=TRUE)
  dist.deepsnow.list[[p]][["uwUF"]] <- UniFrac(temp, weighted=FALSE, normalized=TRUE) 
  
  rm(temp)
}

names(dist.deepsnow.list) <- c("X16S", "X18S", "FITS", "PITS", "protists", "algae", "lichens")

####### ALPHA & BETA DIVERSITY STATISTICS ####
## [1] Significant differences in alpha diversity metrics ####
# Run Levene's test, an ANOVA, and Tukey's post hoc test for each alpha diversity metric,
# for each amplicon, both among populations and among deep-snow caribou herds.
adiv.signif.all <- list() 
adiv.signif.deepsnow <- list() 

for(p in c(1:7)){
  adiv.signif.all[[p]] <- data.frame(matrix(nrow=0, ncol=12))
  adiv.signif.deepsnow[[p]] <- data.frame(matrix(nrow=0, ncol=21))
  
  # For each alpha diversity metric...
  for(i in c("Observed", "Observed_Extrap", "Shannon", "Shannon_Extrap", "PD")){
    # Perform the tests.
    levene <- car::leveneTest(sample_data.mb[[p]][,i] ~ Group, sample_data.mb[[p]])
    anova <- summary(aov(sample_data.mb[[p]][,i] ~ Group, sample_data.mb[[p]]))
    tukey <- TukeyHSD(aov(sample_data.mb[[p]][,i] ~ Group, sample_data.mb[[p]]))
    
    # Save the test results to the data frame.
    adiv.signif.all[[p]][i,1] <- i # variable name
    adiv.signif.all[[p]][i,2] <- levene$`F value`[1] # Levene's F
    adiv.signif.all[[p]][i,3] <- levene$`Pr(>F)`[1] # Levene's p
    adiv.signif.all[[p]][i,4] <- anova[[1]][1,4] # ANOVA F
    adiv.signif.all[[p]][i,5] <- anova[[1]][1,1] # ANOVA df
    adiv.signif.all[[p]][i,6] <- anova[[1]][1,5] # ANOVA p
    adiv.signif.all[[p]][i,7] <- tukey[[1]][1,4] # Tukey 1 DS - SS
    adiv.signif.all[[p]][i,8] <- tukey[[1]][2,4] # Tukey 2 DS - RP
    adiv.signif.all[[p]][i,9] <- tukey[[1]][3,4] # Tukey 3 DS - LARS
    adiv.signif.all[[p]][i,10] <- tukey[[1]][4,4] # Tukey 4 SS - RP
    adiv.signif.all[[p]][i,11] <- tukey[[1]][5,4] # Tukey 5 SS - LARS
    adiv.signif.all[[p]][i,12] <- tukey[[1]][6,4] # Tukey 6 RP - LARS
    
    # Clean workspace.
    rm(levene, anova, tukey)
    colnames(adiv.signif.all[[p]]) <- c("Measure", "Levene.F", "Levene.p", "ANOVA.F", "ANOVA.df", "ANOVA.p",
                                        "DS.SS", "DS.RP", "DS.LARS", "SS.RP", "SS.LARS", "RP.LARS")
    
    # Comparisons among deep snow caribou (separated by herd) for all amplicons except plants.
    if(p != 4){
      test_data <- subset(sample_data.mb[[p]], Group=="deep-snow")
      
      levene <- car::leveneTest(test_data[,i] ~ Herd, test_data)
      anova <- summary(aov(test_data[,i] ~ Herd, test_data))
      tukey <- TukeyHSD(aov(test_data[,i] ~ Herd, test_data))
      
      adiv.signif.deepsnow[[p]][i,1] <- i # variable name
      adiv.signif.deepsnow[[p]][i,2] <- levene$`F value`[1] # Levene's F
      adiv.signif.deepsnow[[p]][i,3] <- levene$`Pr(>F)`[1] # Levene's p
      adiv.signif.deepsnow[[p]][i,4] <- anova[[1]][1,4] # ANOVA F
      adiv.signif.deepsnow[[p]][i,5] <- anova[[1]][1,1] # ANOVA df
      adiv.signif.deepsnow[[p]][i,6] <- anova[[1]][1,5] # ANOVA p
      adiv.signif.deepsnow[[p]][i,7] <- tukey[[1]][1,4] # Tukey 1 BV - Selkirk
      adiv.signif.deepsnow[[p]][i,8] <- tukey[[1]][2,4] # Tukey 2 BV - Hart
      adiv.signif.deepsnow[[p]][i,9] <- tukey[[1]][3,4] # Tukey 3 BV - North Cariboo
      adiv.signif.deepsnow[[p]][i,10] <- tukey[[1]][4,4] # Tukey 4 BV - North Columbia
      adiv.signif.deepsnow[[p]][i,11] <- tukey[[1]][5,4] # Tukey 5 BV - Wells Grey
      adiv.signif.deepsnow[[p]][i,12] <- tukey[[1]][6,4] # Tukey 6 Selkirk - Hart
      adiv.signif.deepsnow[[p]][i,13] <- tukey[[1]][7,4] # Tukey 7 Selkirk - North Cariboo
      adiv.signif.deepsnow[[p]][i,14] <- tukey[[1]][8,4] # Tukey 8 Selkirk - North Columbia
      adiv.signif.deepsnow[[p]][i,15] <- tukey[[1]][9,4] # Tukey 9 Selkirk - Wells Grey
      adiv.signif.deepsnow[[p]][i,16] <- tukey[[1]][10,4] # Tukey 10 Hart - North Cariboo
      adiv.signif.deepsnow[[p]][i,17] <- tukey[[1]][11,4] # Tukey 11 Hart - North Columbia
      adiv.signif.deepsnow[[p]][i,18] <- tukey[[1]][12,4] # Tukey 12 Hart - Wells Grey
      adiv.signif.deepsnow[[p]][i,19] <- tukey[[1]][13,4] # Tukey 13 North Cariboo - North Columbia
      adiv.signif.deepsnow[[p]][i,20] <- tukey[[1]][14,4] # Tukey 14 North Cariboo - Wells Grey 
      adiv.signif.deepsnow[[p]][i,21] <- tukey[[1]][15,4] # Tukey 15 North Columbia - Wells Grey 
      
      rm(levene, anova, tukey, test_data)
      
      colnames(adiv.signif.deepsnow[[p]]) <- c("Measure", "Levene.F", "Levene.p", "ANOVA.F", "ANOVA.df", "ANOVA.p",
                                          "BV.SK", "BV.Hart", "BV.NCoo", "BV.NCmb", "BV.WG",
                                          "SK.Hart", "SK.NCoo", "SK.NCmb", "SK.WG",
                                          "Hart.NCoo", "Hart.NCmb", "Hart.WG",
                                          "NCoo.NCmb", "NCoo.WG", "NCmb.WG")
    }
  }
}

names(adiv.signif.all) <- c("X16S", "X18S", "FITS", "PITS", "protists", "algae", "lichens")
names(adiv.signif.deepsnow) <- c("X16S", "X18S", "FITS", "PITS", "protists", "algae", "lichens")

# Add the name of the amplicon to each data frame.
for(j in c(1:7)){
  adiv.signif.all[[j]]$amplicon <- names(adiv.signif.all)[j]
  if(j !=4){adiv.signif.deepsnow[[j]]$amplicon <- names(adiv.signif.deepsnow)[j]}
}

# Bind all the data frames together for one table of significance tests.
adiv.signif.all <- as.data.frame(dplyr::bind_rows(adiv.signif.all))
adiv.signif.deepsnow <- as.data.frame(dplyr::bind_rows(adiv.signif.deepsnow))

## [2] Distance-based ordinations ####
bdiv.ordinations.all <- list()
bdiv.ordinations.deepsnow <- list()

# Note that ordinations are saved by this hierarchy:
# bdiv.ordinations -> amplicon -> data type (models / points / ellipses) -> 
# distance metric (Bray, Jaccard, Aitchison, wUF, uwUF)

# Define a function for calculating ellipses
veganCovEllipse <-function (cov, center = c(0, 0), scale = 1, npoints = 100){
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

for(p in c(1:7)){ # For every data set...
  bdiv.ordinations.all[[p]] <- list()
  bdiv.ordinations.all[[p]][["models"]] <- list()
  bdiv.ordinations.all[[p]][["points"]] <- list()
  bdiv.ordinations.all[[p]][["ellipses"]] <- list()
  
  bdiv.ordinations.deepsnow[[p]] <- list()
  bdiv.ordinations.deepsnow[[p]][["models"]] <- list()
  bdiv.ordinations.deepsnow[[p]][["points"]] <- list()
  
  # Calculate PCoA for the Bray-Curtis, Jaccard, weighted and unweighted UniFrac distances.
  for(i in c(1:5)){
    # Run the PCoA model for each of the non-Aitchison distances.
    if(i !=3){
      bdiv.ordinations.all[[p]][["models"]][[i]] <- capscale(dist.cb.list[[p]][[i]]~1)
      bdiv.ordinations.deepsnow[[p]][["models"]][[i]] <- capscale(dist.deepsnow.list[[p]][[i]]~1)
    }
    # Or run the model for the Aitchison distances.
    if(i==3){
      bdiv.ordinations.all[[p]][["models"]][[i]] <- 
        rda(as.data.frame(otu_table(microbiome::transform(pseq.clr[[p]], "clr"))))
      bdiv.ordinations.deepsnow[[p]][["models"]][[i]] <- 
        rda(as.data.frame(otu_table(microbiome::transform(subset_samples(pseq.clr[[p]], Group=="deep-snow"), "clr"))))
      
    }
    
    # Extract the first two axes loadings.
    bdiv.ordinations.all[[p]][["points"]][[i]] = 
      merge(scores(bdiv.ordinations.all[[p]][["models"]][[i]],
                   display="sites", choices=c(1:2)),
            subset(sample_data(pseq.clr[[p]]), 
                   SampleID %in% rownames(scores(bdiv.ordinations.all[[p]][["models"]][[i]],
                                                 display="sites", choices=c(1:2))))[,c("SampleID", "Group","Herd")],
            by=0, all=TRUE)
    
    bdiv.ordinations.deepsnow[[p]][["points"]][[i]] = 
      merge(scores(bdiv.ordinations.deepsnow[[p]][["models"]][[i]],
                   display="sites", choices=c(1:2)),
            subset(sample_data(pseq.clr[[p]]), Group=="deep-snow" & 
                     SampleID %in% rownames(scores(bdiv.ordinations.deepsnow[[p]][["models"]][[i]],
                                                   display="sites", choices=c(1:2))))[,c("SampleID", "Group", "Herd")],
            by=0, all=TRUE)
    
    bdiv.ordinations.deepsnow[[p]][["points"]][[i]]$Herd <- 
      factor(bdiv.ordinations.deepsnow[[p]][["points"]][[i]]$Herd)
    
    # Calculate confidence ellipses for the all-types comparison.
    plot.new()
    temp <- ordiellipse(bdiv.ordinations.all[[p]][["models"]][[i]], 
                        bdiv.ordinations.all[[p]][["points"]][[i]]$Group, 
                        display="sites", 
                        kind="sd", 
                        conf=0.90, 
                        label=T)
                        
    bdiv.ordinations.all[[p]][["ellipses"]][[i]] <- data.frame()
    
    for(g in levels(bdiv.ordinations.all[[p]][["points"]][[i]]$Group)){
      bdiv.ordinations.all[[p]][["ellipses"]][[i]] <- rbind(bdiv.ordinations.all[[p]][["ellipses"]][[i]],
                                                            cbind(as.data.frame(with(bdiv.ordinations.all[[p]][["points"]][[i]][bdiv.ordinations.all[[p]][["points"]][[i]]$Group==g,],
                                                                                     veganCovEllipse(temp[[g]]$cov,
                                                                                                     temp[[g]]$center,
                                                                                                     temp[[g]]$scale)))
                                                                  ,Group=g))
    }
    
  }
}

# Name all the objects that were just created.
names(bdiv.ordinations.all) <- c("X16S", "X18S", "FITS", "PITS", "protists", "algae", "lichens")
names(bdiv.ordinations.deepsnow) <- c("X16S", "X18S", "FITS", "PITS", "protists", "algae", "lichens")

for(p in c(1:7)){
  for(j in c(1:3)){
    names(bdiv.ordinations.all[[p]][[j]]) <- c("BrayCurtis", "Jaccard", "Aitchison",
                                               "weightedUF", "unweightedUF")
    
    if(p !=4){
      names(bdiv.ordinations.deepsnow[[p]][[j]]) <- c("BrayCurtis", "Jaccard", "Aitchison",
                                                      "weightedUF", "unweightedUF")}
  }
}

## [3] PERMANOVA and homogeneity of multivariate dispersion ####
bdiv.permanova <- list()
bdiv.betadisper <- list()
bdiv.distances <- list()
bdiv.betadisper.pw <- list()

for(p in c(1:7)){
  # Prepare data frames to receive the PERMANOVA data.
  bdiv.permanova[[p]] <- data.frame(matrix(ncol=6, nrow=10))
  bdiv.permanova[[p]][,1] <- rep(c("Bray-Curtis", "Jaccard", "Aitchison", "wUF", "uwUF"), 2)
  bdiv.permanova[[p]][,2] <- c(rep("all samples", 5), rep("deepsnow only", 5))
  colnames(bdiv.permanova[[p]]) <- c("Metric", "Comparison", "F", "df", "R2", "p")
  
  # Run the PERMANOVA analysis.
  for(i in c(1:10)){
    if(i %in% c(1:5)){ # By study population.
      temp <- adonis(dist.cb.list[[p]][[i]] ~ Group, as.data.frame(cbind(subset(sample_data(pseq.clr[[p]]),
                                                                                SampleID %in% rownames(as.matrix(dist.cb.list[[p]][[i]]))))), permutations=1000)
      bdiv.permanova[[p]][i,3] <- temp[["aov.tab"]]$F.Model[1]
      bdiv.permanova[[p]][i,4] <- temp[["aov.tab"]]$Df[1]
      bdiv.permanova[[p]][i,5] <- temp[["aov.tab"]]$R2[1]
      bdiv.permanova[[p]][i,6] <- temp[["aov.tab"]]$`Pr(>F)`[1]
    }
    if(i %in% c(6:10)){ # By deep-snow caribou herd.
      temp <- adonis(dist.deepsnow.list[[p]][[i-5]] ~ Herd, subset(as.data.frame(cbind(sample_data(pseq.clr[[p]]))), Group=="deep-snow" & 
                                                                     SampleID %in% rownames(as.matrix(dist.cb.list[[p]][[i-5]]))), permuations=1000)
      bdiv.permanova[[p]][i,3] <- temp[["aov.tab"]]$F.Model[1]
      bdiv.permanova[[p]][i,4] <- temp[["aov.tab"]]$Df[1]
      bdiv.permanova[[p]][i,5] <- temp[["aov.tab"]]$R2[1]
      bdiv.permanova[[p]][i,6] <- temp[["aov.tab"]]$`Pr(>F)`[1]
    }
    rm(temp)
  }
  
  # Prepare data frames to receive the betadispersion data.
  # Initial significance tests.
  bdiv.betadisper[[p]] <- data.frame(matrix(ncol=5, nrow=10))
  bdiv.betadisper[[p]][,1] <- rep(c("Bray-Curtis", "Jaccard", "Aitchison", "wUF", "uwUF"), 2)
  bdiv.betadisper[[p]][,2] <- c(rep("all samples", 5), rep("deepsnow only", 5))
  colnames(bdiv.betadisper[[p]]) <- c("Metric", "Comparison", "F", "df", "p")
  
  # Actual distances to centroid for each sample.
  bdiv.distances[[p]] <- data.frame(SampleID = sample_data$SampleID)
  
  # Pairwise comparison results.
  if(p<7){bdiv.betadisper.pw[[p]] <- data.frame(matrix(nrow=6))}
  if(p==7){bdiv.betadisper.pw[[p]] <- data.frame(matrix(nrow=3))}
  
  # Run the beta dispersion analysis.
  for(i in c(1:5)){
    if(i %in% c(1:5)){ # Among study populations.
      disp <- betadisper(dist.cb.list[[p]][[i]], 
                         as.data.frame(cbind(subset(sample_data(pseq.clr[[p]]),
                                                    SampleID %in% rownames(as.matrix(dist.cb.list[[p]][[i]])))))$Group)
      temp <- permutest(disp,
                        permuations=1000, pairwise=TRUE)
      
      # Test for significant differences in multivariate dispersion.
      bdiv.betadisper[[p]][i,3] <- temp[["tab"]][1,4]
      bdiv.betadisper[[p]][i,4] <- temp[["tab"]][1,1]
      bdiv.betadisper[[p]][i,5] <- temp[["tab"]][1,6]
      
      # Save the actual measures of distance to centroid.
      temp2 <- as.data.frame(disp$distances)
      temp2$SampleID <- rownames(temp2)
      bdiv.distances[[p]] <- merge(bdiv.distances[[p]], temp2, by="SampleID", all=TRUE)
      
      # Save the pairwise comparisons for beta dispersion.
      if(i==1){bdiv.betadisper.pw[[p]]$comparison <- names(temp$pairwise$observed)}
      temp.p <- as.data.frame(temp$pairwise$observed)
      temp.p$comparison <- rownames(temp.p)
      temp.stat <- as.data.frame(temp$statistic[c(2:length(temp$statistic))])
      temp.stat$comparison <- rownames(temp.p)
      
      bdiv.betadisper.pw[[p]] <- merge(bdiv.betadisper.pw[[p]], temp.stat,
                                       by="comparison", all=TRUE)
      bdiv.betadisper.pw[[p]] <- merge(bdiv.betadisper.pw[[p]], temp.p,
                                       by="comparison", all=TRUE)
      
    }
    if(i %in% c(6:10)){
      disp <- betadisper(dist.deepsnow.list[[p]][[i-5]], 
                         subset(as.data.frame(cbind(sample_data(pseq.clr[[p]]))),
                                Group=="deep-snow" & SampleID %in% rownames(as.matrix(dist.cb.list[[p]][[i-5]])))$Herd)
      temp <- permutest(disp,
                        permuations=1000)
      
      bdiv.betadisper[[p]][i,3] <- temp[["tab"]][1,4]
      bdiv.betadisper[[p]][i,4] <- temp[["tab"]][1,1]
      bdiv.betadisper[[p]][i,5] <- temp[["tab"]][1,6]
      
      temp2 <- as.data.frame(disp$distances)
      temp2$SampleID <- rownames(temp2)
      bdiv.distances[[p]] <- merge(bdiv.distances[[p]], temp2, by="SampleID", all=TRUE)
    }
    rm(temp, temp2, disp)
  }
}

names(bdiv.permanova) <- c("X16S", "X18S", "FITS", "PITS", "protists", "algae", "lichens")
names(bdiv.betadisper) <- c("X16S", "X18S", "FITS", "PITS", "protists", "algae", "lichens")
names(bdiv.distances) <- c("X16S", "X18S", "FITS", "PITS", "protists", "algae", "lichens")
names(bdiv.betadisper.pw) <- c("X16S", "X18S", "FITS", "PITS", "protists", "algae", "lichens")

# Clean pairwise comparison results into a single data frame.
for(p in c(1:7)){
  bdiv.betadisper.pw[[p]][,2] <- NULL
  t1 <- bdiv.betadisper.pw[[p]][,c(1,2:3)]
  t2 <- bdiv.betadisper.pw[[p]][,c(1,4:5)]
  t3 <- bdiv.betadisper.pw[[p]][,c(1,6:7)]
  t4 <- bdiv.betadisper.pw[[p]][,c(1,8:9)]
  t5 <- bdiv.betadisper.pw[[p]][,c(1,10:11)]
  
  colnames(t1) <- c("comparison", "t", "p")
  colnames(t2) <- c("comparison", "t", "p")
  colnames(t3) <- c("comparison", "t", "p")
  colnames(t4) <- c("comparison", "t", "p")
  colnames(t5) <- c("comparison", "t", "p")
  
  bdiv.betadisper.pw[[p]] <- rbind(t1, t2, t3, t4, t5)
  rm(t1, t2, t3, t4, t5)
  
  if(p<7){bdiv.betadisper.pw[[p]]$distance <- c(rep("BrayCurtis", 6),
                                        rep("Jaccard", 6),
                                        rep("Aitchison", 6),
                                        rep("weightedUF", 6),
                                        rep("unweightedUF", 6))}
  if(p==7){bdiv.betadisper.pw[[p]]$distance <- c(rep("BrayCurtis", 3),
                                                rep("Jaccard", 3),
                                                rep("Aitchison", 3),
                                                rep("weightedUF", 3),
                                                rep("unweightedUF", 3))}
}

temp <- dplyr::bind_cols(bdiv.betadisper.pw[c(1:7)])

## [4] Pairwise PERMANOVA among study populations ####
bdiv.pairwise <- list()

devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

for(p in c(1:8)){
  bdiv.pairwise[[p]] <- list()
  
  for(i in c(1:5)){
    temp <- pairwiseAdonis::pairwise.adonis2(
      dist.cb.list[[p]][[i]] ~ Group, 
      as.data.frame(cbind(subset(sample_data(pseq.clr[[p]]),
                                 SampleID %in% rownames(as.matrix(dist.cb.list[[p]][[i]]))))))
    for(j in c(2:length(temp))){
      temp[[j]] <- temp[[j]][1,]
      temp[[j]]$comparison <- names(temp)[j]
    }
    
    temp <- dplyr::bind_rows(temp[c(2:length(temp))])
    temp$amplicon <- names(pseq.clr.sub)[p]
    
    if(i==1){temp$distance <- "BrayCurtis"}
    if(i==2){temp$distance <- "Jaccard"}
    if(i==3){temp$distance <- "Aitchison"}
    if(i==4){temp$distance <- "weightedUF"}
    if(i==5){temp$distance <- "unweightedUF"}
    
    bdiv.pairwise[[p]][[i]] <- temp
    rm(temp)
  }
}

names(bdiv.pairwise) <- names(pseq.clr.sub)

for(p in c(1:8)){
  bdiv.pairwise[[p]] <- dplyr::bind_rows(bdiv.pairwise[[p]])
}

####### TAXON-BASED ANALYSES ####
## [1] Differential abundance among study populations ####
library(ALDEx2)

diff.abund.GROUP <- list()

for(p in c(1:7)){ # For each data set...
  diff.abund.sets <- list()
  
  # Define the taxonomic ranks being used for comparison.
  if(p %in% c(1,3,5,6,7)){taxranks <- c("Phylum", "Class", "Order", "Family", "Genus")}
  if(p == 2){taxranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")}
  if(p == 4){taxranks <- c("Class", "Order", "Family", "Genus")}
  
  for(i in c(1:(length(taxranks)+1))){ # For each taxonomic rank...
    temp.physeq <- pseq.clr.sub[[p]]
    
    # Agglomerate the phyloseq object to the taxonomic rank of interest.
    if(i < (length(taxranks)+1)){
      temp.physeq <- tax_glom(temp.physeq, taxrank=taxranks[[i]]) 
      taxa_names(temp.physeq) <- as.data.frame(cbind(tax_table(temp.physeq)))[,taxranks[[i]]]
    }
    
    # Extract the OTU table so that taxa are rows
    temp.otu.table <- as.data.frame(otu_table(temp.physeq))
    if(taxa_are_rows(temp.physeq)=="FALSE"){temp.otu.table <- as.data.frame(t(temp.otu.table))}
    
    # Define the covariate (variable of interest for the comparison)
    covariates <- as.character(sample_data(temp.physeq)$Group)
    
    # Run the differential abundance analysis.
    result <- aldex(temp.otu.table, covariates, mc.samples=128, test="kw", effect=FALSE, denom="all")
    
    # Obtain a new OTU table and taxonomy table based on relative abundance.
    temp.physeq <- transform_sample_counts(temp.physeq, function(x) 100*x/sum(x))  
    temp.tax <- as.data.frame(cbind(tax_table(temp.physeq)))
    temp.otu.table <- as.data.frame(otu_table(temp.physeq))
    if(taxa_are_rows(temp.physeq)=="TRUE"){temp.otu.table <- as.data.frame(t(temp.otu.table))}
    
    temp.otu.table <- cbind(covariates, temp.otu.table)
    
    # Calculate the mean (sd) and median for each taxon in each tested covariate group.
    table <- temp.otu.table %>%
      dplyr::group_by(covariates) %>%
      dplyr::summarise_if(is.numeric, .funs=c(mn = mean, sd = sd)) %>%
      tibble::column_to_rownames("covariates") %>%
      t() %>%
      as.data.frame() %>%
      tibble::rownames_to_column("covariates") %>%
      tidyr::separate(., "covariates", into=c("ASV","fnctn"), sep=-3, remove=TRUE, convert=FALSE) %>%
      dplyr::group_split(fnctn) %>%
      purrr::map(. %>%
                   tibble::column_to_rownames("ASV") %>%
                   as.data.frame()) %>%
      as.list()
    
    for(j in c(1:length(table))){
      table[[j]] <- as.data.frame(table[[j]])
      colnames(table[[j]]) <- paste0(colnames(table[[j]]), table[[j]][1,"fnctn"])
      table[[j]][,1] <- NULL
    }
    
    table <- dplyr::bind_cols(table)
    
    final <- transform(merge(result, temp.tax, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
    final <- transform(merge(final, table, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
    
    diff.abund.sets[[i]] <- final
    
    rm(final, result, temp.physeq, temp.tax, temp.otu.table, covariates, table)
  }
  
  names(diff.abund.sets) <- c(taxranks, "ASV")
  diff.abund.GROUP[[p]] <- diff.abund.sets
  rm(diff.abund.sets)
}

names(diff.abund.GROUP) <- c("X16S", "X18S", "FITS", "PITS",
                             "protists", "algae", "lichens")

## [2] Differential abundance among deep-snow caribou herds ####
diff.abund.HERD <- list()

for(p in c(1:7)){ # For each of the four amplicons...
  diff.abund.sets <- list()
  
  if(p %in% c(1,3,5,6,7)){taxranks <- c("Phylum", "Class", "Order", "Family", "Genus")}
  if(p == 2){taxranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")}
  if(p == 4){taxranks <- c("Class", "Order", "Family", "Genus")}
  
  for(i in c(1:(length(taxranks)+1))){ # For each taxonomic rank...
    temp.physeq <- subset_samples(pseq.clr.sub[[p]], Group=="deep-snow")
    temp.physeq <- prune_taxa(taxa_sums(temp.physeq) > 0, temp.physeq)
    
    if(i < (length(taxranks)+1)){
      temp.physeq <- tax_glom(temp.physeq, taxrank=taxranks[[i]]) 
      taxa_names(temp.physeq) <- as.data.frame(cbind(tax_table(temp.physeq)))[,taxranks[[i]]]
    }
    
    # Extract the OTU table so that taxa are rows
    temp.otu.table <- as.data.frame(otu_table(temp.physeq))
    if(taxa_are_rows(temp.physeq)=="FALSE"){temp.otu.table <- as.data.frame(t(temp.otu.table))}
    
    # Define the covariate (variable of interest for the comparison)
    covariates <- as.character(sample_data(temp.physeq)$Herd)
    
    # Run the differential abundance analysis.
    result <- aldex(temp.otu.table, covariates, mc.samples=128, test="kw", effect=FALSE, denom="all")
    
    # Obtain a new OTU table and taxonomy table based on relative abundance.
    temp.physeq <- transform_sample_counts(temp.physeq, function(x) 100*x/sum(x))  
    temp.tax <- as.data.frame(cbind(tax_table(temp.physeq)))
    temp.otu.table <- as.data.frame(otu_table(temp.physeq))
    if(taxa_are_rows(temp.physeq)=="TRUE"){temp.otu.table <- as.data.frame(t(temp.otu.table))}
    
    temp.otu.table <- cbind(covariates, temp.otu.table)
    
    # Calculate the mean (sd) and median for each taxon in each tested covariate group.
    table <- temp.otu.table %>%
      dplyr::group_by(covariates) %>%
      dplyr::summarise_if(is.numeric, .funs=c(mn = mean, sd = sd)) %>%
      tibble::column_to_rownames("covariates") %>%
      t() %>%
      as.data.frame() %>%
      tibble::rownames_to_column("covariates") %>%
      tidyr::separate(., "covariates", into=c("ASV","fnctn"), sep=-3, remove=TRUE, convert=FALSE) %>%
      dplyr::group_split(fnctn) %>%
      purrr::map(. %>%
                   tibble::column_to_rownames("ASV") %>%
                   as.data.frame()) %>%
      as.list()
    
    for(j in c(1:length(table))){
      table[[j]] <- as.data.frame(table[[j]])
      colnames(table[[j]]) <- paste0(colnames(table[[j]]), table[[j]][1,"fnctn"])
      table[[j]][,1] <- NULL
    }
    
    table <- dplyr::bind_cols(table)
    
    final <- transform(merge(result, temp.tax, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
    final <- transform(merge(final, table, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
    
    diff.abund.sets[[i]] <- final
    
    rm(final, result, temp.physeq, temp.tax, temp.otu.table, covariates, table)
  }
  
  names(diff.abund.sets) <- c(taxranks, "ASV")
  diff.abund.GROUP[[p]] <- diff.abund.sets
  rm(diff.abund.sets)
}

names(diff.abund.HERD) <- c("X16S", "X18S", "FITS", "PITS",
                            "protists", "algae", "lichens")

## [3] Random forest models predicting study population assignments ####
rf.results.all <- list()

# This code is adapted from https://rpubs.com/michberr/randomforestmicrobe.
for(p in c(1:7)){
  rf.results.all[[p]] <- list()
  rf.results.all[[p]][["Models"]] <- list()
  rf.results.all[[p]][["Gini_scores"]] <- list()
  
  # Identify taxa with prevalence > 4 across the data and mean relative abundance >0.005%
  goodTaxa <- subset(prev.data[[p]]$Genus, Prevalence > 4 & Abundance > 0.005)
  
  go <- tax_glom(pseq.clr.sub[[p]], taxrank="Genus")
  taxa_names(go) <- as.data.frame(cbind(tax_table(go)))$Genus
  
  sites.prune <- subset_taxa(go, taxa_names(go) %in% goodTaxa$Taxon)
  sites.prune <- microbiome::transform(sites.prune, "clr")
  
  # Define the model predictors as ASV abundances.
  predictors <- otu_table(sites.prune)
  
  # Define the model response as the variable of interest (segment, individual identity, or intestinal site).
  response <- sample_data(sites.prune)$Group
  
  # Combine the predictors and response into a single data frame.
  rf.data <- data.frame(response, predictors)
  
  # Perform random forest model.
  set.seed(2)
  rf.raw.results <- randomForest(response~., data = rf.data, ntree = 1000)
  
  # Extract the 100 most important ASVs, based on their mean decrease in the Gini coefficient.
  rf.importance <- randomForest::importance(rf.raw.results)
  rf.importance <- data.frame(predictors = rownames(rf.importance), rf.importance)
  rf.importance <- dplyr::arrange(rf.importance, plyr::desc(MeanDecreaseGini))
  rf.importance$predictors <- factor(rf.importance$predictors,
                                     levels = rf.importance$predictors)
  if(ntaxa(sites.prune) >= 100){rf.importance <- rf.importance[1:100, ]}
  rownames(rf.importance) <- rf.importance$predictors
  rf.importance$predictors <- NULL
  
  # Add taxonomy information to the predictors.
  temp <- as.data.frame(cbind(tax_table(sites.prune)))
  temp$Taxon <- rownames(temp)
  
  rf.importance$Taxon <- rownames(rf.importance)
  rf.importance$Taxon <- gsub(".", "/", rf.importance$Taxon, fixed=TRUE)
  
  temp <- subset(temp, temp$Taxon %in% rf.importance$Taxon)
  rf.importance <- merge(rf.importance, temp, by="Taxon", all=TRUE)
  
  if(p < 7){rf.importance <- merge(rf.importance, 
                                   prev.data[[p]]$Genus[,c(1:7)],
                                   by="Taxon", all=FALSE)}
  if(p==7){rf.importance <- merge(rf.importance, 
                                  prev.data[[p]]$Genus[,c(1:6)],
                                  by="Taxon", all=FALSE)}
  
  # Organize by Gini score.  
  rf.importance <- dplyr::arrange(rf.importance, plyr::desc(MeanDecreaseGini))
  
  # Save results and clean workspace.
  rf.results.all[[p]][["Models"]] <- rf.raw.results
  rf.results.all[[p]][["Gini_scores"]] <- rf.importance
  
  rm(sites.prune, predictors, 
     response, rf.data, temp, rf.raw.results, rf.importance, go)
}


## [4] Indicator species analysis ####
indicator.species <- list()

for(p in c(1:7)){
  indicator.species[[p]] <- list()
  
  # Agglomerate taxa to the genus level and convert to relative abundances.
    temp.otu <- tax_glom(pseq.clr.sub[[p]], taxrank="Genus")
    temp.otu <- transform_sample_counts(temp.otu, function(x) 100*x/sum(x))
    taxa_names(temp.otu) <- as.data.frame(cbind(tax_table(temp.otu)))$Genus
    
    # Subset to only genera with prevalence > 4 and mean relative abundance > 0.005%. 
    goodTaxa <- subset(prev.data[[p]]$Genus, Prevalence > 4 & Abundance > 0.005)
    temp.otu <- subset_taxa(temp.otu, taxa_names(temp.otu) %in% goodTaxa$Taxon)
    
    # Extract OTU table and group labels.
    temp.otu <- as.data.frame(cbind(otu_table(temp.otu)))
    grps <- sample_data(pseq.clr.sub[[p]])$Group
    
    # Run indicator species analysis.
    indval <- indicspecies::multipatt(temp.otu, grps, control=how(nperm=999))
  
    # Save tables of 'A', 'B', and 'stat' values.
    tableA <- indval$A
    tableB <- indval$B
    tableS <- indval$sign
    
    # Create separate tables for each indicator targets.
    # (e.g., a table for all taxa indicating deep-snow caribou, another table for all indicators of shallow-snow...
    # ... as well as tables for multi-target indicators (e.g., deep-snow & shallow-snow or deep-snow & Rev. pen)
    tableS <- subset(tableS, is.na(p.value)=="FALSE")
    
    if(p < 7){
      tableS.CM <- subset(tableS, tableS$`s.shallow-snow`==1 &
                            tableS$`s.LARS`==0 &
                            tableS$`s.deep-snow`==0 &
                            tableS$`s.Revelstoke pen`==0)
      tableS.SM <- subset(tableS, tableS$`s.shallow-snow`==0  &
                            tableS$`s.LARS`==0 &
                            tableS$`s.deep-snow`==1 &
                            tableS$`s.Revelstoke pen`==0)
      tableS.RP <- subset(tableS, tableS$`s.shallow-snow`==0  &
                            tableS$`s.LARS`==0 &
                            tableS$`s.deep-snow`==0 &
                            tableS$`s.Revelstoke pen`==1)
      tableS.wild <- subset(tableS, tableS$`s.shallow-snow`==1  &
                              tableS$`s.LARS`==0 &
                              tableS$`s.deep-snow`==1 &
                              tableS$`s.Revelstoke pen`==0)
      tableS.deepsnow <- subset(tableS, tableS$`s.shallow-snow`==0  &
                                  tableS$`s.LARS`==0 &
                                  tableS$`s.deep-snow`==1 &
                                  tableS$`s.Revelstoke pen`==1)
      tableS.LARS <- subset(tableS, tableS$`s.shallow-snow`==0 &
                              tableS$`s.LARS`==1 &
                              tableS$`s.deep-snow`==0 &
                              tableS$`s.Revelstoke pen`==0)
    }
    
    if(p==7){
      tableS.CM <- subset(tableS, tableS$`s.shallow-snow`==1 &
                            tableS$`s.deep-snow`==0 &
                            tableS$`s.Revelstoke pen`==0)
      tableS.SM <- subset(tableS, tableS$`s.shallow-snow`==0  &
                            tableS$`s.deep-snow`==1 &
                            tableS$`s.Revelstoke pen`==0)
      tableS.RP <- subset(tableS, tableS$`s.shallow-snow`==0  &
                            tableS$`s.deep-snow`==0 &
                            tableS$`s.Revelstoke pen`==1)
      tableS.wild <- subset(tableS, tableS$`s.shallow-snow`==1  &
                              tableS$`s.deep-snow`==1 &
                              tableS$`s.Revelstoke pen`==0)
      tableS.deepsnow <- subset(tableS, tableS$`s.shallow-snow`==0  &
                                  tableS$`s.deep-snow`==1 &
                                  tableS$`s.Revelstoke pen`==1)
    }
    
    # For each of these indicator tables, add the 'A' and 'B' values for the target groups. ####
    # Shallow-snow caribou
    if(nrow(tableS.CM) > 0){
      tableS.CM$indicatee <- "shallow-snow"
      tableS.CM$taxon <- rownames(tableS.CM)
      
      tableA.CM <- as.data.frame(subset(tableA, rownames(tableA) %in% rownames(tableS.CM)))
      tableA.CM$taxon <- rownames(tableA.CM)
      
      tableS.CM <- merge(tableS.CM[,c("indicatee", "stat", "p.value", "taxon")],
                         as.data.frame(tableA.CM[,c("taxon", "shallow-snow")]),
                         by="taxon", all=TRUE)
      colnames(tableS.CM) <- c("taxon", "indicatee", "stat", "p.value", "A")
      
      tableB.CM <- as.data.frame(subset(tableB, rownames(tableB) %in% tableS.CM$taxon))
      tableB.CM$taxon <- rownames(tableB.CM)
      
      tableS.CM <- merge(tableS.CM,
                         tableB.CM[,c("taxon", "shallow-snow")],
                         by="taxon", all=TRUE)
      colnames(tableS.CM) <- c("taxon", "indicatee", "stat", "p.value", "A", "B")
      
      rm(tableA.CM, tableB.CM)
    }
    
    # Deep-snow
    if(nrow(tableS.SM) > 0){
      tableS.SM$indicatee <- "deep-snow"
      tableS.SM$taxon <- rownames(tableS.SM)
      
      tableA.SM <- as.data.frame(subset(tableA, rownames(tableA) %in% rownames(tableS.SM)))
      tableA.SM$taxon <- rownames(tableA.SM)
      
      tableS.SM <- merge(tableS.SM[,c("indicatee", "stat", "p.value", "taxon")],
                         as.data.frame(tableA.SM[,c("taxon", "deep-snow")]),
                         by="taxon", all=TRUE)
      colnames(tableS.SM) <- c("taxon", "indicatee", "stat", "p.value", "A")
      
      tableB.SM <- as.data.frame(subset(tableB, rownames(tableB) %in% tableS.SM$taxon))
      tableB.SM$taxon <- rownames(tableB.SM)
      
      tableS.SM <- merge(tableS.SM,
                         tableB.SM[,c("taxon", "deep-snow")],
                         by="taxon", all=TRUE)
      colnames(tableS.SM) <- c("taxon", "indicatee", "stat", "p.value", "A", "B")
      
      rm(tableA.SM, tableB.SM)
    }
    
    # Revelstoke pen
    if(nrow(tableS.RP) > 0){
      tableS.RP$indicatee <- "Revelstoke pen"
      tableS.RP$taxon <- rownames(tableS.RP)
      
      tableA.RP <- as.data.frame(subset(tableA, rownames(tableA) %in% rownames(tableS.RP)))
      tableA.RP$taxon <- rownames(tableA.RP)
      
      tableS.RP <- merge(tableS.RP[,c("indicatee", "stat", "p.value", "taxon")],
                         as.data.frame(tableA.RP[,c("taxon", "Revelstoke pen")]),
                         by="taxon", all=TRUE)
      colnames(tableS.RP) <- c("taxon", "indicatee", "stat", "p.value", "A")
      
      tableB.RP <- as.data.frame(subset(tableB, rownames(tableB) %in% tableS.RP$taxon))
      tableB.RP$taxon <- rownames(tableB.RP)
      
      tableS.RP <- merge(tableS.RP,
                         tableB.RP[,c("taxon", "Revelstoke pen")],
                         by="taxon", all=TRUE)
      colnames(tableS.RP) <- c("taxon", "indicatee", "stat", "p.value", "A", "B")
      
      rm(tableA.RP, tableB.RP)
    }
    
    # LARS
    if(p < 7){
      if(nrow(tableS.LARS) > 0){
      tableS.LARS$indicatee <- "LARS"
      tableS.LARS$taxon <- rownames(tableS.LARS)
      
      tableA.LARS <- as.data.frame(subset(tableA, rownames(tableA) %in% rownames(tableS.LARS)))
      tableA.LARS$taxon <- rownames(tableA.LARS)
      
      tableS.LARS <- merge(tableS.LARS[,c("indicatee", "stat", "p.value", "taxon")],
                           as.data.frame(tableA.LARS[,c("taxon", "LARS")]),
                           by="taxon", all=TRUE)
      colnames(tableS.LARS) <- c("taxon", "indicatee", "stat", "p.value", "A")
      
      tableB.LARS <- as.data.frame(subset(tableB, rownames(tableB) %in% tableS.LARS$taxon))
      tableB.LARS$taxon <- rownames(tableB.LARS)
      
      tableS.LARS <- merge(tableS.LARS,
                           tableB.LARS[,c("taxon", "LARS")],
                           by="taxon", all=TRUE)
      colnames(tableS.LARS) <- c("taxon", "indicatee", "stat", "p.value", "A", "B")
      
      rm(tableA.LARS, tableB.LARS)
      }
    }
    
    # wild (either deep-snow or shallow-snow)
    if(nrow(tableS.wild) > 0){
      tableS.wild$indicatee <- "wild"
      tableS.wild$taxon <- rownames(tableS.wild)
      
      tableA.wild <- as.data.frame(subset(tableA, rownames(tableA) %in% rownames(tableS.wild)))
      tableA.wild$taxon <- rownames(tableA.wild)
      
      tableS.wild <- merge(tableS.wild[,c("indicatee", "stat", "p.value", "taxon")],
                           as.data.frame(tableA.wild[,c("taxon", "deep-snow+shallow-snow")]),
                           by="taxon", all=TRUE)
      colnames(tableS.wild) <- c("taxon", "indicatee", "stat", "p.value", "A")
      
      tableB.wild <- as.data.frame(subset(tableB, rownames(tableB) %in% tableS.wild$taxon))
      tableB.wild$taxon <- rownames(tableB.wild)
      
      tableS.wild <- merge(tableS.wild,
                           tableB.wild[,c("taxon", "deep-snow+shallow-snow")],
                           by="taxon", all=TRUE)
      colnames(tableS.wild) <- c("taxon", "indicatee", "stat", "p.value", "A", "B")
      
      rm(tableA.wild, tableB.wild)
    }
    
    # both deep snow (either wild deep-snow caribou or Revelstoke pen)
    if(nrow(tableS.deepsnow) > 0){
      tableS.deepsnow$indicatee <- "bothdeepsnow"
      tableS.deepsnow$taxon <- rownames(tableS.deepsnow)
      
      tableA.deepsnow <- as.data.frame(subset(tableA, rownames(tableA) %in% rownames(tableS.deepsnow)))
      tableA.deepsnow$taxon <- rownames(tableA.deepsnow)
      
      tableS.deepsnow <- merge(tableS.deepsnow[,c("indicatee", "stat", "p.value", "taxon")],
                               as.data.frame(tableA.deepsnow[,c("taxon", "deep-snow+Revelstoke pen")]),
                               by="taxon", all=TRUE)
      colnames(tableS.deepsnow) <- c("taxon", "indicatee", "stat", "p.value", "A")
      
      tableB.deepsnow <- as.data.frame(subset(tableB, rownames(tableB) %in% tableS.deepsnow$taxon))
      tableB.deepsnow$taxon <- rownames(tableB.deepsnow)
      
      tableS.deepsnow <- merge(tableS.deepsnow,
                               tableB.deepsnow[,c("taxon", "deep-snow+Revelstoke pen")],
                               by="taxon", all=TRUE)
      colnames(tableS.deepsnow) <- c("taxon", "indicatee", "stat", "p.value", "A", "B")
      
      rm(tableA.deepsnow, tableB.deepsnow)
    }
    rm(tableA, tableB)
    
    # Stick all these tables together ####
    if(p < 7){
      tableS <- rbind(tableS.SM, tableS.CM,
                      tableS.RP, tableS.LARS,
                      tableS.wild, tableS.deepsnow)
    }
    if(p==7){
      tableS <- rbind(tableS.SM, tableS.CM,
        tableS.RP, tableS.wild, tableS.deepsnow)
    }
    
    rm(tableS.SM, tableS.CM, tableS.RP, tableS.LARS, tableS.wild, tableS.deepsnow)
    
    indicator.species[[p]] <- tableS
    rm(tableS, temp.otu, grps)
}

names(indicator.species) <- names(pseq.clr)

## [5] Create a single data frame that shows Gini scores and indicator results for each genus. ####
indicator.Gini.synthesis <- list()
for(p in c(1:7)){
  temp <- rf.results.all[[p]]$Gini_scores
  rownames(temp) <- temp$Taxon
  temp$Taxon <- NULL
  
  rownames(indicator.species[[p]]) <- indicator.species[[p]]$taxon
  temp <- merge(indicator.species[[p]], temp, by=0, all=TRUE)
  rownames(temp) <- temp$Row.names
  temp$Row.names <- NULL
  
  indicator.Gini.synthesis[[p]] <- temp
  rm(temp)
}

####### CORRELATION-BASED NETWORK ANALYSIS ####
## [1] Prepare the data sets. ####
# Create a list of genus-level phyloseq objects.
supertax <- list(
  X16S = tax_glom(pseq.clr.sub[[1]], taxrank="Genus"),
  X18S = tax_glom(pseq.clr.sub[[2]], taxrank="Genus"),
  FITS = tax_glom(pseq.clr.sub[[3]], taxrank="Genus"),
  PITS = tax_glom(pseq.clr.sub[[4]], taxrank="Genus"),
  protist = tax_glom(pseq.clr.sub$protists, taxrank="Genus"),
  algae = tax_glom(pseq.clr.sub$algae, taxrank="Genus"),
  lichen = tax_glom(pseq.clr.sub$lichens, taxrank="Genus")
)

groups <- c("deep-snow", "shallow-snow", "Revelstoke pen", "LARS")

# Create a separate file where each caribou population has its own phyloseq object for each amplicon.
# supertax.group -> four populations -> eight amplicons
supertax.group <- list()
for(i in c(1:4)){
  supertax.group[[i]] <- list()
  
  if(i < 4){
    for(p in c(1:7)){
      supertax.group[[i]][[p]] <- subset_samples(supertax[[p]], Group==groups[i])
      supertax.group[[i]][[p]] <- prune_taxa(taxa_sums(supertax.group[[i]][[p]]) > 0, supertax.group[[i]][[p]])
    }
  }
  if(i==4){
    for(p in c(1:6)){
      supertax.group[[i]][[p]] <- subset_samples(supertax[[p]], Group==groups[i])
      supertax.group[[i]][[p]] <- prune_taxa(taxa_sums(supertax.group[[i]][[p]]) > 0, supertax.group[[i]][[p]])
    }
  }
  
  if(i < 4){names(supertax.group[[i]]) <- names(pseq.clr)}
  if(i==4){names(supertax.group[[i]]) <- names(pseq.clr)[c(1:6)]}
}
names(supertax.group) <- c("deep-snow", "shallow-snow", "Revelstoke pen", "LARS")

# Convert taxon names to genus names and convert all objects to relative abundances.
for(p in c(1:length(supertax))){
  supertax[[p]] <- transform_sample_counts(supertax[[p]], function(x) 100*x/sum(x))
  taxa_names(supertax[[p]]) <- as.data.frame(cbind(tax_table(supertax[[p]])))$Genus
  
  if(p < 7){
    for(i in c(1:4)){
      supertax.group[[i]][[p]] <- transform_sample_counts(supertax.group[[i]][[p]], function(x) 100*x/sum(x))
      taxa_names(supertax.group[[i]][[p]]) <- as.data.frame(cbind(tax_table(supertax.group[[i]][[p]])))$Genus
    }
  }
  if(p==7){
    for(i in c(1:3)){
      supertax.group[[i]][[p]] <- transform_sample_counts(supertax.group[[i]][[p]], function(x) 100*x/sum(x))
      taxa_names(supertax.group[[i]][[p]]) <- as.data.frame(cbind(tax_table(supertax.group[[i]][[p]])))$Genus
    }
  }
}

# Subset to only taxa present with >0.5% abundance in each amplicon.
for(p in c(1:7)){
  taxa_list <- subset(prev.data[[p]]$Genus, Abundance > 0.5)$Taxon
  supertax[[p]] <- subset_taxa(supertax[[p]], 
                               taxa_names(supertax[[p]]) %in% taxa_list)
  
  for(i in c(1:4)){
    if(p==7 & i==4){print("Can't do it")} else{
      supertax.group[[i]][[p]] <- subset_taxa(supertax.group[[i]][[p]], 
                                              taxa_names(supertax.group[[i]][[p]]) %in% taxa_list)
    }
  }
}

## [2] Run the cross-amplicon correlations. ####
groups <- c("BAC", "EUK", "FUNG", "PLANT", "PROT", "ALG", "LICH")

# BAC=1, PROT=5, FUNG=3, PLANT=4, ALGAE=6, LICHEN=7

allway.corr <- list()

# Define 14 different pairwise comparisons based on the 6 data sets of interest.
for(j in c(1:14)){
  if(j==1){a=1
  b=5} # BAC_PROT
  if(j==2){a=1
  b=3} # BAC_FUNG
  if(j==3){a=1
  b=4} # BAC_PLANT
  if(j==4){a=1
  b=6} # BAC_ALG
  if(j==5){a=1
  b=7} # BAC_LICH
  if(j==6){a=5
  b=3} # PROT_FUNG
  if(j==7){a=5
  b=4} # PROT_PLANT
  if(j==8){a=5
  b=6} # PROT_ALG
  if(j==9){a=5
  b=7} # PROT_LICH
  if(j==10){a=3
  b=4} # FUNG_PLANT
  if(j==11){a=3
  b=6} # FUNG_ALG
  if(j==12){a=4
  b=6} # PLANT_ALG
  if(j==13){a=4
  b=7} # PLANT_LICH
  if(j==14){a=6
  b=7} # ALG_LICH
  
  # Create a single data frame with the relative abundances of all taxa in both comparison groups.
  df <- as.data.frame(otu_table(supertax[[a]]))
  df <- merge(df, as.data.frame(otu_table(supertax[[b]])), by=0, all=TRUE)
  rownames(df) <- df$Row.names
  df$Row.names <- NULL
  
  # Run taxa-by-taxa correlations.
  supertax.corr <- Hmisc::rcorr(as.matrix(df), type="spearman")
  
  # For the 'r' and 'p' values, subset to only cross-amplicon comparisons.
  for(k in c(1:length(supertax.corr))){
    supertax.corr[[k]] <- subset(supertax.corr[[k]], rownames(supertax.corr[[k]]) %in% taxa_names(supertax[[a]]))
    supertax.corr[[k]] <- as.data.frame(t(supertax.corr[[k]]))
    supertax.corr[[k]] <- subset(supertax.corr[[k]], rownames(supertax.corr[[k]]) %in% taxa_names(supertax[[b]]))
    supertax.corr[[k]] <- as.data.frame(t(supertax.corr[[k]]))
    supertax.corr[[k]]$TaxName <- rownames(supertax.corr[[k]])
    supertax.corr[[k]] <- reshape2::melt(supertax.corr[[k]])
  }
  
  colnames(supertax.corr$r) <- c("Tax1", "Tax2", "R")
  colnames(supertax.corr$P) <- c("Tax1", "Tax2", "p")
  
  # Merge the 'r' and 'p' values into a single data frame.
  supertax.corr <- merge(supertax.corr$r, supertax.corr$P, by=c("Tax1", "Tax2"), all=TRUE)
  supertax.corr <- subset(supertax.corr, is.na(R)=="FALSE")
  
  # Adjust the p-values for multiple comparisons.
  supertax.corr$p.adj <- p.adjust(supertax.corr$p, method = "fdr")
  
  # Add data on the two amplicon data sets being compared.
  supertax.corr$Tax1Group <- groups[a]
  supertax.corr$Tax2Group <- groups[b]
  
  # Save data.
  allway.corr[[j]] <- supertax.corr
  rm(supertax.corr, df)
}

## [3] Run same-domain correlations (BAC-BAC, FUNG-FUNG, PLANT-PLANT, PROT-PROT, ALG-ALG, LICH-LICH) ####
for(j in c(15:20)){
  if(j==15){a=1} # BAC_BAC
  if(j==16){a=3} # FUNG-FUNG
  if(j==17){a=4} # PLANT_PLANT
  if(j==18){a=5} # PROT_PROT
  if(j==19){a=6} # ALG_ALG
  if(j==20){a=7} # LICH-LICH
  
  df <- as.data.frame(otu_table(supertax[[a]]))
  supertax.corr <- Hmisc::rcorr(as.matrix(df), type="spearman")
  
  for(k in c(1:length(supertax.corr))){
    supertax.corr[[k]] <- as.data.frame(supertax.corr[[k]])
    supertax.corr[[k]][upper.tri(supertax.corr[[k]])] <- NA
    supertax.corr[[k]]$TaxName <- rownames(supertax.corr[[k]])
    supertax.corr[[k]] <- reshape2::melt(supertax.corr[[k]])
  }
  
  colnames(supertax.corr$r) <- c("Tax1", "Tax2", "R")
  colnames(supertax.corr$P) <- c("Tax1", "Tax2", "p")
  
  supertax.corr <- merge(supertax.corr$r, supertax.corr$P, by=c("Tax1", "Tax2"), all=TRUE)
  supertax.corr <- subset(supertax.corr, is.na(p)=="FALSE")
  supertax.corr$p.adj <- p.adjust(supertax.corr$p, method = "fdr")
  
  supertax.corr$Tax1Group <- groups[a]
  supertax.corr$Tax2Group <- groups[a]
  
  allway.corr[[j]] <- supertax.corr
  
  rm(supertax.corr, df)
}

## [4] Isolate positive correlations for network analysis in CytoScape ####
# Subset to only positive correlations and export for analysis in CytoScape.
allway.corr.signif <- allway.corr
for(j in c(1:20)){
  allway.corr.signif[[j]] <- subset(allway.corr.signif[[j]], R > 0)
}

allway.corr.signif2 <- dplyr::bind_rows(allway.corr.signif[c(1,4,5,8,9,14,15,18,19,20)])
allway.corr.signif2 <- subset(allway.corr.signif2, p.adj < 0.01 & R > 0.6)

write.csv(allway.corr.signif2, "~/network_import.csv")

# Calculate the percentage of positive correlations as a function of the total number of correlations.
compares <- c("BAC_PROT", "BAC_FUNG", "BAC_PLANT", "BAC_ALG", "BAC_LICH", 
              "PROT_FUNG", "PROT_PLANT", "PROT_ALG", "PROT_LICH", "FUNG_PLANT", 
              "FUNG_ALG", "PLANT_ALG", "PLANT_LICH", "ALG_LICH", "BAC_BAC", "FUNG_FUNG",
              "PLANT_PLANT", "PROT_PROT", "ALG_ALG", "LICH_LICH")

allway.corr.group.binary <- allway.corr.group
for(j in c(1:20)){
  allway.corr.group.binary[[j]]$R[allway.corr.group.binary[[j]]$R > 0] <- 1
  allway.corr.group.binary[[j]]$R[allway.corr.group.binary[[j]]$R < 0] <- 0
  
  allway.corr.group.binary[[j]] <- allway.corr.group.binary[[j]] %>%
    dplyr::group_by(StudyGroup) %>%
    dplyr::summarise(total = dplyr::n(),
                     positives = sum(R))
  allway.corr.group.binary[[j]]$compare <- compares[j]
}

allway.corr.group.binary <- as.data.frame(dplyr::bind_rows(allway.corr.group.binary))
allway.corr.group.binary$percent <- 100*allway.corr.group.binary$positives/allway.corr.group.binary$total

####### MANTEL TESTS #####
## [1] Aitchison distance ####
mantel.test.results <- data.frame(matrix(ncol=5))

mantel <- list()

for(i in c(1:7)){
  for(j in c(1:7)){
    if(i != j){
      mantel[[i]] <- pseq.clr.sub[[i]]
      
      if(i==4){mantel[[i]] <- subset_samples(mantel[[i]], Group !="deep-snow")}
      if(i==7){mantel[[i]] <- subset_samples(mantel[[i]], Group !="LARS")}
      
      mantel[[j]] <- pseq.clr.sub[[j]]
      
      if(j==4){mantel[[j]] <- subset_samples(mantel[[j]], Group !="deep-snow")}
      if(j==7){mantel[[j]] <- subset_samples(mantel[[j]], Group !="LARS")}
      
      # Subset so the two frames have the same number of samples.
      if(nsamples(mantel[[j]]) > nsamples(mantel[[i]])){
        mantel[[j]] <- subset_samples(mantel[[j]], sample_names(mantel[[j]]) %in%
                                        sample_names(mantel[[i]]))}
      if(nsamples(mantel[[j]]) < nsamples(mantel[[i]])){
        mantel[[i]] <- subset_samples(mantel[[i]], sample_names(mantel[[i]]) %in%
                                        sample_names(mantel[[j]]))}
      if(nsamples(mantel[[j]]) > nsamples(mantel[[i]])){
        mantel[[j]] <- subset_samples(mantel[[j]], sample_names(mantel[[j]]) %in%
                                        sample_names(mantel[[i]]))}
      
      number <- nsamples(mantel[[i]])
      
      mantel[[i]] <- vegdist(as.data.frame(otu_table(microbiome::transform(mantel[[i]], "clr"))),
                                  method="euclidean")
      mantel[[j]] <- vegdist(as.data.frame(otu_table(microbiome::transform(mantel[[j]], "clr"))),
                             method="euclidean")
      
      x <- vegan::mantel(mantel[[i]], 
                         mantel[[j]],
                         method="spearman")
      
      vector <- c(names(pseq.clr.sub)[i], names(pseq.clr.sub)[j], number, x$statistic, x$signif)
      mantel.test.results <- rbind(mantel.test.results, vector)
    }
  }
}
 
mantel.test.results$dup <- duplicated(mantel.test.results$X4)
mantel.test.results <- subset(mantel.test.results, dup=="FALSE")
mantel.test.results$X4 <- as.numeric(as.character(mantel.test.results$X4))

write.csv(mantel.test.results, "~/mantel.csv")

## [2] Other distance metrics ####
mantel.test.results <- data.frame(matrix(ncol=11))
mantel <- list()

for(i in c(1:7)){
  for(j in c(1:7)){
    if(i != j){
      mantel[[i]] <- pseq.clr.sub[[i]]
      
      if(i==4){mantel[[i]] <- subset_samples(mantel[[i]], Group !="deep-snow")}
      if(i==7){mantel[[i]] <- subset_samples(mantel[[i]], Group !="LARS")}
      
      mantel[[j]] <- pseq.clr.sub[[j]]
      
      if(j==4){mantel[[j]] <- subset_samples(mantel[[j]], Group !="deep-snow")}
      if(j==7){mantel[[j]] <- subset_samples(mantel[[j]], Group !="LARS")}
      
      # Subset so the two frames have the same number of samples.
      if(nsamples(mantel[[j]]) > nsamples(mantel[[i]])){
        mantel[[j]] <- subset_samples(mantel[[j]], sample_names(mantel[[j]]) %in%
                                        sample_names(mantel[[i]]))}
      if(nsamples(mantel[[j]]) < nsamples(mantel[[i]])){
        mantel[[i]] <- subset_samples(mantel[[i]], sample_names(mantel[[i]]) %in%
                                        sample_names(mantel[[j]]))}
      if(nsamples(mantel[[j]]) > nsamples(mantel[[i]])){
        mantel[[j]] <- subset_samples(mantel[[j]], sample_names(mantel[[j]]) %in%
                                        sample_names(mantel[[i]]))}
      
      number <- nsamples(mantel[[i]])
      
      mantel[[i]] <- transform_sample_counts(mantel[[i]], function(x) 100*x/sum(x))
      mantel[[j]] <- transform_sample_counts(mantel[[j]], function(x) 100*x/sum(x))
      
      dist.list <- list()
      dist.list[[i]] <- list()
      dist.list[[i]][[1]] <- vegdist(as.data.frame(otu_table(mantel[[i]])), method="bray")
      dist.list[[i]][[2]] <- vegdist(as.data.frame(otu_table(mantel[[i]])), method="jaccard")
      dist.list[[i]][[3]] <- UniFrac(mantel[[i]], weighted=TRUE, normalized=TRUE)
      dist.list[[i]][[4]] <- UniFrac(mantel[[i]], weighted=FALSE, normalized=TRUE)
      
      dist.list[[j]] <- list()
      dist.list[[j]][[1]] <- vegdist(as.data.frame(otu_table(mantel[[j]])), method="bray")
      dist.list[[j]][[2]] <- vegdist(as.data.frame(otu_table(mantel[[j]])), method="jaccard")
      dist.list[[j]][[3]] <- UniFrac(mantel[[j]], weighted=TRUE, normalized=TRUE)
      dist.list[[j]][[4]] <- UniFrac(mantel[[j]], weighted=FALSE, normalized=TRUE)
      
      x.bc <- vegan::mantel(dist.list[[i]][[1]], dist.list[[j]][[1]],
                            method="spearman")
      x.jc <- vegan::mantel(dist.list[[i]][[2]], dist.list[[j]][[2]],
                            method="spearman")
      x.wUF <- vegan::mantel(dist.list[[i]][[3]], dist.list[[j]][[3]],
                            method="spearman")
      x.uwUF <- vegan::mantel(dist.list[[i]][[4]], dist.list[[j]][[4]],
                            method="spearman")
      
      vector <- c(names(pseq.clr.sub)[i], names(pseq.clr.sub)[j], number, 
                  x.bc$statistic, x.bc$signif,
                  x.jc$statistic, x.jc$signif,
                  x.wUF$statistic, x.wUF$signif,
                  x.uwUF$statistic, x.uwUF$signif)
      mantel.test.results <- rbind(mantel.test.results, vector)
      
      rm(vector, number, x.bc, x.jc, x.wUF, x.uwUF)
    }
  }
}

mantel.test.results$dup <- duplicated(mantel.test.results$X4)
mantel.test.results <- subset(mantel.test.results, dup=="FALSE")
mantel.test.results$X4 <- as.numeric(as.character(mantel.test.results$X4))
mantel.test.results$X6 <- as.numeric(as.character(mantel.test.results$X6))
mantel.test.results$X8 <- as.numeric(as.character(mantel.test.results$X8))
mantel.test.results$X10 <- as.numeric(as.character(mantel.test.results$X10))

write.csv(mantel.test.results, "~/mantel.csv")

#### MURIBACULACEAE GENOME ANALYSIS ####
# Download genomes from NCBI using scripts in the terminal.
# Annotate genomes using dbCAN using scripts in the terminal.

## [1] Look at a phylogenetic tree of ASV sequences. ####
temp <- subset_taxa(pseq.clr$X16S, Family=="Muribaculaceae")
temp2 <- as.data.frame(cbind(tax_table(temp)))
temp2$name <- paste0(temp2$Genus, "_", rownames(temp2))
taxa_names(temp) <- temp2$name
rm(temp2)

db_out <- data.frame(
  ids = taxa_names(temp),
  seqs = refseq(temp)
)

fasta <- ShortRead(sread = DNAStringSet(db_out$seqs),
                   id = BStringSet(db_out$ids))
writeFasta(fasta, "asv_concantenate_seqs.fasta")

# Create a tree using MAFFT.
# Compare to the tree generated by phangorn:
tree <- phy_tree(temp)
ape::write.tree(tree, "phyloseq_tree.tree")

## [2] Import Muribaculaceae genome annotations into R for visualization. ####
# For each species, read the annotations file.
annotations <- list(
  # Paramuribaculum intestinale
  p.intestinale = read.csv("GCF_003024815.1_P_intestinale/fna_output/overview.txt", sep="\t"),
  # Muribaculum intestinale
  m.intestinale = read.csv("GCF_001688845.2_M_intestinale/fna_output/overview.txt", sep="\t"),
  # Duncaniella muris
  d.muris = read.csv("GCF_003024805.1_D_muris/fna_output/overview.txt", sep="\t"),
  # Muribaculum caecicola
  m.caecicola = read.csv("GCF_004801635.1_M_caecicola/fna_output/overview.txt", sep="\t"),
  # Muribaculum gordoncarteri
  m.gordoncarteri = read.csv("GCF_004803695.1_M_gordoncarteri//fna_output/overview.txt", sep="\t"),
  # Duncaniella dubosii
  d.dubosii = read.csv("GCF_004803915.1_D_dubosii/fna_output/overview.txt", sep="\t"),
  # Duncaniella murioclitica
  d.muricolitica = read.csv("GCF_910574735.1_D_muricolitica/fna_output/overview.txt", sep="\t"),
  # Duncaniella freteri
  d.freteri = read.csv("GCF_004766125.1_D_freteri/fna_output/overview.txt", sep="\t"),
  # Phocaeicola vulgatus
  p.vulgatus = read.csv("GCF_020885855.1_P_vulgatus/fna_output/overview.txt", sep="\t"),
  # Phocaeicola dorei
  p.dorei = read.csv("GCF_902387545.1_P_dorei/fna_output/overview.txt", sep="\t"),
  # Mediterranea massiliensis
  med.massiliensis = read.csv("GCF_900128475.1_Med_massiliensis/fna_output/overview.txt", sep="\t")
)

# Clean the annotations file.
for(i in c(1:length(annotations))){
  
  temp <- annotations[[i]] %>%
    tidyr::separate(., "HMMER", into=c("HMMER_1", "HMMER_2"), sep="\\+", remove=TRUE) %>%
    tidyr::separate(., "dbCAN_sub", into=c("dbCAN_1", "dbCAN_2"), sep="\\+", remove=TRUE) %>%
    tidyr::separate(., "DIAMOND", into=c("DIAMOND_1", "DIAMOND_2"), sep="\\+", remove=TRUE) %>%
    dplyr::select(-EC.) %>%
    subset(., X.ofTools > 1) %>%
    tidyr::separate(., "HMMER_1", into=c("HMMER_1", "junk1"), sep="\\(", remove=TRUE) %>%
    tidyr::separate(., "HMMER_2", into=c("HMMER_2", "junk2"), sep="\\(", remove=TRUE) %>%
    tidyr::separate(., "HMMER_1", into=c("HMMER_1", "junk3"), sep="_", remove=TRUE) %>%
    tidyr::separate(., "HMMER_2", into=c("HMMER_2", "junk4"), sep="_", remove=TRUE) %>%
    tidyr::separate(., "dbCAN_1", into=c("dbCAN_1", "junk5"), sep="_", remove=TRUE) %>%
    tidyr::separate(., "dbCAN_2", into=c("dbCAN_2", "junk6"), sep="_", remove=TRUE) %>%
    tidyr::separate(., "DIAMOND_1", into=c("DIAMOND_1", "junk7"), sep="_", remove=TRUE) %>%
    tidyr::separate(., "DIAMOND_2", into=c("DIAMOND_2", "junk8"), sep="_", remove=TRUE) %>%
    dplyr::select(-junk1, -junk2, -junk3, -junk4, -junk5, -junk6, -junk7, -junk8)
  temp[temp=="-"] <- NA
  
  # For genes that are double-annotated, count each annotation.
  temp1 <- subset(temp, is.na(HMMER_2)=="TRUE" & is.na(dbCAN_2)=="TRUE" & is.na(DIAMOND_2)=="TRUE")
  temp2 <- subset(temp, !Gene.ID %in% temp1$Gene.ID)
  temp1$Final <- temp1$HMMER_1
  temp1$Final[is.na(temp1$Final)=="TRUE"] <- temp1$dbCAN_1[is.na(temp1$Final)=="TRUE"]
  temp1[,c("HMMER_1", "HMMER_2", "dbCAN_1", "dbCAN_2", "DIAMOND_1", "DIAMOND_2", "X.ofTools")] <- NULL
  
  temp2 <- data.frame(
    Gene.ID = rep(temp2$Gene.ID, 6),
    Final = c(temp2[,2], temp2[,3], temp2[,4], temp2[,5], temp2[,6], temp2[,7])
  )
  temp2 <- subset(temp2, is.na(Final)=="FALSE")
  temp2 <- temp2 %>% 
    dplyr::group_by(Gene.ID, Final) %>%
    dplyr::count()
  temp2$n <- NULL
  
  # Put data frames back together.
  annotations[[p]] <- rbind(temp1, temp2)
  annotations[[p]] <- annotations[[p]] %>% dplyr::group_by(Final) %>% dplyr::count()
  
  rm(temp, temp1, temp2)
}

## [3] Summarize genome annotations ####
for(i in c(1:length(annotations))){
  colnames(annotations[[i]]) <- c("final", names(annotations)[i])
}

annotations_tab <- merge(annotations[[1]], annotations[[2]], by="final", all=TRUE)
annotations_tab <- merge(annotations_tab, annotations[[3]], by="final", all=TRUE)
annotations_tab <- merge(annotations_tab, annotations[[4]], by="final", all=TRUE)
annotations_tab <- merge(annotations_tab, annotations[[5]], by="final", all=TRUE)
annotations_tab <- merge(annotations_tab, annotations[[6]], by="final", all=TRUE)
annotations_tab <- merge(annotations_tab, annotations[[7]], by="final", all=TRUE)
annotations_tab <- merge(annotations_tab, annotations[[8]], by="final", all=TRUE)
annotations_tab <- merge(annotations_tab, annotations[[9]], by="final", all=TRUE)
annotations_tab <- merge(annotations_tab, annotations[[10]], by="final", all=TRUE)
annotations_tab <- merge(annotations_tab, annotations[[11]], by="final", all=TRUE)

annotations_tab <- subset(annotations_tab, final !="NA")

rownames(annotations_tab) <- annotations_tab$final
annotations_tab$final <- NULL

# Look at a hierarchical clustering.
annotations_tab[is.na(annotations_tab)=="TRUE"] <- 0
dist_mat <- vegan::vegdist(as.data.frame(t(annotations_tab)), method="bray")
clust <- hclust(dist_mat)
plot(clust)

plot_tab <- annotations_tab
plot_tab$cluster <- rownames(plot_tab)
plot_tab <- reshape2::melt(plot_tab)
plot_tab$facet <- substr(plot_tab$cluster, start = 1, stop = 2)
plot_tab$facet[plot_tab$facet=="CB"] <- "CBM"

plot_tab$number <- substr(plot_tab$cluster, start = 3, stop = 6)
plot_tab$number[plot_tab$facet=="CBM"] <- substr(plot_tab$cluster[plot_tab$facet=="CBM"], start = 4, stop = 6)

plot_tab$number <- as.numeric(as.character(plot_tab$number))
plot_tab <- plot_tab[order(plot_tab$facet, plot_tab$number), ]

plot_tab$cluster <- forcats::fct_inorder(plot_tab$cluster)
