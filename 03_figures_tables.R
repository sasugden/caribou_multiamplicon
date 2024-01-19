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

#### *FIG. 2: LICHEN PANELS ####
# (a) Alpha diversity scatter plot ####
alpha.means <- list()
for(p in c(1:7)){
  alpha.means[[p]] <- sample_data.mb[[p]] %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise(m.Observed_Extrap = mean(Observed_Extrap, na.rm=TRUE),
                     sd.Observed_Extrap = sd(Observed_Extrap, na.rm=TRUE),
                     m.Shannon_Extrap = mean(Shannon_Extrap, na.rm=TRUE),
                     sd.Shannon_Extrap = sd(Shannon_Extrap, na.rm=TRUE))
  alpha.means[[p]]$Group <- factor(alpha.means[[p]]$Group,
                                   levels=c("deep-snow", "shallow-snow", "Revelstoke pen", "LARS"))
  
}

p1 <- ggplot() +
  geom_errorbar(data = alpha.means[[7]], aes(x=m.Observed_Extrap,
                                             ymin=m.Shannon_Extrap - sd.Shannon_Extrap,
                                             ymax=m.Shannon_Extrap + sd.Shannon_Extrap),
                color="grey", width=0) +
  geom_errorbarh(data = alpha.means[[7]], aes(y=m.Shannon_Extrap,
                                              xmin=m.Observed_Extrap - sd.Observed_Extrap,
                                              xmax=m.Observed_Extrap + sd.Observed_Extrap),
                 color="grey", height=0) +
  geom_point(data = sample_data.mb[[7]], aes(x=Observed_Extrap, y=Shannon_Extrap, color=Group), size=0.5) +
  geom_point(data = alpha.means[[7]], aes(x=m.Observed_Extrap, y=m.Shannon_Extrap,
                                          color=Group, shape=Group), size=3.5) +
  scale_shape_manual(values=c(16,17,15,18)) +
  theme_bw() + plot_theme +
  guides(shape=guide_legend(ncol=1, byrow=TRUE),
         color=guide_legend(ncol=1, byrow=TRUE, override.aes = list(size = 2, shape=c(16,17,15,18)))) +
  labs(x="ASV richness", y="Shannon diversity", tag="a") +
  scale_color_manual(values=c("chocolate1", "chartreuse3",
                              "darkorchid1", "gray62")) +
  theme(legend.position=c(0.76,0.15), legend.title=element_blank(),
        legend.key.height=unit(0.75,"line"),
        legend.key.width=unit(0.75,"line"),
        legend.text=element_text(size=7),
        legend.spacing.y=unit(0.002,'cm'))

# (b) Ordination ####
temp <- bdiv.ordinations.all[[7]][["points"]][[3]]
temp$Herd <- factor(temp$Herd, levels=c("Hart Ranges", "North Cariboo", "Barkerville",
                                        "Wells Grey North", "Columbia North", "Central Selkirks",
                                        "Tonquin", "Brazeau", "Revelstoke pen"))

temp$GroupHerd <- paste0(temp$Group, temp$Herd)
temp$GroupHerd <- factor(temp$GroupHerd,
                         levels=c("deep-snowHart Ranges", "deep-snowNorth Cariboo",
                                  "deep-snowBarkerville", "deep-snowWells Grey North",
                                  "deep-snowColumbia North", "deep-snowCentral Selkirks",
                                  "shallow-snowTonquin", "shallow-snowBrazeau", "Revelstoke penRevelstoke pen"))

p2 <- ggplot(data =  temp, aes(x=PC1, y=PC2, color=GroupHerd, shape=GroupHerd)) + 
  geom_point(size=1.5) +
  stat_ellipse(aes(group=Group), level=0.95, linewidth=0.5, color="darkgrey") +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=1)) +
  labs(x="PC1 (32.0%)", y="PC2 (15.9%)", tag="b") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.text=element_text(size=7),
        legend.position="right",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0),
        legend.box="horizontal",
        legend.key.height=unit(0.75,"line"),
        legend.key.width=unit(0.75,"line"),
        legend.spacing.y=unit(0.01,"cm")) +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)], 
                                             brewer.pal(4, "Paired")[c(3:4)], "darkorchid1"),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen")) +
  scale_shape_manual(name="Guide1", values=c(16,16,16,16,16,16,17,17,15),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen"))

# (c) Heat map ####
# Identify taxa that are significantly differentially abundant among study populations or herds.
a.df.lichens <- subset(diff.abund.GROUP[["lichens"]][["Genus"]], glm.eBH < 0.05)
a.df.lichens.hrd <- subset(diff.abund.HERD[["lichens"]][["Genus"]], glm.eBH < 0.05)

df.lichens <- subset(diff.abund.GROUP[["lichens"]][["Genus"]], glm.eBH < 0.05 | Genus %in% a.df.lichens.hrd$Genus)
df.lichens.hrd <- subset(diff.abund.HERD[["lichens"]][["Genus"]], glm.eBH < 0.05 | Genus %in% a.df.lichens$Genus)

rm(a.df.lichens, a.df.lichens.hrd, df.lichens.hrd)

# Extract abundance information for these taxa.
df.lichens2 <- subset(prev.data$lichens$Genus, Taxon %in% df.lichens$Genus)
df.lichens2$Genus <- df.lichens2$Taxon
df.lichens2$Taxon <- NULL
df.lichens2 <- subset(df.lichens2, Abundance > 0.1)

df.lichens2[,c("Prevalence", "Abundance", "LARS")] <- NULL
rownames(df.lichens2) <- df.lichens2$Genus
df.lichens2$Genus <- NULL
df.lichens2 <- as.data.frame(t(df.lichens2))
df.lichens2$Group <- rownames(df.lichens2)
df.lichens2 <- reshape2::melt(df.lichens2)

# Prepare faceting variables for plotting.
df.lichens2$Facet <- rep("population")
df.lichens2$Facet[df.lichens2$Group %in% c("Barkerville", "Central Selkirks", "Hart Ranges",
                                           "North Cariboo", "Columbia North", "Wells Grey North")] <- "deep-snow"
df.lichens2$Facet[df.lichens2$Group %in% c("Brazeau", "Tonquin")] <- "shallow-snow"

df.lichens2$Group <- factor(df.lichens2$Group, levels=c(
  "deep-snow", "shallow-snow", "Revelstoke pen",
  "Hart Ranges", "North Cariboo", "Barkerville", "Wells Grey North", "Columbia North", "Central Selkirks",
  "Tonquin", "Brazeau"
))

df.lichens2$Facet <- factor(df.lichens2$Facet, levels=c(
  "population", "deep-snow", "shallow-snow"
))

# Replace zero abundances with '-' and replace abundances >50 with 50.
df.lichens2$rounded <- as.numeric(df.lichens2$value)
df.lichens2$rounded[df.lichens2$rounded==0] <- "-"
df.lichens2$rounded[df.lichens2$rounded != "-"] <- NA

df.lichens2$mean_color <- as.numeric(as.character(df.lichens2$value))
df.lichens2$mean_color[df.lichens2$mean_color > 50] <- 50

# Add full taxonomy information.
df.lichens2$Genus <- df.lichens2$variable
df.lichens2$variable <- NULL
df.lichens2 <- merge(df.lichens2,
                     df.lichens[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")],
                     by="Genus", all.x=TRUE, all.y=FALSE)

df.lichens2$Family[df.lichens2$Family !="Parmeliaceae" & df.lichens2$Class=="Lecanoromycetes"] <- "Other Lecanoromycetes"
df.lichens2$Family[df.lichens2$Family !="Parmeliaceae" & df.lichens2$Family !="Other Lecanoromycetes"] <- "Other"
df.lichens2$Family <- factor(df.lichens2$Family, levels=c("Parmeliaceae", "Other Lecanoromycetes", "Other"))

# Add indicator species information.
indic.lichen <- indicator.Gini.synthesis[[7]]
indic.lichen <- subset(indic.lichen, Genus %in% df.lichens2$Genus)

indic.lichen2 <- subset(indic.lichen, indicatee=="bothdeepsnow")
indic.lichen2$indicatee <- rep("deep-snow")
indic.lichen3 <- subset(indic.lichen, indicatee=="bothdeepsnow")
indic.lichen3$indicatee <- rep("Revelstoke pen")
indic.lichen <- subset(indic.lichen, indicatee !="bothdeepsnow" | is.na(indicatee=="TRUE"))

indic.lichen <- rbind(indic.lichen, indic.lichen2, indic.lichen3)
rm(indic.lichen2, indic.lichen3)

indic.lichen$symbol <- rep("B")
indic.lichen$symbol[indic.lichen$p.value < 0.05] <- "A"
indic.lichen$Group <- indic.lichen$indicatee

df.lichens2 <- merge(df.lichens2, indic.lichen[,c("Genus","A","B","symbol","Group")],
                     by=c("Group","Genus"), all=TRUE)

df.lichens2 <- subset(df.lichens2, is.na(Group)=="FALSE")
df.lichens2$symbol[is.na(df.lichens2$symbol)=="TRUE"] <- "B"

# Add line breaks for plotting.
levels(df.lichens2$Facet)[levels(df.lichens2$Facet)=="shallow-snow"] <- "shallow-\nsnow"

# Plot
p3 <- ggplot(df.lichens2, aes(x=Genus, y=Group, fill=sqrt(mean_color))) +
  geom_tile() +
  geom_text(aes(label=rounded), size=2.25) +
  geom_point(data=subset(df.lichens2, symbol=="A"), 
             aes(x=Genus, y=Group)) +
  facet_grid(Facet~Family, scales="free", space="free") +
  scale_fill_gradient(low="white", high="firebrick3") +
  scale_color_manual(values=c("black", "darkgrey")) +
  theme_bw() + plot_theme + 
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0), limits=rev) +
  theme(axis.text.x = element_text(angle=90, vjust=0.25, hjust=1, size=7),
        strip.text = element_text(size=7),
        axis.text.y = element_text(size=7),
        legend.position="none") +
  labs(tag="c", x="Genus", y="Population / Herd\n")

legend <- cowplot::get_legend(ggplot(df.lichens2, aes(x=Genus, y=Group, fill=sqrt(mean_color))) +
                                geom_tile() +
                                facet_grid(Facet~Family, scales="free", space="free") +
                                scale_fill_gradient(low="white", high="firebrick3", 
                                                    breaks=c(0,2.236068,3.162278,4.472136,5.477226,6.324555,7.061068),
                                                    labels=c("0","5","10","20","30","40",">50"),
                                                    name="Relative\nabundance (%)") +
                                scale_color_manual(values=c("black", "darkgrey")) +
                                theme_bw() + plot_theme + 
                                scale_x_discrete(expand=c(0,0)) +
                                scale_y_discrete(expand=c(0,0), limits=rev) +
                                theme(axis.text.x = element_text(angle=90, vjust=0.25, hjust=1, size=8),
                                      strip.text = element_text(size=7),
                                      legend.position="right", legend.text=element_text(size=7),
                                      legend.title=element_text(hjust=0.5, size=7)) +
                                labs(tag="c", x="Genus", y="Population / Herd\n"))

# (-) Assemble ####
grid.arrange(p1, p2, p3, layout_matrix=rbind(c(1,1,1,2,2,2,2),c(3,3,3,3,3,3,3)), heights=c(0.4,0.6))

p1 <- ggplotGrob(p1)
p2 <- ggplotGrob(p2)
p3 <- ggplotGrob(p3)

plot <- arrangeGrob(p1, p2, p3, layout_matrix=rbind(c(1,1,1,2,2,2,2),c(3,3,3,3,3,3,3)), heights=c(0.45,0.55))

ggsave("~/caribou/figures/raw_fig2_lichen_panels.jpg", plot, width=6.5, height=7, units="in", dpi=300)
ggsave("~/caribou/figures/raw_fig2_legend.jpg", legend, width=6.5, height=4, units="in", dpi=300)

rm(p1, p2, p3, plot, temp, df.lichens, df.lichens2, indic.lichen)

#### *FIG. 3: PLANT/ALGAE PANELS ####
# (a) Algae ordination ####
temp <- bdiv.ordinations.all$algae$points[[3]]
temp$Herd <- factor(temp$Herd, levels=c("Hart Ranges", "North Cariboo", "Barkerville",
                                        "Wells Grey North", "Columbia North", "Central Selkirks",
                                        "Tonquin", "Brazeau", "Revelstoke pen", "LARS"))

temp$GroupHerd <- paste0(temp$Group, temp$Herd)
temp$GroupHerd <- factor(temp$GroupHerd,
                         levels=c("deep-snowHart Ranges", "deep-snowNorth Cariboo",
                                  "deep-snowBarkerville", "deep-snowWells Grey North",
                                  "deep-snowColumbia North", "deep-snowCentral Selkirks",
                                  "shallow-snowTonquin", "shallow-snowBrazeau", "Revelstoke penRevelstoke pen",
                                  "LARSLARS"))

summary(bdiv.ordinations.all$algae$models[[3]])

p1 <- ggplot(data =  temp, aes(PC1, PC2)) + 
  geom_point(aes(color = GroupHerd, shape = GroupHerd), size=1.25) +
  stat_ellipse(aes(group=Group), level=0.95, linewidth=0.5, color="darkgrey") +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=2), shape=FALSE) +
  labs(x="PC1 (20.2%)", y="PC2 (19.5%)", tag="a") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.text=element_text(size=7),
        legend.position="none",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0),
        legend.box="horizontal",
        legend.key.height=unit(0.75,"line"),
        legend.key.width=unit(0.75,"line"),
        legend.spacing.y=unit(0.01,"cm")) +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)], 
                                             brewer.pal(4, "Paired")[c(3:4)], 
                                             "darkorchid1", "gray62"),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS")) +
  scale_shape_manual(name="Guide1", values=c(16,16,16,16,16,16,17,17,15,18),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS"))

legend <- cowplot::get_legend(ggplot(data =  temp, aes(PC1, PC2)) + 
                                geom_point(aes(color = GroupHerd, shape = GroupHerd), size=1.25) +
                                stat_ellipse(aes(group=Group), level=0.95, linewidth=0.5, color="darkgrey") +
                                theme_bw() + plot_theme +
                                guides(color=guide_legend(ncol=2)) +
                                labs(x="PC1 (20.2%)", y="PC2 (19.5%)", tag="a") +
                                theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
                                      axis.line = element_line(colour = "black", size=0.5),
                                      legend.title=element_blank(),
                                      legend.text=element_text(size=7),
                                      legend.position="bottom",
                                      legend.box="horizontal",
                                      legend.key.height=unit(0.75,"line"),
                                      legend.key.width=unit(0.75,"line"),
                                      legend.spacing.y=unit(0.01,"cm")) +
                                scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)], 
                                                                           brewer.pal(4, "Paired")[c(3:4)], 
                                                                           "darkorchid1", "gray62"),
                                                   labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                                                            "Wells Grey North", "Columbia North", "Central Selkirks",
                                                            "Tonquin", "Brazeau", "Revelstoke pen","LARS")) +
                                scale_shape_manual(name="Guide1", values=c(16,16,16,16,16,16,17,17,15,18),
                                                   labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                                                            "Wells Grey North", "Columbia North", "Central Selkirks",
                                                            "Tonquin", "Brazeau", "Revelstoke pen","LARS")))

# (b) Algae heat map ####
a.df.algae <- subset(diff.abund.GROUP[["algae"]][["Genus"]], glm.eBH < 0.05)
a.df.algae.hrd <- subset(diff.abund.HERD[["algae"]][["Genus"]], glm.eBH < 0.05)

df.algae <- subset(diff.abund.GROUP[["algae"]][["Genus"]], glm.eBH < 0.05 | Genus %in% a.df.algae.hrd$Genus)
df.algae.hrd <- subset(diff.abund.HERD[["algae"]][["Genus"]], glm.eBH < 0.05 | Genus %in% a.df.algae$Genus)

rm(a.df.algae, a.df.algae.hrd, df.algae.hrd)

# Extract abundance information for these taxa.
df.algae2 <- subset(prev.data$algae$Genus, Taxon %in% df.algae$Genus)
df.algae2$Genus <- df.algae2$Taxon
df.algae2$Taxon <- NULL
df.algae2 <- subset(df.algae2, Abundance > 1 | `deep-snow` > 1)

# List in order of abundance.
df.algae2 <- df.algae2[order(-df.algae2$Abundance), ]
df.algae2$Genus <- forcats::fct_inorder(df.algae2$Genus)

df.algae2[,c("Prevalence", "Abundance")] <- NULL
rownames(df.algae2) <- df.algae2$Genus
df.algae2$Genus <- NULL
df.algae2 <- as.data.frame(t(df.algae2))
df.algae2$Group <- rownames(df.algae2)

# Simplify abundances for plotting.
df.algae2$Uncl_Chlorophyceae <- df.algae2$Chlorophyceae + df.algae2$Uncl_Chlorophyceae
df.algae2$Chlorophyceae <- NULL

df.algae2 <- reshape2::melt(df.algae2)

df.algae2$Facet <- rep("population")
df.algae2$Facet[df.algae2$Group %in% c("Barkerville", "Central Selkirks", "Hart Ranges",
                                       "North Cariboo", "Columbia North", "Wells Grey North")] <- "deep-snow"
df.algae2$Facet[df.algae2$Group %in% c("Brazeau", "Tonquin")] <- "shallow-snow"
df.algae2$Group <- factor(df.algae2$Group, levels=c(
  "deep-snow", "shallow-snow", "Revelstoke pen", "LARS",
  "Hart Ranges", "North Cariboo", "Barkerville", "Wells Grey North", "Columbia North", "Central Selkirks",
  "Tonquin", "Brazeau"
))
df.algae2$Facet <- factor(df.algae2$Facet, levels=c(
  "population", "deep-snow", "shallow-snow"
))

df.algae2$rounded <- as.numeric(df.algae2$value)
df.algae2$rounded[df.algae2$rounded==0] <- "-"
df.algae2$rounded[df.algae2$rounded != "-"] <- NA

df.algae2$mean_color <- as.numeric(as.character(df.algae2$value))
df.algae2$mean_color[df.algae2$mean_color > 50] <- 50

df.algae2$Genus <- df.algae2$variable
df.algae2$variable <- NULL
df.algae2 <- merge(df.algae2,
                   df.algae[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")],
                   by="Genus", all.x=TRUE, all.y=FALSE)

# Add indicator species.
indic.algae <- indicator.Gini.synthesis[[6]]
indic.algae <- subset(indic.algae, Genus %in% df.algae2$Genus)

indic.algae2 <- subset(indic.algae, indicatee=="bothdeepsnow")
indic.algae2$indicatee <- rep("deep-snow")
indic.algae3 <- subset(indic.algae, indicatee=="bothdeepsnow")
indic.algae3$indicatee <- rep("Revelstoke pen")
indic.algae <- subset(indic.algae, indicatee !="bothdeepsnow" | is.na(indicatee=="TRUE"))

indic.algae <- rbind(indic.algae, indic.algae2, indic.algae3)
rm(indic.algae2, indic.algae3)

indic.algae$symbol <- rep("B")
indic.algae$symbol[indic.algae$p.value < 0.05] <- "A"
indic.algae$Group <- indic.algae$indicatee

df.algae2 <- merge(df.algae2, indic.algae[,c("Genus","A","B","symbol","Group")],
                   by=c("Group","Genus"), all=TRUE)
df.algae2 <- subset(df.algae2, is.na(Group)=="FALSE")
df.algae2$symbol[is.na(df.algae2$symbol)=="TRUE"] <- "B"

levels(df.algae2$Facet)[levels(df.algae2$Facet)=="shallow-snow"] <- "shallow-\nsnow"

df.algae2$symbol[df.algae2$Genus=="Trebouxia" & df.algae2$Group=="deep-snow"] <- "A"
df.algae2$symbol[df.algae2$Genus=="Chloroidium" & df.algae2$Group=="LARS"] <- "A"

df.algae2$Genus <- factor(df.algae2$Genus, 
                          levels=levels(df.algae2$Genus)[c(1:8,10:11,9)])

p2 <- ggplot(df.algae2, aes(x=Group, y=Genus, fill=sqrt(mean_color))) +
  geom_tile() +
  geom_text(aes(label=rounded), size=2.25) +
  geom_point(data=subset(df.algae2, symbol=="A"), 
             aes(x=Group, y=Genus)) +
  facet_grid(Phylum~Facet, scales="free", space="free") +
  scale_fill_gradient(low="white", high="firebrick3") +
  theme_bw() + plot_theme + 
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0), limits=rev) +
  theme(
    axis.text.x=element_blank(),
    #axis.text.x = element_text(angle=90, vjust=0.25, hjust=1, size=7),
    strip.text = element_text(size=7),
    axis.text.y = element_text(size=7),
    axis.title.x=element_blank(),
    legend.position="none") +
  labs(tag="b", y="Genus")

algae_legend <- cowplot::get_legend(
  ggplot(df.algae2, aes(x=Group, y=Genus, fill=sqrt(mean_color))) +
    geom_tile() +
    geom_text(aes(label=rounded), size=2.25) +
    geom_point(data=subset(df.algae2, symbol=="A"), 
               aes(x=Group, y=Genus)) +
    facet_grid(Phylum~Facet, scales="free", space="free") +
    scale_fill_gradient(low="white", high="firebrick3", 
                        breaks=c(0,2.236068, 3.162278,4.472136,5.477226,6.324555,7.061068),
                        labels=c("0","5", "10","20","30","40",">50"),
                        name="Relative\nabundance (%)") +
    guides(fill=guide_colorbar(title.position="top")) +
    theme_bw() + plot_theme + 
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0), limits=rev) +
    theme(
      axis.text.x=element_blank(),
      #axis.text.x = element_text(angle=90, vjust=0.25, hjust=1, size=7),
      strip.text = element_text(size=7),
      axis.text.y = element_text(size=7),
      axis.title.x=element_blank(),
      legend.title=element_text(size=7, hjust=0.5),
      legend.position="right", legend.text=element_text(size=7)) +
    labs(tag="b", y="Genus")
)

# (c) Plant ordination ####
temp <- bdiv.ordinations.all[[4]][["points"]][[3]]
temp$Herd <- factor(temp$Herd, levels=c("Hart Ranges", "North Cariboo", "Barkerville",
                                        "Wells Grey North", "Columbia North", "Central Selkirks",
                                        "Tonquin", "Brazeau", "Revelstoke pen", "LARS"))

temp$GroupHerd <- paste0(temp$Group, temp$Herd)
temp$GroupHerd <- factor(temp$GroupHerd,
                         levels=c("deep-snowHart Ranges", "deep-snowNorth Cariboo",
                                  "deep-snowBarkerville", "deep-snowWells Grey North",
                                  "deep-snowColumbia North", "deep-snowCentral Selkirks",
                                  "shallow-snowTonquin", "shallow-snowBrazeau", "Revelstoke penRevelstoke pen",
                                  "LARSLARS"))

p3 <- ggplot(data =  temp, aes(x=PC1, y=PC2, color=GroupHerd, shape=GroupHerd)) + 
  geom_point(size=1.5) +
  stat_ellipse(aes(group=Group), level=0.95, linewidth=0.5, color="darkgrey") +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=2)) +
  labs(x="PC1 (28.9%)", y="PC2 (22.5%)", tag="c") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.text=element_text(size=7),
        legend.position="none",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0),
        legend.box="horizontal",
        legend.key.height=unit(0.75,"line"),
        legend.key.width=unit(0.75,"line"),
        legend.spacing.y=unit(0.01,"cm")) +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,4,8)], 
                                             brewer.pal(4, "Paired")[c(3:4)], 
                                             "darkorchid1", "gray62"),
                     labels=c("Hart Ranges",  "Barkerville",
                              "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS")) +
  scale_shape_manual(name="Guide1", values=c(16,16,16,17,17,15,18),
                     labels=c("Hart Ranges",  "Barkerville",
                              "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS"))

# (d) Plant heat map ####
df.plants <- subset(diff.abund.GROUP[["PITS"]][["Genus"]], glm.eBH < 0.05)

df.plants2 <- subset(prev.data$PITS$Genus, Taxon %in% df.plants$Genus | Taxon %in% c("Hedlundia", "Sorbus"))
df.plants2$Genus <- df.plants2$Taxon
df.plants2$Taxon <- NULL
df.plants2 <- subset(df.plants2, Abundance > 1 | Genus %in% c("Hedlundia", "Sorbus"))

df.plants2 <- df.plants2[order(-df.plants2$Abundance), ]
df.plants2$Genus <- forcats::fct_inorder(df.plants2$Genus)

df.plants2[,c("Prevalence", "Abundance")] <- NULL
rownames(df.plants2) <- df.plants2$Genus
df.plants2$Genus <- NULL
df.plants2 <- as.data.frame(t(df.plants2))
df.plants2$Group <- rownames(df.plants2)
df.plants2 <- reshape2::melt(df.plants2)

df.plants2$Facet <- rep("population")
df.plants2$Facet[df.plants2$Group %in% c("Barkerville", "Central Selkirks", "Hart Ranges",
                                         "North Cariboo", "Columbia North", "Wells Grey North")] <- "deep-snow"
df.plants2$Facet[df.plants2$Group %in% c("Brazeau", "Tonquin")] <- "shallow-snow"
df.plants2$Group <- factor(df.plants2$Group, levels=c(
  "deep-snow", "shallow-snow", "Revelstoke pen", "LARS",
  "Hart Ranges", "North Cariboo", "Barkerville", "Wells Grey North", "Columbia North", "Central Selkirks",
  "Tonquin", "Brazeau"
))
df.plants2$Facet <- factor(df.plants2$Facet, levels=c(
  "population", "deep-snow", "shallow-snow"
))

df.plants2$rounded <- as.numeric(df.plants2$value)
df.plants2$rounded[df.plants2$rounded==0] <- "-"
df.plants2$rounded[df.plants2$rounded != "-"] <- NA

df.plants2$mean_color <- as.numeric(as.character(df.plants2$value))
df.plants2$mean_color[df.plants2$mean_color > 50] <- 50

df.plants2$Genus <- df.plants2$variable
df.plants2$variable <- NULL
df.plants2 <- merge(df.plants2,
                    df.plants[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")],
                    by="Genus", all.x=TRUE, all.y=FALSE)

indic.plants <- indicator.Gini.synthesis[[4]]
indic.plants <- subset(indic.plants, Genus %in% df.plants2$Genus)

indic.plants2 <- subset(indic.plants, indicatee=="bothdeepsnow")
indic.plants2$indicatee <- rep("deep-snow")
indic.plants3 <- subset(indic.plants, indicatee=="bothdeepsnow")
indic.plants3$indicatee <- rep("Revelstoke pen")
indic.plants <- subset(indic.plants, indicatee !="bothdeepsnow" | is.na(indicatee=="TRUE"))
indic.plants4 <- subset(indic.plants, indicatee=="wild")
indic.plants4$indicatee <- rep("deep-snow")
indic.plants5 <- subset(indic.plants, indicatee=="wild")
indic.plants5$indicatee <- rep("shallow-snow")
indic.plants <- subset(indic.plants, indicatee !="wild" | is.na(indicatee=="TRUE"))

indic.plants <- rbind(indic.plants, indic.plants2, indic.plants3, indic.plants4, indic.plants5)
rm(indic.plants2, indic.plants3, indic.plants4, indic.plants5)

indic.plants$symbol <- rep("B")
indic.plants$symbol[indic.plants$p.value < 0.05] <- "A"
indic.plants$Group <- indic.plants$indicatee
indic.plants$Group[indic.plants$Group=="Revelstokepen"] <- "Revelstoke pen"

df.plants2 <- merge(df.plants2, indic.plants[,c("Genus","A","B","symbol","Group")],
                    by=c("Group","Genus"), all=TRUE)
df.plants2 <- subset(df.plants2, is.na(Group)=="FALSE")
df.plants2$symbol[is.na(df.plants2$symbol)=="TRUE"] <- "B"

levels(df.plants2$Facet)[levels(df.plants2$Facet)=="shallow-snow"] <- "shallow-\nsnow"

df.plants2$Phylum <- rep("Plants")

# Add empty data for the missing herds.
mock.data <- df.plants2
mock.data <- subset(df.plants2, Facet=="deep-snow")
mock.data$Group[mock.data$Group=="Barkerville"] <- "North Cariboo"
mock.data$Group[mock.data$Group=="Hart Ranges"] <- "Wells Grey North"
mock.data$Group[mock.data$Group=="Central Selkirks"] <- "Columbia North"
mock.data$rounded = rep("-")
mock.data$value = rep(0)
mock.data$mean_color = rep(0)

df.plants2 <- rbind(df.plants2, mock.data)

p4 <- ggplot(df.plants2, aes(x=Group, y=Genus, fill=sqrt(mean_color))) +
  geom_tile() +
  geom_text(aes(label=rounded), size=2.25) +
  geom_point(data=subset(df.plants2, symbol=="A"), 
             aes(x=Group, y=Genus)) +
  facet_grid(Phylum~Facet, scales="free", space="free") +
  scale_fill_gradient(low="white", high="firebrick3") +
  theme_bw() + plot_theme + 
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0), limits=rev) +
  theme(
    axis.text.x = element_text(angle=90, vjust=0.25, hjust=1, size=7),
    strip.text.x = element_blank(),
    strip.text.y=element_text(size=7), 
    strip.background.x=element_blank(), 
    axis.text.y = element_text(size=7, vjust=0.5),
    legend.position="none") +
  labs(tag="d", x="Population / Herd", y="Genus")

# (-) Assemble ####
grid.arrange(p1, p2, p3, p4, ncol=2, widths=c(0.3, 0.7), heights=c(0.35,0.65))

p1 <- ggplotGrob(p1)
p2 <- ggplotGrob(p2)
p3 <- ggplotGrob(p3)
p4 <- ggplotGrob(p4)

p4$widths <- p2$widths
p1$widths <- p3$widths

grid.arrange(p1, p2, p3, p4, legend,
             layout_matrix=rbind(c(1,2),c(3,4),c(5,4)),
             widths=c(0.3, 0.7),
             heights=c(0.4,0.4,0.15))

plot <- arrangeGrob(p1, p2, p3, p4, legend,
                    layout_matrix=rbind(c(1,2),c(3,4),c(5,4)),
                    widths=c(0.3, 0.7),
                    heights=c(0.4,0.4,0.15))

ggsave("~/caribou/figures/raw_fig3_plant_algae_panels.jpg", plot, width=6.5, height=5.5, units="in", dpi=300)
ggsave("~/fig3_legend.jpg", algae_legend, width=6.5, height=4, units="in", dpi=300)

rm(p1, p2, p3, p4, plot, df.algae2, df.algae.hrd, df.algae,
   df.plants, df.plants2)

#### *FIG. 4: BACTERIA PANELS ####
# (a) Bacteria scatter plot ####
p1 <- ggplot() +
  geom_errorbar(data = alpha.means[[1]], aes(x=m.Observed_Extrap,
                                             ymin=m.Shannon_Extrap - sd.Shannon_Extrap,
                                             ymax=m.Shannon_Extrap + sd.Shannon_Extrap),
                color="grey", width=0) +
  geom_errorbarh(data = alpha.means[[1]], aes(y=m.Shannon_Extrap,
                                              xmin=m.Observed_Extrap - sd.Observed_Extrap,
                                              xmax=m.Observed_Extrap + sd.Observed_Extrap),
                 color="grey", height=0) +
  geom_point(data = sample_data.mb[[1]], aes(x=Observed_Extrap, y=Shannon_Extrap, color=Group), size=0.5) +
  
  geom_point(data = alpha.means[[1]], aes(x=m.Observed_Extrap, y=m.Shannon_Extrap,
                                          color=Group, shape=Group), size=3.5) +
  scale_shape_manual(values=c(16,17,15,18)) +
  theme_bw() + plot_theme +
  guides(shape=guide_legend(ncol=1, byrow=TRUE),
         color=guide_legend(ncol=1, byrow=TRUE, override.aes = list(size = 2, shape=c(16,17,15,18)))) +
  labs(x="ASV richness", y="Shannon diversity", tag="a") +
  scale_color_manual(values=c("chocolate1", "chartreuse3",
                              "darkorchid1", "gray62")) +
  theme(legend.position=c(0.76,0.15), legend.title=element_blank(),
        legend.key.height=unit(0.75,"line"),
        legend.key.width=unit(0.75,"line"),
        legend.text=element_text(size=7),
        legend.spacing.y=unit(0.002,'cm'))

# (b) Bacteria ordination ####
temp <- bdiv.ordinations.all[[1]][["points"]][[3]]
temp$Herd <- factor(temp$Herd, levels=c("Hart Ranges", "North Cariboo", "Barkerville",
                                        "Wells Grey North", "Columbia North", "Central Selkirks",
                                        "Tonquin", "Brazeau", "Revelstoke pen", "LARS"))

temp$GroupHerd <- paste0(temp$Group, temp$Herd)
temp$GroupHerd <- factor(temp$GroupHerd,
                         levels=c("deep-snowHart Ranges", "deep-snowNorth Cariboo",
                                  "deep-snowBarkerville", "deep-snowWells Grey North",
                                  "deep-snowColumbia North", "deep-snowCentral Selkirks",
                                  "shallow-snowTonquin", "shallow-snowBrazeau", "Revelstoke penRevelstoke pen",
                                  "LARSLARS"))

summary(bdiv.ordinations.all$X16S$models[[3]])

p2 <- ggplot(data =  temp, aes(PC1, PC2)) + 
  geom_point(aes(color = GroupHerd, shape = GroupHerd), size=1.25) +
  stat_ellipse(aes(group=Group), level=0.95, linewidth=0.5, color="darkgrey") +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=1)) +
  labs(x="PC1 (31.3%)", y="PC2 (11.6%)", tag="b") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.text=element_text(size=7),
        legend.position="right",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0),
        legend.box="horizontal",
        legend.key.height=unit(0.75,"line"),
        legend.key.width=unit(0.75,"line"),
        legend.spacing.y=unit(0.01,"cm")) +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)], 
                                             brewer.pal(4, "Paired")[c(3:4)], 
                                             "darkorchid1", "gray62"),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS")) +
  scale_shape_manual(name="Guide1", values=c(16,16,16,16,16,16,17,17,15,18),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS"))

# (c) Heat map ####
a.df.16S <- subset(diff.abund.GROUP[["X16S"]][["Genus"]], glm.eBH < 0.05)
a.df.16S.hrd <- subset(diff.abund.HERD[["X16S"]][["Genus"]], glm.eBH < 0.05)

df.16S <- subset(diff.abund.GROUP[["X16S"]][["Genus"]], glm.eBH < 0.05 | Genus %in% a.df.16S.hrd$Genus)
df.16S.hrd <- subset(diff.abund.HERD[["X16S"]][["Genus"]], glm.eBH < 0.05 | Genus %in% a.df.16S$Genus)

rm(a.df.16S, a.df.16S.hrd, df.16S.hrd)

df.16S2 <- subset(prev.data$X16S$Genus, Taxon %in% df.16S$Genus)
df.16S2$Genus <- df.16S2$Taxon
df.16S2$Taxon <- NULL
df.16S2 <- subset(df.16S2, Abundance > 0.5 | `deep-snow` > 1)

df.16S2[,c("Prevalence", "Abundance")] <- NULL
rownames(df.16S2) <- df.16S2$Genus
df.16S2$Genus <- NULL
df.16S2 <- as.data.frame(t(df.16S2))
df.16S2$Group <- rownames(df.16S2)
df.16S2 <- reshape2::melt(df.16S2)

df.16S2$Facet <- rep("population")
df.16S2$Facet[df.16S2$Group %in% c("Barkerville", "Central Selkirks", "Hart Ranges",
                                   "North Cariboo", "Columbia North", "Wells Grey North")] <- "deep-snow"
df.16S2$Facet[df.16S2$Group %in% c("Brazeau", "Tonquin")] <- "shallow-snow"
df.16S2$Group <- factor(df.16S2$Group, levels=c(
  "deep-snow", "shallow-snow", "Revelstoke pen", "LARS",
  "Hart Ranges", "North Cariboo", "Barkerville", "Wells Grey North", "Columbia North", "Central Selkirks",
  "Tonquin", "Brazeau"
))

df.16S2$Facet <- factor(df.16S2$Facet, levels=c(
  "population", "deep-snow", "shallow-snow"
))

df.16S2$rounded <- as.numeric(df.16S2$value)
df.16S2$rounded[df.16S2$rounded==0] <- "-"
df.16S2$rounded[df.16S2$rounded != "-"] <- NA

df.16S2$mean_color <- as.numeric(as.character(df.16S2$value))

df.16S2$Genus <- df.16S2$variable
df.16S2$variable <- NULL
df.16S2 <- merge(df.16S2,
                 df.16S[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")],
                 by="Genus", all.x=TRUE, all.y=FALSE)

indic.16S <- indicator.Gini.synthesis[[1]]
indic.16S <- subset(indic.16S, taxon %in% df.16S2$Genus)

indic.16S2 <- subset(indic.16S, indicatee=="deepsnow")
#indic.16S2$indicatee <- rep("deep-snow")
indic.16S3 <- subset(indic.16S, indicatee=="deepsnow")
#indic.16S3$indicatee <- rep("Revelstoke pen")
indic.16S <- subset(indic.16S, indicatee !="deepsnow" | is.na(indicatee=="TRUE"))

indic.16S <- rbind(indic.16S, indic.16S2, indic.16S3)
rm(indic.16S2, indic.16S3)

indic.16S$symbol <- rep("B")
indic.16S$symbol[indic.16S$p.value < 0.05] <- "A"
indic.16S$Group <- indic.16S$indicatee

df.16S2 <- merge(df.16S2, indic.16S[,c("Genus","A","B","symbol","Group")],
                 by=c("Group","Genus"), all=TRUE)
df.16S2 <- subset(df.16S2, is.na(Group)=="FALSE")
df.16S2$symbol[is.na(df.16S2$symbol)=="TRUE"] <- "B"

df.16S2 <- subset(df.16S2, Genus !="Treponema")
df.16S2 <- subset(df.16S2, Genus !="Uncl_Firmicutes")

levels(df.16S2$Facet)[levels(df.16S2$Facet)=="shallow-snow"] <- "shallow-\nsnow"

df.16S2$symbol[df.16S2$Genus=="Paramuribaculum" & df.16S2$Group=="deep-snow"] <- "A"
df.16S2$symbol[df.16S2$Genus=="Paramuribaculum" & df.16S2$Group=="Revelstoke pen"] <- "A"
df.16S2$symbol[df.16S2$Genus=="Prevotella" & df.16S2$Group=="LARS"] <- "A"

levels(df.16S2$Genus)[levels(df.16S2$Genus)=="Uncl_Ruminococcaceae"] <- "Ruminococcaceae spp."
levels(df.16S2$Genus)[levels(df.16S2$Genus)=="Uncl_Eubacteriaceae"] <- "Eubacteriaceae spp."
levels(df.16S2$Genus)[levels(df.16S2$Genus)=="Uncl_Clostridia"] <- "Uncl. Clostridia"
levels(df.16S2$Genus)[levels(df.16S2$Genus)=="Uncl_Lachnospiraceae"] <- "Lachnospiraceae spp."
levels(df.16S2$Genus)[levels(df.16S2$Genus)=="Uncl_Rikenellaceae"] <- "Rikenellaceae spp."
levels(df.16S2$Genus)[levels(df.16S2$Genus)=="Uncl_Muribaculaceae"] <- "Uncl. Muribaculaceae"
levels(df.16S2$Genus)[levels(df.16S2$Genus)=="Uncl_Bacteroidales"] <- "Uncl. Bacteroidales"

p3 <- ggplot(df.16S2, aes(x=Genus, y=Group, fill=sqrt(mean_color))) +
  geom_tile() +
  geom_text(aes(label=rounded), size=2.25) +
  geom_point(data=subset(df.16S2, symbol=="A"), 
             aes(x=Genus, y=Group)) +
  facet_grid(Facet~Phylum, scales="free", space="free") +
  scale_fill_gradient(low="white", high="firebrick3") +
  scale_color_manual(values=c("black", "darkgrey")) +
  theme_bw() + plot_theme + 
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0), limits=rev) +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1, size=7),
        strip.text = element_text(size=7),
        axis.text.y = element_text(size=7),
        legend.position="none") +
  labs(tag="c", x="Genus", y="Population / Herd\n")

legend <- cowplot::get_legend(
  ggplot(df.16S2, aes(x=Genus, y=Group, fill=sqrt(mean_color))) +
    geom_tile() +
    geom_text(aes(label=rounded), size=2.25) +
    geom_point(data=subset(df.16S2, symbol=="A"), 
               aes(x=Genus, y=Group)) +
    facet_grid(Facet~Phylum, scales="free", space="free") +
    scale_fill_gradient(low="white", high="firebrick3", 
                        breaks=c(0, 1, 2.236068, 3.162278, 4.1833, 5),
                        labels=c("0", "1", "5", "10", "17.5", "25"),
                        name="Relative\nabundance (%)") +
    scale_color_manual(values=c("black", "darkgrey")) +
    theme_bw() + plot_theme + 
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0), limits=rev) +
    theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1, size=7),
          strip.text = element_text(size=7),
          axis.text.y = element_text(size=7),
          legend.title=element_text(hjust=0.5, size=7),
          legend.position="right") +
    labs(tag="c", x="Genus", y="Population / Herd\n")
)


# (-) Assemble ####
grid.arrange(p1, p2, p3, layout_matrix=rbind(c(1,1,1,2,2,2,2),c(3,3,3,3,3,3,3)), heights=c(0.4,0.6))

p1 <- ggplotGrob(p1)
p2 <- ggplotGrob(p2)
p3 <- ggplotGrob(p3)

plot <- arrangeGrob(p1, p2, p3, layout_matrix=rbind(c(1,1,1,2,2,2,2),c(3,3,3,3,3,3,3)), heights=c(0.45,0.55))

ggsave("~/caribou/figures/raw_fig4_microbiome_panels.jpg", plot, width=6.5, height=7, units="in", dpi=300)
ggsave("~/caribou/figures/raw_fig4_legend.jpg", legend, width=6.5, height=4, units="in", dpi=300)

rm(p1, p2, p3, df.16S, df.16S2, plot, indic.16S)

#### *FIG S1: NEGATIVE CONTROLS AND SEQUENCING OVERVIEW ####
# (a) Read count distribution ####
read_counts.melt <- transform(merge(sample_data[,c("Group", "Herd")],
                                    read_counts.final, by=0, all=TRUE),
                              row.names=Row.names, Row.names=NULL)
read_counts.melt <- reshape2::melt(read_counts.melt)
read_counts.melt <- subset(read_counts.melt, variable %in% c("X16S", "X18S", "FITS", "PITS"))

levels(read_counts.melt$variable)[levels(read_counts.melt$variable)=="X16S"] <- "16S"
levels(read_counts.melt$variable)[levels(read_counts.melt$variable)=="X18S"] <- "18S"
levels(read_counts.melt$variable)[levels(read_counts.melt$variable)=="FITS"] <- "fungal ITS2"
levels(read_counts.melt$variable)[levels(read_counts.melt$variable)=="PITS"] <- "plant ITS2"

p1 <- ggplot(read_counts.melt, aes(x=value, fill=Group)) +
  geom_histogram() +
  facet_wrap(~variable, scales="free_x") +
  theme_bw() + plot_theme +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     breaks=c(0,2,4,6,8)) +
  geom_vline(xintercept = 2000, linetype="dashed") +
  labs(x="Number of reads", y="Number of samples", tag="a") +
  scale_fill_manual(values=c("chocolate1", "chartreuse3",
                             "darkorchid1", "gray62"))

# (b) Abundances of ASVs in the negative controls. ####
# Extract OTU tables with the negative control samples.
neg.ctrl <- list(
  X16S = subset(ctrl_data$extraction.blanks$X16S$otu_table,
                !Row.names %in% c("20-1328-additional", "20-1334", "MOCK")),
  X18S = subset(ctrl_data$extraction.blanks$X18S$otu_table,
                !Row.names %in% c("20-1328-additional", "20-1334-additional", "MOCK"))
)

neg.ctrl.plots <- list()

# Identify the mean and SD relative abundance for each OTU in the negative control samples. 
for(i in c(1:2)){
  neg.ctrl[[i]]$Row.names[neg.ctrl[[i]]$Row.names=="20-1328-additional"] <- "20-1328"
  neg.ctrl[[i]]$Row.names[neg.ctrl[[i]]$Row.names=="20-1334-additional"] <- "20-1334"
  
  temp <- neg.ctrl[[i]] %>%
    dplyr::group_by(Herd) %>%
    dplyr::select(., -Kit_16S, -Kit_18S, -Kit_FITS, -Kit_PITS, -Control, -Year, -Month, -Day, -Sex) %>%
    dplyr::summarise(dplyr::across(where(is.numeric), list(mean, sd), na.rm=TRUE))
  temp <- subset(temp, is.na(Herd)=="FALSE")
  
  if(i==1){
    temp1 <- temp[,c(1:3)]
    temp2 <- temp[,c(1,4:5)]
    
    temp1$ASVid <- strsplit(colnames(temp1)[2], "_")[[1]][1]
    temp2$ASVid <- strsplit(colnames(temp2)[2], "_")[[1]][1]
    
    colnames(temp1) <- c("Herd", "mean", "sd", "ASVid")
    colnames(temp2) <- c("Herd", "mean", "sd", "ASVid")
    
    temp <- rbind(temp1, temp2)
    rm(temp1, temp2)
  }
  
  if(i==2){
    temp$ASVid <- strsplit(colnames(temp)[2], "_")[[1]][1]
    colnames(temp) <- c("Herd", "mean", "sd", "ASVid")
  }
  
  tax <- ctrl_data$extraction.blanks[[i]]$tax_table
  tax$ASVid <- rownames(tax)
  tax <- subset(tax, ASVid %in% temp$ASVid)
  
  temp <- merge(tax, temp, by="ASVid", all=TRUE)
  
  if(i==1){
    temp$Genus[temp$Genus=="Escherichia/Shigella"] <- "Escherichia"
    temp$Genus[temp$Phylum=="Actinobacteria"] <- "Uncl. Micrococcales"}
  if(i==2){temp$Genus[temp$Family=="Magnoliophyta"] <- "Uncl. Magnoliophyta"}
  
  temp$Herd <- factor(temp$Herd, 
                      levels=c("Hart Ranges", "North Cariboo", "Barkerville", "Wells Grey North",
                               "Columbia North", "Central Selkirks", "Tonquin", "Brazeau",
                               "Revelstoke pen", "LARS"))
  
  neg.ctrl.plots[[i]] <- temp
  
  rm(temp, tax)
}

neg.ctrl.plots[[2]][,c("Domain", "Species")] <- NULL

neg.ctrl.plots <- rbind(neg.ctrl.plots[[1]], neg.ctrl.plots[[2]])
neg.ctrl.plots$mean[neg.ctrl.plots$Herd=="Columbia North" & neg.ctrl.plots$ASVid=="temp57"] <- 0.3882078
neg.ctrl.plots$ASV <- neg.ctrl.plots$Genus

neg.ctrl.plots$ASV <- factor(neg.ctrl.plots$ASV,
                             levels=c("Escherichia", "Uncl. Micrococcales",
                                      "Uncl. Magnoliophyta"))

p2 <- ggplot(neg.ctrl.plots,
             aes(x=mean, y=Herd)) +
  geom_bar(stat="identity") +
  facet_grid(~Genus, scales="free_x") +
  scale_y_discrete(limits=rev) +
  theme_bw() + plot_theme +
  theme(panel.spacing = unit(1.25, "lines")) +
  labs(tag="b", x="mean relative abundance (%)") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1)))

# (-) Assemble ####
grid.arrange(p1, p2, ncol=1, heights=c(0.6,0.4))

p1 <- ggplotGrob(p1)
p2 <- ggplotGrob(p2)

plot <- arrangeGrob(p1, p2, ncol=1, heights=c(0.6,0.4))

ggsave("~/FigS1_extraction_blanks.jpg", plot, width=6.5, height=6.5, units="in", dpi=300)

rm(p1, p2, neg.ctrl.plots, neg.ctrl, read_counts.melt)


#### *FIG S2: MOCK COMMUNITY ####
mock <- ctrl_data$mock.community$pseq.objects
mock.sub <- list()

for(p in c(1:length(mock))){
  mock.sub[[p]] <- mock[[p]]
  taxa_names(mock.sub[[p]]) <- refseq(mock.sub[[p]])
  
  for(i in 1:nrow(tax_table(mock.sub[[p]]))){
    for(j in 2:6){
      if(is.na(tax_table(mock.sub[[p]])[i,j])==TRUE){
        if(substr(tax_table(mock.sub[[p]])[i,j-1], 1, 4)=="Uncl"){
          tax_table(mock.sub[[p]])[i,j] <- tax_table(mock.sub[[p]])[i,j-1]}
        else {
          tax_table(mock.sub[[p]])[i,j] <- paste0("Uncl_", tax_table(mock.sub[[p]])[i,j-1])}}
    }}
}

names(mock.sub) <- names(mock)

# Extract OTU table and sequences
mock.comm.data <- list()

for(p in c(1:3)){
  mock.comm.data[[p]] <- as.data.frame(t(otu_table(transform_sample_counts(mock.sub[[p]], function(x) 100*x/sum(x)))))
  mock.comm.data[[p]] <- cbind(as.data.frame(cbind(tax_table(mock.sub[[p]]))),
                               mock.comm.data[[p]])
  mock.comm.data[[p]]$Kingdom <- NULL
  
  mock.comm.data[[p]] <- mock.comm.data[[p]][order(-mock.comm.data[[p]]$MOCK), ]
  mock.comm.data[[p]]$Taxon <- paste0("MockCom", seq(nrow(mock.comm.data[[p]])))
  levels(mock.comm.data[[p]]$Genus)[levels(mock.comm.data[[p]]$Genus)=="Escherichia/Shigella"] <- "Escherichia"
  mock.comm.data[[p]]$Taxon <- paste0(mock.comm.data[[p]]$Taxon, ": ", mock.comm.data[[p]]$Genus)
}

# Export sequence data for a MAFFT alignment and BLAST search.
fasta1 <- ShortRead(sread = DNAStringSet(rownames(mock.comm.data[[1]])),
                   id = BStringSet(mock.comm.data[[1]]$Taxon))
writeFasta(fasta1, file = "mock.16S.fasta")

fasta2 <- ShortRead(sread = DNAStringSet(rownames(mock.comm.data[[2]])),
                    id = BStringSet(mock.comm.data[[2]]$Taxon))
writeFasta(fasta2, file = "mock.18S.fasta")

fasta3 <- ShortRead(sread = DNAStringSet(rownames(mock.comm.data[[3]])),
                    id = BStringSet(mock.comm.data[[3]]$Taxon))
writeFasta(fasta3, file = "mock.FITS.fasta")

write.csv(mock.comm.data, "mock.community.csv")
rm(fasta)

# Run MAFFT alignment between mock community sequences and downloaded reference genomes

# (a) Phylogenetic tree for 16S sequences ####
mock.comm.tree.16S <- ape::read.tree(text='(((((
MC-16S-ASV1_Bacillus
:0.0000,
Bacillus_subtilis
:0.0000):0.0040,
MC-16S-ASV16_Bacillus
:0.0041):0.0053,
MC-16S-ASV13_Bacillus
:0.0028):0.0215,((
MC-16S-ASV3_Listeria
:0.0000,
Listeria_monocytogenes
:0.0000):0.0128,
MC-16S-ASV15_Uncl_Bacilli
:0.0063):0.0037):0.0035,(((
MC-16S-ASV6_Enterococcus
:0.0000,
Enterococcus_faecalis
:0.0000):0.0027,
MC-16S-ASV12_Enterococcus
:0.0000):0.0154,((
MC-16S-ASV4_Lactobacillus
:0.0000,
Lactobacillus_fermentum
:0.0000):0.0686,(
MC-16S-ASV7_Uncl_Clostridiales
:0.1101,((((
MC-16S-ASV5_Escherichia
:0.0000,
Escherichia_coli
:0.0000):0.0010,
MC-16S-ASV11_Escherichia
:0.0017):0.0169,(
MC-16S-ASV8_Salmonella
:0.0000,(
MC-16S-ASV10_Salmonella
:0.0000,
Salmonella_enterica
:0.0000):0.0032):0.0134):0.0620,((
MC-16S-ASV9_Pseudomonas
:0.0000,
Pseudomonas_aeruginosa
:0.0000):0.0248,
MC-16S-ASV14_Pseudomonas
:0.0138):0.0652):0.0804):0.0428):0.0151):0.0060,(
MC-16S-ASV2_Staphylococcus
:0.0000,
Staphylococcus_aureus
:0.0000):0.0312);')

library(ggtree)

# Visualize tree
ggtree(mock.comm.tree.16S) + 
  geom_tiplab() + 
  xlim(0,1) + 
  labs(tag="a") + 
  plot_theme + 
  theme(panel.border=element_blank()) #+ geom_text(aes(label=node))

mock.comm.tree.16S$tip.label
mock.comm.tree.16S$tip.label[mock.comm.tree.16S$tip.label=="\nMC-16S-ASV1_Bacillus\n"] <- "MC-16S-ASV1 Bacillus"
mock.comm.tree.16S$tip.label[mock.comm.tree.16S$tip.label=="\nBacillus_subtilis\n"] <- "Bacillus subtilis 16S"
mock.comm.tree.16S$tip.label[mock.comm.tree.16S$tip.label=="\nMC-16S-ASV16_Bacillus\n"] <- "MC-16S-ASV16 Bacillus"
mock.comm.tree.16S$tip.label[mock.comm.tree.16S$tip.label=="\nMC-16S-ASV13_Bacillus\n"] <- "MC-16S-ASV13 Bacillus"
mock.comm.tree.16S$tip.label[mock.comm.tree.16S$tip.label=="\nMC-16S-ASV3_Listeria\n"] <- "MC-16S-ASV3 Listeria"

mock.comm.tree.16S$tip.label[mock.comm.tree.16S$tip.label=="\nListeria_monocytogenes\n"] <- "Listeria monocytogenes"
mock.comm.tree.16S$tip.label[mock.comm.tree.16S$tip.label=="\nMC-16S-ASV15_Uncl_Bacilli\n"] <- "MC-16S-ASV15 Uncl. Bacilli"
mock.comm.tree.16S$tip.label[mock.comm.tree.16S$tip.label=="\nMC-16S-ASV6_Enterococcus\n"] <- "MC-16S-ASV6 Enterococcus"
mock.comm.tree.16S$tip.label[mock.comm.tree.16S$tip.label=="\nEnterococcus_faecalis\n"] <- "Enterococcus faecalis"
mock.comm.tree.16S$tip.label[mock.comm.tree.16S$tip.label=="\nMC-16S-ASV12_Enterococcus\n"] <- "MC-16S-ASV12 Enterococcus"

mock.comm.tree.16S$tip.label[mock.comm.tree.16S$tip.label=="\nMC-16S-ASV4_Lactobacillus\n"] <- "MC-16S-ASV4 Lactobacillus"
mock.comm.tree.16S$tip.label[mock.comm.tree.16S$tip.label=="\nLactobacillus_fermentum\n"] <- "Lactobacillus fermentum"
mock.comm.tree.16S$tip.label[mock.comm.tree.16S$tip.label=="\nMC-16S-ASV7_Uncl_Clostridiales\n"] <- "MC-16S-ASV7 Uncl. Clostridiales"
mock.comm.tree.16S$tip.label[mock.comm.tree.16S$tip.label=="\nMC-16S-ASV5_Escherichia\n"] <- "MC-16S-ASV5 Escherichia"
mock.comm.tree.16S$tip.label[mock.comm.tree.16S$tip.label=="\nEscherichia_coli\n"] <- "Escherichia coli"

mock.comm.tree.16S$tip.label[mock.comm.tree.16S$tip.label=="\nMC-16S-ASV11_Escherichia\n"] <- "MC-16S-ASV11 Escherichia"
mock.comm.tree.16S$tip.label[mock.comm.tree.16S$tip.label=="\nMC-16S-ASV8_Salmonella\n"] <- "MC-16S-ASV8 Salmonella"
mock.comm.tree.16S$tip.label[mock.comm.tree.16S$tip.label=="\nMC-16S-ASV10_Salmonella\n"] <- "MC-16S-ASV10 Salmonella"
mock.comm.tree.16S$tip.label[mock.comm.tree.16S$tip.label=="\nSalmonella_enterica\n"] <- "Salmonella enterica"
mock.comm.tree.16S$tip.label[mock.comm.tree.16S$tip.label=="\nMC-16S-ASV9_Pseudomonas\n"] <- "MC-16S-ASV9 Pseudomonas"

mock.comm.tree.16S$tip.label[mock.comm.tree.16S$tip.label=="\nPseudomonas_aeruginosa\n"] <- "Pseudomonas aeruginosa"
mock.comm.tree.16S$tip.label[mock.comm.tree.16S$tip.label=="\nMC-16S-ASV14_Pseudomonas\n"] <- "MC-16S-ASV14 Pseudomonas"
mock.comm.tree.16S$tip.label[mock.comm.tree.16S$tip.label=="\nMC-16S-ASV2_Staphylococcus\n"] <- "MC-16S-ASV2 Staphylococcus"
mock.comm.tree.16S$tip.label[mock.comm.tree.16S$tip.label=="\nStaphylococcus_aureus\n"] <- "Staphylococcus aureus"

p1 <- ggtree(mock.comm.tree.16S) + 
  geom_tiplab(size=3) + 
  xlim(0,0.5) + 
  labs(tag="a") + 
  plot_theme + 
  theme(panel.border=element_blank())

# (b) Relative abundances for 16S sequences ####
mock.bars.16S <- ctrl_data$mock.community$summary.data$protists
mock.bars.16S <- mock.bars.16S %>% 
  dplyr::group_by(Genus) %>% 
  dplyr::summarise_if(is.numeric, sum) %>%
  as.data.frame()
mock.bars.16S <- subset(mock.bars.16S, is.na(Genus)=="FALSE")

mock.bars.16S$MOCK <- as.numeric(as.character(mock.bars.16S$MOCK))
mock.bars.16S$Type <- rep("Actual")

temp <- ctrl_data$mock.community$actual.values
colnames(temp) <- c("Genus", "MOCK")
temp$Type <- rep("Expected")

mock.bars.16S <- rbind(mock.bars.16S, temp)

mock.bars.16S$Genus[mock.bars.16S$Genus=="Escherichia/Shigella"] <- "Escherichia"
mock.bars.16S$Genus[mock.bars.16S$Genus=="Limosilactobacillus"] <- "Lactobacillus"
mock.bars.16S$Genus[mock.bars.16S$Genus=="Metabacillus"] <- "Bacillus"

p2 <- ggplot(mock.bars.16S, aes(x=Type, y=MOCK, fill=Genus)) + geom_bar(stat="identity", color="black") +
  theme_bw() + plot_theme + 
  guides(fill=guide_legend(ncol=1, title.position="top", title.hjust=0.5)) +
  theme(axis.title.x=element_blank(),
        legend.position="right") +
  labs(x="Mock community", y="Relative abundance (%)", tag="b") +
  scale_y_continuous(expand=c(0,0))

grid.arrange(p1, p2, ncol=2)

# (c) Tree for 18S ####
mock.comm.tree.18S <- ape::read.tree(text='((
MC-18S-ASV1
:0.0000,(
MC-18S-ASV3
:0.0011,(
MC-18S-ASV6
:0.0000,((
MC-18S-ASV2
:0.0000,
Cryptococcusneoformans
:0.0000):0.0093,
MC-18S-ASV5
:0.0000):0.1626):0.0419):0.0016):0.0001,
Saccharomycescerevisiae
:0.0000,
MC-18S-ASV4
:0.0028);')

mock.comm.tree.18S$tip.label
mock.comm.tree.18S$tip.label[mock.comm.tree.18S$tip.label=="\nMC-18S-ASV1\n"] <- "MC-18S-ASV1 Saccharomyces"
mock.comm.tree.18S$tip.label[mock.comm.tree.18S$tip.label=="\nMC-18S-ASV3\n"] <- "MC-18S-ASV3 Saccharomyces"
mock.comm.tree.18S$tip.label[mock.comm.tree.18S$tip.label=="\nMC-18S-ASV6\n"] <- "MC-18S-ASV6 Saccharomyces"
mock.comm.tree.18S$tip.label[mock.comm.tree.18S$tip.label=="\nMC-18S-ASV2\n"] <- "MC-18S-ASV2 Cryptococcus"
mock.comm.tree.18S$tip.label[mock.comm.tree.18S$tip.label=="\nCryptococcusneoformans\n"] <- "Cryptococcus neoformans"
mock.comm.tree.18S$tip.label[mock.comm.tree.18S$tip.label=="\nMC-18S-ASV5\n"] <- "MC-18S-ASV5 Cryptococcus"
mock.comm.tree.18S$tip.label[mock.comm.tree.18S$tip.label=="\nSaccharomycescerevisiae\n"] <- "Saccharomyces cerevisiae"
mock.comm.tree.18S$tip.label[mock.comm.tree.18S$tip.label=="\nMC-18S-ASV4\n"] <- "MC-18S-ASV4 Saccharomyces"

p3 <- ggtree(mock.comm.tree.18S) + 
  geom_tiplab(size=3) + 
  xlim(0,1) + 
  labs(tag="c") + 
  plot_theme + 
  theme(panel.border=element_blank())

# (d) Relative abundances for 18S sequences ####
mock.bars.18S <- ctrl_data$mock.community$summary.data$FITS
mock.bars.18S <- mock.bars.18S %>% 
  dplyr::group_by(Genus) %>% 
  dplyr::summarise_if(is.numeric, sum) %>%
  as.data.frame()
mock.bars.18S <- subset(mock.bars.18S, is.na(Genus)=="FALSE")

mock.bars.18S$MOCK <- as.numeric(as.character(mock.bars.18S$MOCK))
mock.bars.18S$Type <- rep("Actual")

temp <- data.frame(
  Genus=c("Saccharomyces", "Cryptococcus"),
  MOCK=c(50, 50)
)
temp$Type <- rep("Expected")

mock.bars.18S <- rbind(mock.bars.18S, temp)
mock.bars.18S <- subset(mock.bars.18S, Genus !="Elaphomyces")

p4 <- ggplot(mock.bars.18S, aes(x=Type, y=MOCK, fill=Genus)) + geom_bar(stat="identity", color="black") +
  theme_bw() + plot_theme + 
  guides(fill=guide_legend(ncol=1, title.position="top", title.hjust=0.5)) +
  theme(axis.title.x=element_blank(),
        legend.position="right") +
  labs(x="Mock community", y="Relative abundance (%)", tag="d") +
  scale_y_continuous(expand=c(0,0))

grid.arrange(p3, p4, ncol=2)

# (-) Assemble ####
grid.arrange(p1, p2, p3, p4, ncol=2)

p1 <- ggplotGrob(p1)
p2 <- ggplotGrob(p2)
p3 <- ggplotGrob(p3)
p4 <- ggplotGrob(p4)

plot <- arrangeGrob(p1, p2, p3, p4, ncol=2)

ggsave("~/FigS2_mock_community.jpg", plot, width=6.5, height=8, units="in", dpi=300)

rm(p1, p2, p3, p4, mock.bars.16S, mock.bars.18S, mock.comm.tree.16S, mock.comm.tree.18S, plot)

#### *FIG S3: ALPHA DIVERSITY PANEL, FOUR AMPLICONS ####
# Prepare data ####
# 16S data
p1.data <- sample_data.mb[[1]] %>%
  dplyr::select(., -Kit_16S, -Kit_18S, -Kit_FITS, -Kit_PITS, -Year, -Month, -Day, -Sex)

p1.data$Group <- as.character(p1.data$Group)
p1.data$Herd <- as.character(p1.data$Herd)

p1.data$Facet <- p1.data$Group
p1.data$Facet[p1.data$Herd=="LARS"] <- "population"
p1.data$Facet[p1.data$Herd=="Revelstoke pen"] <- "population"

temp1 <- subset(p1.data, Group=="deep-snow")
temp1$Herd <- "deep-snow"
temp1$Facet <- "population"

temp2 <- subset(p1.data, Group=="shallow-snow")
temp2$Herd <- "shallow-snow"
temp2$Facet <- "population"

p1.data <- rbind(temp1, temp2, p1.data)
p1.data$Amplicon <- rep("16S")

# 18S data
p2.data <- sample_data.mb[[2]] %>%
  dplyr::select(., -Kit_16S, -Kit_18S, -Kit_FITS, -Kit_PITS, -Year, -Month, -Day, -Sex)

p2.data$Group <- as.character(p2.data$Group)
p2.data$Herd <- as.character(p2.data$Herd)

p2.data$Facet <- p2.data$Group
p2.data$Facet[p2.data$Herd=="LARS"] <- "population"
p2.data$Facet[p2.data$Herd=="Revelstoke pen"] <- "population"

temp1 <- subset(p2.data, Group=="deep-snow")
temp1$Herd <- "deep-snow"
temp1$Facet <- "population"

temp2 <- subset(p2.data, Group=="shallow-snow")
temp2$Herd <- "shallow-snow"
temp2$Facet <- "population"

p2.data <- rbind(temp1, temp2, p2.data)
p2.data$Amplicon <- rep("18S")

# FITS data
p3.data <- sample_data.mb[[3]] %>%
  dplyr::select(., -Kit_16S, -Kit_18S, -Kit_FITS, -Kit_PITS, -Year, -Month, -Day, -Sex)

p3.data$Group <- as.character(p3.data$Group)
p3.data$Herd <- as.character(p3.data$Herd)

p3.data$Facet <- p3.data$Group
p3.data$Facet[p3.data$Herd=="LARS"] <- "population"
p3.data$Facet[p3.data$Herd=="Revelstoke pen"] <- "population"

temp1 <- subset(p3.data, Group=="deep-snow")
temp1$Herd <- "deep-snow"
temp1$Facet <- "population"

temp2 <- subset(p3.data, Group=="shallow-snow")
temp2$Herd <- "shallow-snow"
temp2$Facet <- "population"

p3.data <- rbind(temp1, temp2, p3.data)
p3.data$Amplicon <- rep("fungal ITS2")

# PITS data
p4.data <- sample_data.mb[[4]] %>%
  dplyr::select(., -Kit_16S, -Kit_18S, -Kit_FITS, -Kit_PITS, -Year, -Month, -Day, -Sex)

p4.data$Group <- as.character(p4.data$Group)
p4.data$Herd <- as.character(p4.data$Herd)
p4.data$Facet <- p4.data$Group

p4.data$Facet[p4.data$Herd=="LARS"] <- "population"
p4.data$Facet[p4.data$Herd=="Revelstoke pen"] <- "population"

temp1 <- subset(p4.data, Group=="deep-snow")
temp1$Herd <- "deep-snow"
temp1$Facet <- "population"

temp2 <- subset(p4.data, Group=="shallow-snow")
temp2$Herd <- "shallow-snow"
temp2$Facet <- "population"

p4.data <- rbind(temp1, temp2, p4.data)
p4.data$Amplicon <- rep("plant ITS2")

# Assemble figure ####
p.data <- rbind(p1.data, p2.data, p3.data, p4.data)

p.data$Herd <- factor(p.data$Herd,
                       levels=c("deep-snow", "shallow-snow", "Revelstoke pen", "LARS",
                                "Hart Ranges", "North Cariboo", "Barkerville",
                                "Wells Grey North", "Columbia North", "Central Selkirks",
                                "Tonquin", "Brazeau"))
p.data$Group <- factor(p.data$Group,
                        levels=c("deep-snow", "shallow-snow", "Revelstoke pen", "LARS"))
p.data$Facet <- factor(p.data$Facet,
                        levels=c("population", "deep-snow", "shallow-snow"))

levels(p.data$Facet)[levels(p.data$Facet)=="shallow-snow"] <- "shallow-\nsnow"

p1 <- ggplot(p.data, aes(x=Herd, y=Observed_Extrap, fill=Group)) +
  geom_boxplot() +
  facet_grid(Amplicon~Facet, scales="free", space="free_x") +
  theme_bw() + plot_theme +
  scale_fill_manual(values=c("chocolate1", "chartreuse3",
                             "darkorchid1", "gray62")) +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        legend.position="none") +
  labs(x="Population / Herd", y="ASV richness", tag="a")

p2 <- ggplot(p.data, aes(x=Herd, y=Shannon_Extrap, fill=Group)) +
  geom_boxplot() +
  facet_grid(Amplicon~Facet, scales="free", space="free_x") +
  theme_bw() + plot_theme +
  scale_fill_manual(values=c("chocolate1", "chartreuse3",
                             "darkorchid1", "gray62")) +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        legend.position="none") +
  labs(x="Population / Herd", y="Shannon diversity", tag="b")

legend <- cowplot::get_legend(
  ggplot(p.data, aes(x=Herd, y=Shannon_Extrap, fill=Group)) +
    geom_boxplot() +
    facet_grid(Amplicon~Facet, scales="free", space="free_x") +
    theme_bw() + plot_theme +
    guides(fill=guide_legend(ncol=2, title.hjust=0.5)) +
    scale_fill_manual(values=c("chocolate1", "chartreuse3",
                               "darkorchid1", "gray62")) +
    theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1),
          legend.position="right",
          legend.title=element_blank()) +
    labs(x="Population / Herd", y="Shannon diversity", tag="b")
)

grid.arrange(p1, p2, legend,
             layout_matrix=rbind(c(1,2),c(3,3)),
             heights=c(0.9, 0.1))

p1 <- ggplotGrob(p1)
p2 <- ggplotGrob(p2)

plot <- arrangeGrob(p1, p2, legend,
                    layout_matrix=rbind(c(1,2),c(3,3)),
                    heights=c(0.9, 0.1))

ggsave("~/caribou/figures/raw_FigS3_adiv.overview.jpg", plot, width=9, height=6, units="in", dpi=300)

rm(p1, p2, p1.data, p2.data, p3.data, p4.data,
   temp1, temp2, p.data, plot)

#### *FIG S4: AITCHISON BETA DIVERSITY PANEL, FOUR AMPLICONS ####
# Prepare data ####
p.data <- list()
x.data <- list()

for(i in c(1:4)){
  temp <- bdiv.ordinations.all[[i]][["points"]][[3]]
  temp$Herd <- factor(temp$Herd, levels=c("Hart Ranges", "North Cariboo", "Barkerville",
                                                "Wells Grey North", "Columbia North", "Central Selkirks",
                                                "Tonquin", "Brazeau", "Revelstoke pen", "LARS"))
  
  temp$GroupHerd <- paste0(temp$Group, temp$Herd)
  temp$GroupHerd <- factor(temp$GroupHerd,
                              levels=c("deep-snowHart Ranges", "deep-snowNorth Cariboo",
                                       "deep-snowBarkerville", "deep-snowWells Grey North",
                                       "deep-snowColumbia North", "deep-snowCentral Selkirks",
                                       "shallow-snowTonquin", "shallow-snowBrazeau", "Revelstoke penRevelstoke pen",
                                       "LARSLARS"))
  
  if(i==1){temp$Amplicon <- "16S"}
  if(i==2){temp$Amplicon <- "18S"}
  if(i==3){temp$Amplicon <- "fungal ITS2"}
  if(i==4){temp$Amplicon <- "plant ITS2"}
  
  p.data[[i]] <- temp
  rm(temp)
  
  temp <- bdiv.ordinations.all[[i]][["models"]][[3]]
  x.data[[i]] <- c(100*temp[["CA"]][["eig"]][[1]]/temp[["CA"]]$tot.chi, 
                   100*temp[["CA"]][["eig"]][[2]]/temp[["CA"]]$tot.chi)
}

p1 <- ggplot(data =  p.data[[1]], aes(x=PC1, y=PC2, color=GroupHerd, shape=GroupHerd)) + 
  geom_point(size=2) +
  facet_wrap(~Amplicon) +
  stat_ellipse(aes(group=Group), level=0.95, linewidth=0.5, color="darkgrey") +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (31.3%)", y="PC2 (11.7%)") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)], 
                                             brewer.pal(4, "Paired")[c(3:4)], 
                                             "darkorchid1", "gray62"),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS")) +
  scale_shape_manual(name="Guide1", values=c(16,16,16,16,16,16,17,17,15,18),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS"))

p2 <- ggplot(data =  p.data[[2]], aes(x=PC1, y=PC2, color=GroupHerd, shape=GroupHerd)) + 
  geom_point(size=2) +
  facet_wrap(~Amplicon) +
  stat_ellipse(aes(group=Group), level=0.95, linewidth=0.5, color="darkgrey") +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (19.7%)", y="PC2 (15.3%)") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)], 
                                             brewer.pal(4, "Paired")[c(3:4)], 
                                             "darkorchid1", "gray62"),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS")) +
  scale_shape_manual(name="Guide1", values=c(16,16,16,16,16,16,17,17,15,18),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS"))

p3 <- ggplot(data =  p.data[[3]], aes(x=PC1, y=PC2, color=GroupHerd, shape=GroupHerd)) + 
  geom_point(size=2) +
  facet_wrap(~Amplicon) +
  stat_ellipse(aes(group=Group), level=0.95, linewidth=0.5, color="darkgrey") +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (22.5%)", y="PC2 (17.9%)") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)], 
                                             brewer.pal(4, "Paired")[c(3:4)], 
                                             "darkorchid1", "gray62"),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS")) +
  scale_shape_manual(name="Guide1", values=c(16,16,16,16,16,16,17,17,15,18),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS"))

p4 <- ggplot(data =  p.data[[4]], aes(x=PC1, y=PC2, color=GroupHerd, shape=GroupHerd)) + 
  geom_point(size=2) +
  facet_wrap(~Amplicon) +
  stat_ellipse(aes(group=Group), level=0.95, linewidth=0.5, color="darkgrey") +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (28.9%)", y="PC2 (22.5%)") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,4,8)], 
                                             brewer.pal(4, "Paired")[c(3:4)], 
                                             "darkorchid1", "gray62"),
                     labels=c("Hart Ranges", # "North Cariboo", 
                              "Barkerville",
                              # "Wells Grey North", "Columbia North", 
                              "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS")) +
  scale_shape_manual(name="Guide1", values=c(16,16,16,17,17,15,18),
                     labels=c("Hart Ranges", # "North Cariboo", 
                              "Barkerville",
                              # "Wells Grey North", "Columbia North", 
                              "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS"))

legend <- cowplot::get_legend(
  ggplot(data =  p.data[[2]], aes(x=PC1, y=PC2, color=GroupHerd, shape=GroupHerd)) + 
  geom_point(size=2) +
  facet_wrap(~Amplicon) +
  stat_ellipse(aes(group=Group), level=0.95, linewidth=0.5, color="darkgrey") +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="bottom",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)], 
                                             brewer.pal(4, "Paired")[c(3:4)], 
                                             "darkorchid1", "gray62"),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS")) +
  scale_shape_manual(name="Guide1", values=c(16,16,16,16,16,16,17,17,15,18),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS"))
)

# Assemble figure ####
grid.arrange(p1, p2, p3, p4, legend,
             layout_matrix=rbind(c(1,2),c(3,4),c(5,5)), heights=c(0.4,0.4,0.15))

plot <- arrangeGrob(p1, p2, p3, p4, legend,
                    layout_matrix=rbind(c(1,2),c(3,4),c(5,5)), heights=c(0.4,0.4,0.15))

ggsave("~/caribou/figures/raw_figS4_bdiv.overview.jpg", plot, width=6.5, height=7.5, units="in", dpi=300)

rm(p1, p2, p3, p4, legend, p.data, x.data, plot)

#### *FIG S5: READ DISTRIBUTIONS: 18S MAMMALS/PLANTS/ALGAE/PROTIST, FITS LICHENS ####
# (a) 18S distribution bar plots ####
pseq.18S.filter <- pseq_rchiv$raw$X18S
pseq.18S.filter.tax <- as.data.frame(cbind(tax_table(pseq_rchiv$raw$X18S)))

temp <- pseq.18S.filter.tax
for(i in c(1:ncol(temp))){temp[,i] <- as.character(temp[,i])}

# Clean taxonomy annotations
temp$Kingdom[temp$Class=="Mammalia"] <- "Mammalia"
temp$Kingdom[temp$Kingdom=="Phragmoplastophyta"] <- "Plant"
temp$Kingdom[temp$Genus=="Blastocystis"] <- "Protista"
temp$Kingdom[temp$Genus=="Dicksoniaceae"] <- "Plant"

tax_table(pseq.18S.filter) <- as.matrix(temp)

pseq.neg.ctrl <- subset_samples(pseq.18S.filter, sample_names(pseq.18S.filter) %in% c("K1-NC", "NC-K2"))
pseq.neg.ctrl <- prune_taxa(taxa_sums(pseq.neg.ctrl) > 0, pseq.neg.ctrl)

pseq.18S.filter <- subset_taxa(pseq.18S.filter, !(taxa_names(pseq.18S.filter) %in% taxa_names(pseq.neg.ctrl)))

pseq.18S.filter.glom <- tax_glom(pseq.18S.filter, taxrank="Kingdom")
pseq.18S.filter.glom <- subset_taxa(pseq.18S.filter.glom, Kingdom !="Eukaryota")
pseq.18S.filter.glom <- subset_taxa(pseq.18S.filter.glom, Kingdom !="Unassigned")

pseq.18S.filter.glom <- transform_sample_counts(pseq.18S.filter.glom, function(x) 100*x/sum(x))
taxa_names(pseq.18S.filter.glom) <- tax_table(pseq.18S.filter.glom)[,1]

temp <- as.data.frame(cbind(otu_table(pseq.18S.filter.glom)))
temp <- subset(temp, !rownames(temp) %in% c("20-1328-additional", "20-1334-additional"))
temp <- subset(temp, !rownames(temp) %in% c("K1-NC", "MOCK"))

rownames(sample_data.mb$X18S) <- sample_data.mb$X18S$SampleID

temp2 <- transform(merge(sample_data.mb$X18S[,c("Group", "Herd")], temp, by=0, all=TRUE),
                   row.names=Row.names, Row.names=NULL)

temp3 <- temp2 %>% dplyr::group_by(Group, Herd) %>%
  dplyr::summarise_if(is.numeric, mean, na.rm=TRUE)

temp3 <- reshape2::melt(temp3)
temp3$Herd <- factor(temp3$Herd,
                     levels=c("Hart Ranges", "North Cariboo", "Barkerville", "Wells Grey North", "Columbia North",
                              "Central Selkirks", "Tonquin", "Brazeau", "Revelstoke pen", "LARS"))

temp3$Group <- factor(temp3$Group, levels=c("deep-snow", "shallow-snow", "Revelstoke pen", "LARS"))
temp3$variable <- factor(temp3$variable, levels=c("Mammalia", "Animalia", "Fungi", "Plant", "Algae", "Protista"))

levels(temp3$Group)[levels(temp3$Group)=="Revelstoke pen"] <- "Rev.\npen"

p1 <- ggplot(temp3, aes(x=Herd, y=value, fill=variable)) +
  geom_bar(stat="identity", position="stack") +
  theme_bw() + plot_theme +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        legend.position="right") +
  guides(fill=guide_legend(title="Taxonomic group", ncol=1, title.hjust=0.5, title.position="top")) +
  facet_grid(~Group, scales="free_x", space="free") +
  scale_fill_manual(values=c("gray55", "darkorchid2", "gold2", 
                             "forestgreen", "darkseagreen1", "firebrick1")) +
  labs(x="Herd", y="Relative abundance (%)", tag="b")
  
# (b) Lichens as a proportion of fungi ####
temp.lich <- pseq.clr.sub$lichens
temp.fung <- pseq.clr.sub$FITS

temp.tax <- as.data.frame(cbind(tax_table(temp.fung)))
temp.tax$Kingdom[rownames(temp.tax) %in% taxa_names(temp.lich)] <- "Lichen"
tax_table(temp.fung) <- as.matrix(temp.tax)
temp.fung <- tax_glom(temp.fung, taxrank="Kingdom")
temp.fung <- transform_sample_counts(temp.fung, function(x) 100*x/sum(x))
taxa_names(temp.fung) <- tax_table(temp.fung)[,1]

rownames(sample_data.mb$FITS) <- sample_data.mb$FITS$SampleID
temp4 <- as.data.frame(cbind(otu_table(temp.fung)))
temp4 <- transform(merge(sample_data.mb$FITS[,c("Group", "Herd")], temp4, by=0, all=TRUE),
                   row.names=Row.names, Row.names=NULL)
temp4$Fungi <- NULL

temp4 <- temp4 %>% dplyr::group_by(Group, Herd) %>%
  dplyr::summarise_if(is.numeric, .funs=list(mean=mean, sd=sd), na.rm=TRUE)

temp4 <- reshape2::melt(temp4)

temp4a <- subset(temp4, variable=="mean")
temp4b <- subset(temp4, variable=="sd")
temp4a$variable <- NULL
temp4b$variable <- NULL
colnames(temp4a) <- c("Group", "Herd", "mean")
colnames(temp4b) <- c("Group", "Herd", "sd")
temp4 <- merge(temp4a, temp4b, by=c("Group", "Herd"), all=TRUE)

temp4$mean[temp4$Herd %in% c("Hart Ranges", "Columbia North")] <- 
  10 * temp4$mean[temp4$Herd %in% c("Hart Ranges", "Columbia North")]

temp4$sd.min <- 0
for(j in c(1:nrow(temp4))){
  if(temp4[j,"sd"] < temp4[j,"mean"]){
    temp4[j,"sd.min"] <- temp4[j,"mean"] - temp4[j,"sd"]
  }
}

temp4$Herd <- factor(temp4$Herd,
                     levels=c("Hart Ranges", "North Cariboo", "Barkerville", "Wells Grey North", "Columbia North",
                              "Central Selkirks", "Tonquin", "Brazeau", "Revelstoke pen", "LARS"))

temp4$Group <- as.character(temp4$Group)
temp4$Group[temp4$Herd=="Revelstoke pen"] <- "Rev.\npen"
temp4$Group <- factor(temp4$Group, levels=c("deep-snow", "shallow-snow", "Rev.\npen", "LARS"))

p2 <- ggplot(temp4, aes(x=Herd, y=mean)) +
  geom_bar(stat="identity", position="stack", fill="grey") +
  geom_errorbar(aes(x=Herd, ymin=sd.min, ymax=mean+sd),
                width=0.2) +
  theme_bw() + plot_theme +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
  guides(fill=guide_legend(title="Taxonomic group")) +
  facet_grid(~Group, scales="free_x", space="free") +
  labs(x="Herd", y="Lichens (as % of fungal ITS2 reads)", tag="a")

p1 <- ggplotGrob(p1)
p2 <- ggplotGrob(p2)
p2$widths <- p1$widths

grid.arrange(p2, p1, ncol=1)

plot <- arrangeGrob(p2, p1, ncol=1)
ggsave("~/caribou/figures/FigS5_lichen_protist_reads.jpg", plot, width=6, height=6.5, units="in", dpi=300)

rm(p1, p2, plot, temp4, temp4a, temp4b, pseq.18S.filter, temp, pseq.18S.filter.glom, pseq.18S.filter.tax,
   temp2, temp3, pseq.neg.ctrl)

#### *FIG. S6: LICHENS / MULTI-DISTANCE ORDINATION PANEL ####
# Prepare data ####
p.data <- list()
x.data <- list()

for(i in c(1,2,4,5)){
  temp <- bdiv.ordinations.all$lichens$points[[i]]
  temp$Herd <- factor(temp$Herd, levels=c("Hart Ranges", "North Cariboo", "Barkerville",
                                          "Wells Grey North", "Columbia North", "Central Selkirks",
                                          "Tonquin", "Brazeau", "Revelstoke pen", "LARS"))
  
  temp$GroupHerd <- paste0(temp$Group, temp$Herd)
  temp$GroupHerd <- factor(temp$GroupHerd,
                           levels=c("deep-snowHart Ranges", "deep-snowNorth Cariboo",
                                    "deep-snowBarkerville", "deep-snowWells Grey North",
                                    "deep-snowColumbia North", "deep-snowCentral Selkirks",
                                    "shallow-snowTonquin", "shallow-snowBrazeau", "Revelstoke penRevelstoke pen",
                                    "LARSLARS"))
  
  if(i==1){temp$distance <- "Bray-Curtis"}
  if(i==2){temp$distance <- "Jaccard"}
  if(i==4){temp$distance <- "weighted UniFrac"}
  if(i==5){temp$distance <- "unweighted UniFrac"}
  
  p.data[[i]] <- temp
  rm(temp)
  
  temp <- bdiv.ordinations.all$lichens$models[[i]]
  x.data[[i]] <- c(100*temp[["CA"]][["eig"]][[1]]/temp[["CA"]]$tot.chi, 
                   100*temp[["CA"]][["eig"]][[2]]/temp[["CA"]]$tot.chi)
}

p1 <- ggplot(data =  p.data[[1]], aes(x=MDS1, y=MDS2, color=GroupHerd, shape=GroupHerd)) + 
  geom_point(size=2) +
  facet_wrap(~distance) +
  stat_ellipse(aes(group=Group), level=0.95, linewidth=0.5, color="darkgrey") +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (24.7%)", y="PC2 (24.2%)") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)], 
                                             brewer.pal(4, "Paired")[c(3:4)], 
                                             "darkorchid1", "gray62"),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS")) +
  scale_shape_manual(name="Guide1", values=c(16,16,16,16,16,16,17,17,15,18),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS"))

p2 <- ggplot(data =  p.data[[2]], aes(x=MDS1, y=MDS2, color=GroupHerd, shape=GroupHerd)) + 
  geom_point(size=2) +
  facet_wrap(~distance) +
  stat_ellipse(aes(group=Group), level=0.95, linewidth=0.5, color="darkgrey") +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (18.1%)", y="PC2 (17.1%)") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)], 
                                             brewer.pal(4, "Paired")[c(3:4)], 
                                             "darkorchid1", "gray62"),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS")) +
  scale_shape_manual(name="Guide1", values=c(16,16,16,16,16,16,17,17,15,18),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS"))

p3 <- ggplot(data =  p.data[[4]], aes(x=MDS1, y=MDS2, color=GroupHerd, shape=GroupHerd)) + 
  geom_point(size=2) +
  facet_wrap(~distance) +
  stat_ellipse(aes(group=Group), level=0.95, linewidth=0.5, color="darkgrey") +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (40.0%)", y="PC2 (15.3%)") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)], 
                                             brewer.pal(4, "Paired")[c(3:4)], 
                                             "darkorchid1", "gray62"),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS")) +
  scale_shape_manual(name="Guide1", values=c(16,16,16,16,16,16,17,17,15,18),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS"))

p4 <- ggplot(data =  p.data[[5]], aes(x=MDS1, y=MDS2, color=GroupHerd, shape=GroupHerd)) + 
  geom_point(size=2) +
  facet_wrap(~distance) +
  stat_ellipse(aes(group=Group), level=0.95, linewidth=0.5, color="darkgrey") +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (25.5%)", y="PC2 (17.4%)") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)], 
                                             brewer.pal(4, "Paired")[c(3:4)], 
                                             "darkorchid1", "gray62"),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS")) +
  scale_shape_manual(name="Guide1", values=c(16,16,16,16,16,16,17,17,15,18),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS"))

legend <- cowplot::get_legend(
  ggplot(data =  p.data[[2]], aes(x=MDS1, y=MDS2, color=GroupHerd, shape=GroupHerd)) + 
    geom_point(size=2) +
    facet_wrap(~distance) +
    stat_ellipse(aes(group=Group), level=0.95, linewidth=0.5, color="darkgrey") +
    theme_bw() + plot_theme +
    guides(color=guide_legend(ncol=3)) +
    labs(x="PC1 (28.9%)", y="PC2 (22.5%)") +
    theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
          axis.line = element_line(colour = "black", size=0.5),
          legend.title=element_blank(),
          legend.position="bottom",
          legend.box="horizontal") +
    scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)], 
                                               brewer.pal(4, "Paired")[c(3:4)], 
                                               "darkorchid1", "gray62"),
                       labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                                "Wells Grey North", "Columbia North", "Central Selkirks",
                                "Tonquin", "Brazeau", "Revelstoke pen","LARS")) +
    scale_shape_manual(name="Guide1", values=c(16,16,16,16,16,16,17,17,15,18),
                       labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                                "Wells Grey North", "Columbia North", "Central Selkirks",
                                "Tonquin", "Brazeau", "Revelstoke pen","LARS"))
)

# Assemble figure ####
grid.arrange(p1, p2, p3, p4, legend,
             layout_matrix=rbind(c(1,2),c(3,4),c(5,5)), heights=c(0.4,0.4,0.15))

plot <- arrangeGrob(p1, p2, p3, p4, legend,
                    layout_matrix=rbind(c(1,2),c(3,4),c(5,5)), heights=c(0.4,0.4,0.15))

ggsave("~/caribou/figures/raw_figS6_lichen_bdiv_multidistance.jpg",
       plot, width=6.5, height=7.5, units="in", dpi=300)

rm(p1, p2, p3, p4, legend, p.data, x.data)

#### *FIG. S7: LICHEN AMONG-HERD VARIATION ####
p.data <- list()
x.data <- list()

for(i in c(1:5)){
  temp <- bdiv.ordinations.deepsnow$lichens$points[[i]]
  
  temp$Herd <- as.character(temp$Herd)
  temp$Herd <- factor(temp$Herd, levels=c("Hart Ranges", "North Cariboo", "Barkerville",
                                          "Wells Grey North", "Columbia North", "Central Selkirks"))
  if(i==1){temp$distance <- "Bray-Curtis"}
  if(i==2){temp$distance <- "Jaccard"}
  if(i==3){temp$distance <- "Aitchison"}
  if(i==4){temp$distance <- "weighted UniFrac"}
  if(i==5){temp$distance <- "unweighted UniFrac"}
  
  p.data[[i]] <- temp
  rm(temp)
  
  temp <- bdiv.ordinations.deepsnow$lichens$models[[i]]
  
  x.data[[i]] <- c(100*temp[["CA"]][["eig"]][[1]]/temp[["CA"]]$tot.chi, 
                   100*temp[["CA"]][["eig"]][[2]]/temp[["CA"]]$tot.chi)
}

p1 <- ggplot(p.data[[1]], aes(x=MDS1, y=MDS2, color=Herd)) + 
  geom_point(size=2, shape=16) +
  facet_wrap(~distance) +
  ggforce::geom_mark_ellipse(aes(color = Herd)) +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (20.9%)", y="PC2 (13.6%)") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)]),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks"))

p2 <- ggplot(p.data[[2]], aes(x=MDS1, y=MDS2, color=Herd)) + 
  geom_point(size=2, shape=16) +
  facet_wrap(~distance) +
  ggforce::geom_mark_ellipse(aes(color = Herd)) +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (14.5%)", y="PC2 (9.8%)") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)]),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks"))

p3 <- ggplot(p.data[[3]], aes(x=PC1, y=PC2, color=Herd)) + 
  geom_point(size=2, shape=16) +
  facet_wrap(~distance) +
  ggforce::geom_mark_ellipse(aes(color = Herd)) +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (20.0%)", y="PC2 (11.8%)") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)]),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks"))

p4 <- ggplot(p.data[[4]], aes(x=MDS1, y=MDS2, color=Herd)) + 
  geom_point(size=2, shape=16) +
  facet_wrap(~distance) +
  ggforce::geom_mark_ellipse(aes(color = Herd)) +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (31.1%)", y="PC2 (17.4%)") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)]),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks"))

p5 <- ggplot(p.data[[5]], aes(x=MDS1, y=MDS2, color=Herd)) + 
  geom_point(size=2, shape=16) +
  facet_wrap(~distance) +
  ggforce::geom_mark_ellipse(aes(color = Herd)) +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (32.8%)", y="PC2 (20.3%)") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)]),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks"))

p6 <- cowplot::get_legend(
  ggplot(p.data[[5]], aes(x=MDS1, y=MDS2, color=Herd)) + 
  geom_point(size=2, shape=16) +
  facet_wrap(~distance) +
  ggforce::geom_mark_ellipse(aes(color = Herd)) +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=1), line=FALSE) +
  labs(x="PC1 (32.8%)", y="PC2 (20.3%)") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="right") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)]),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks"))
)


grid.arrange(p1, p2, p3, p4, p5, p6, ncol=2)

plot <- arrangeGrob(p1, p2, p3, p4, p5, p6, ncol=2)

ggsave("~/caribou/figures/raw_figS7_lichen_deepsnow_herds.jpg",
       plot, width=6.5, height=7.5, units="in", dpi=300)

rm(p1, p2, p3, p4, p5, p6, plot, temp, p.data, x.data)

#### *FIG. S8: ALGAE ALPHA/BETA DIVERSITY AMONG POPULATIONS #### 
# (a) Scatter plot of alpha diversity ####
p5 <- ggplot() +
  geom_errorbar(data = alpha.means[[6]], aes(x=m.Observed_Extrap,
                                             ymin=m.Shannon_Extrap - sd.Shannon_Extrap,
                                             ymax=m.Shannon_Extrap + sd.Shannon_Extrap),
                color="grey", width=0) +
  geom_errorbarh(data = alpha.means[[6]], aes(y=m.Shannon_Extrap,
                                              xmin=m.Observed_Extrap - sd.Observed_Extrap,
                                              xmax=m.Observed_Extrap + sd.Observed_Extrap),
                 color="grey", height=0) +
  geom_point(data = sample_data.mb[[6]], aes(x=Observed_Extrap, y=Shannon_Extrap, color=Group), size=0.5) +
  geom_point(data = alpha.means[[6]], aes(x=m.Observed_Extrap, y=m.Shannon_Extrap,
                                          color=Group, shape=Group), size=3.5) +
  scale_shape_manual(values=c(16,17,15,18)) +
  theme_bw() + plot_theme +
  guides(shape=guide_legend(ncol=1, byrow=TRUE),
         color=guide_legend(ncol=1, byrow=TRUE, override.aes = list(size = 2, shape=c(16,17,15,18)))) +
  labs(x="ASV richness", y="Shannon diversity", tag="a") +
  scale_color_manual(values=c("chocolate1", "chartreuse3",
                              "darkorchid1", "gray62")) +
  theme(legend.position=c(0.78,0.23), 
        legend.title=element_blank(),
        legend.key.height=unit(0.75,"line"),
        legend.key.width=unit(0.75,"line"),
        #legend.margin=margin(0,0,0,0),
        #legend.box.margin=margin(0,0,0,0),
        legend.text=element_text(size=7),
        legend.spacing.y=unit(0.002,'cm'))

# (b) Series of ordinations excluding Aitchison ####
p.data <- list()
x.data <- list()

for(i in c(1,2,4,5)){
  temp <- bdiv.ordinations.all$algae$points[[i]]
  temp$Herd <- factor(temp$Herd, levels=c("Hart Ranges", "North Cariboo", "Barkerville",
                                          "Wells Grey North", "Columbia North", "Central Selkirks",
                                          "Tonquin", "Brazeau", "Revelstoke pen", "LARS"))
  
  temp$GroupHerd <- paste0(temp$Group, temp$Herd)
  temp$GroupHerd <- factor(temp$GroupHerd,
                           levels=c("deep-snowHart Ranges", "deep-snowNorth Cariboo",
                                    "deep-snowBarkerville", "deep-snowWells Grey North",
                                    "deep-snowColumbia North", "deep-snowCentral Selkirks",
                                    "shallow-snowTonquin", "shallow-snowBrazeau", "Revelstoke penRevelstoke pen",
                                    "LARSLARS"))
  
  if(i==1){temp$distance <- "Bray-Curtis"}
  if(i==2){temp$distance <- "Jaccard"}
  if(i==4){temp$distance <- "weighted UniFrac"}
  if(i==5){temp$distance <- "unweighted UniFrac"}
  
  p.data[[i]] <- temp
  
  temp <- bdiv.ordinations.deepsnow$algae$models[[i]]
  
  x.data[[i]] <- c(100*temp[["CA"]][["eig"]][[1]]/temp[["CA"]]$tot.chi, 
                   100*temp[["CA"]][["eig"]][[2]]/temp[["CA"]]$tot.chi)
  
  rm(temp)
}

p1 <- ggplot(data =  p.data[[1]], aes(x=MDS1, y=MDS2, color=GroupHerd, shape=GroupHerd)) + 
  geom_point(size=2) +
  facet_wrap(~distance) +
  stat_ellipse(aes(group=Group), level=0.95, linewidth=0.5, color="darkgrey") +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (35.4%)", y="PC2 (14.8%)", tag="b") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)], 
                                             brewer.pal(4, "Paired")[c(3:4)], 
                                             "darkorchid1", "gray62"),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS")) +
  scale_shape_manual(name="Guide1", values=c(16,16,16,16,16,16,17,17,15,18),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS"))

p2 <- ggplot(data =  p.data[[2]], aes(x=MDS1, y=MDS2, color=GroupHerd, shape=GroupHerd)) + 
  geom_point(size=2) +
  facet_wrap(~distance) +
  stat_ellipse(aes(group=Group), level=0.95, linewidth=0.5, color="darkgrey") +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (29.2%)", y="PC2 (11.6%)") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)], 
                                             brewer.pal(4, "Paired")[c(3:4)], 
                                             "darkorchid1", "gray62"),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS")) +
  scale_shape_manual(name="Guide1", values=c(16,16,16,16,16,16,17,17,15,18),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS"))

p3 <- ggplot(data =  p.data[[4]], aes(x=MDS1, y=MDS2, color=GroupHerd, shape=GroupHerd)) + 
  geom_point(size=2) +
  facet_wrap(~distance) +
  stat_ellipse(aes(group=Group), level=0.95, linewidth=0.5, color="darkgrey") +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (60.2%)", y="PC2 (19.8%)") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)], 
                                             brewer.pal(4, "Paired")[c(3:4)], 
                                             "darkorchid1", "gray62"),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS")) +
  scale_shape_manual(name="Guide1", values=c(16,16,16,16,16,16,17,17,15,18),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS"))

p4 <- ggplot(data =  p.data[[5]], aes(x=MDS1, y=MDS2, color=GroupHerd, shape=GroupHerd)) + 
  geom_point(size=2) +
  facet_wrap(~distance) +
  stat_ellipse(aes(group=Group), level=0.95, linewidth=0.5, color="darkgrey") +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (46.6%)", y="PC2 (15.6%)") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)], 
                                             brewer.pal(4, "Paired")[c(3:4)], 
                                             "darkorchid1", "gray62"),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS")) +
  scale_shape_manual(name="Guide1", values=c(16,16,16,16,16,16,17,17,15,18),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS"))

legend <- cowplot::get_legend(
  ggplot(data =  p.data[[2]], aes(x=MDS1, y=MDS2, color=GroupHerd, shape=GroupHerd)) + 
    geom_point(size=2) +
    facet_wrap(~distance) +
    stat_ellipse(aes(group=Group), level=0.95, linewidth=0.5, color="darkgrey") +
    theme_bw() + plot_theme +
    guides(color=guide_legend(ncol=2)) +
    labs(x="PC1 (28.9%)", y="PC2 (22.5%)") +
    theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
          axis.line = element_line(colour = "black", size=0.5),
          legend.title=element_blank(),
          legend.position="bottom",
          legend.box="horizontal") +
    scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)], 
                                               brewer.pal(4, "Paired")[c(3:4)], 
                                               "darkorchid1", "gray62"),
                       labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                                "Wells Grey North", "Columbia North", "Central Selkirks",
                                "Tonquin", "Brazeau", "Revelstoke pen","LARS")) +
    scale_shape_manual(name="Guide1", values=c(16,16,16,16,16,16,17,17,15,18),
                       labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                                "Wells Grey North", "Columbia North", "Central Selkirks",
                                "Tonquin", "Brazeau", "Revelstoke pen","LARS"))
)

# Assemble figure ####
grid.arrange(p5, legend, p1, p2, p3, p4, ncol=2)

p1 <- ggplotGrob(p1)
p2 <- ggplotGrob(p2)
p3 <- ggplotGrob(p3)
p4 <- ggplotGrob(p4)
p5 <- ggplotGrob(p5)

p2$widths <- p1$widths
p3$widths <- p1$widths
p4$widths <- p1$widths
p2$heights <- p1$heights
p3$heights <- p1$heights
p4$heights <- p1$heights

plot <- arrangeGrob(p5, legend, p1, p2, p3, p4, ncol=2)

ggsave("~/caribou/figures/raw_figS8_algae_population.jpg",
       plot, width=6.5, height=7.5, units="in", dpi=300)

rm(p1, p2, p3, p4, legend, p5, plot, p.data, x.data)

#### *FIG. S9: ALGAE ALPHA/BETA DIVERSITY AMONG HERDS ####
# (a) Richness box plot ####
p1.data <- subset(sample_data.mb$algae, Group=="deep-snow")
p1.data$Herd <- factor(p1.data$Herd)
p1.data$Herd <- factor(p1.data$Herd, levels=c("Hart Ranges", "North Cariboo", "Barkerville",
                                              "Wells Grey North", "Columbia North", "Central Selkirks"))

p1.data <- p1.data[,c("Herd", "Observed_Extrap", "Shannon_Extrap")]
p1.data <- reshape2::melt(p1.data)

levels(p1.data$variable)[levels(p1.data$variable)=="Observed_Extrap"] <- "ASV richness"
levels(p1.data$variable)[levels(p1.data$variable)=="Shannon_Extrap"] <- "Shannon diversity"

p1 <- ggplot(p1.data, aes(x=Herd, y=value, fill=Herd)) + geom_boxplot() +
  facet_wrap(~variable, scales="free_y") +
  theme_bw() + plot_theme +
  labs(tag="a") +
  scale_fill_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)])) +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        legend.position="none")

# (b) Series of ordinations ####
p.data <- list()
x.data <- list()

for(i in c(1:5)){
  temp <- bdiv.ordinations.deepsnow$algae$points[[i]]
  
  temp$Herd <- as.character(temp$Herd)
  temp$Herd <- factor(temp$Herd, levels=c("Hart Ranges", "North Cariboo", "Barkerville",
                                          "Wells Grey North", "Columbia North", "Central Selkirks"))
  if(i==1){temp$distance <- "Bray-Curtis"}
  if(i==2){temp$distance <- "Jaccard"}
  if(i==3){temp$distance <- "Aitchison"}
  if(i==4){temp$distance <- "weighted UniFrac"}
  if(i==5){temp$distance <- "unweighted UniFrac"}
  
  p.data[[i]] <- temp
  
  temp <- bdiv.ordinations.deepsnow$algae$models[[i]]
  x.data[[i]] <- c(100*temp[["CA"]][["eig"]][[1]]/temp[["CA"]]$tot.chi, 
                   100*temp[["CA"]][["eig"]][[2]]/temp[["CA"]]$tot.chi)
  
  rm(temp)
}

p2 <- ggplot(p.data[[1]], aes(x=MDS1, y=MDS2, color=Herd)) + 
  geom_point(size=2, shape=16) +
  facet_wrap(~distance) +
  ggforce::geom_mark_ellipse(aes(color = Herd)) +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (35.4%)", y="PC2 (14.8%)", tag="b") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)]),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks"))

p3 <- ggplot(p.data[[2]], aes(x=MDS1, y=MDS2, color=Herd)) + 
  geom_point(size=2, shape=16) +
  facet_wrap(~distance) +
  ggforce::geom_mark_ellipse(aes(color = Herd)) +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (29.2%)", y="PC2 (11.6%)") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)]),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks"))

p4 <- ggplot(p.data[[3]], aes(x=PC1, y=PC2, color=Herd)) + 
  geom_point(size=2, shape=16) +
  facet_wrap(~distance) +
  ggforce::geom_mark_ellipse(aes(color = Herd)) +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (25.1%)", y="PC2 (15.9%)") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)]),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks"))

p5 <- ggplot(p.data[[4]], aes(x=MDS1, y=MDS2, color=Herd)) + 
  geom_point(size=2, shape=16) +
  facet_wrap(~distance) +
  ggforce::geom_mark_ellipse(aes(color = Herd)) +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (60.2%)", y="PC2 (19.8%)") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)]),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks"))

p6 <- ggplot(p.data[[5]], aes(x=MDS1, y=MDS2, color=Herd)) + 
  geom_point(size=2, shape=16) +
  facet_wrap(~distance) +
  ggforce::geom_mark_ellipse(aes(color = Herd)) +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (46.6%)", y="PC2 (15.7%)") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)]),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks"))

p7 <- cowplot::get_legend(
  ggplot(p.data[[5]], aes(x=MDS1, y=MDS2, color=Herd)) + 
    geom_point(size=2, shape=16) +
    facet_wrap(~distance) +
    ggforce::geom_mark_ellipse(aes(color = Herd)) +
    theme_bw() + plot_theme +
    guides(color=guide_legend(ncol=1), line=FALSE) +
    labs(x="PC1 (32.8%)", y="PC2 (20.3%)") +
    theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
          axis.line = element_line(colour = "black", size=0.5),
          legend.title=element_blank(),
          legend.position="right") +
    scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)]),
                       labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                                "Wells Grey North", "Columbia North", "Central Selkirks"))
)

# (-) Assemble figure ####
grid.arrange(p1, p2, p3, p4, p5, p6, p7,
             layout_matrix=rbind(c(1,1,1),c(2,3,4),c(5,6,7)),
             heights=c(0.4,0.3,0.3))

p2 <- ggplotGrob(p2)
p3 <- ggplotGrob(p3)
p4 <- ggplotGrob(p4)
p5 <- ggplotGrob(p5)
p6 <- ggplotGrob(p6)

p2$widths <- p5$widths
p3$widths <- p5$widths
p4$widths <- p5$widths
p6$widths <- p5$widths

p3$heights <- p2$heights
p4$heights <- p2$heights
p5$heights <- p2$heights
p6$heights <- p2$heights

plot <- arrangeGrob(p1, p2, p3, p4, p5, p6, p7,
                    layout_matrix=rbind(c(1,1,1),c(2,3,4),c(5,6,7)),
                    heights=c(0.4,0.3,0.3))

ggsave("~/caribou/figures/raw_figS9_algae_deepsnow.jpg",
       plot, width=6.5, height=7.5, units="in", dpi=300)

rm(p1, p2, p3, p4, p5, p6, p7, plot, p.data, x.data)

#### *FIG. S10: PLANT ALPHA/BETA DIVERSITY AMONG POPULATIONS #### 
# (a) Scatter plot of alpha diversity ####
p5 <- ggplot() +
  geom_errorbar(data = alpha.means[[4]], aes(x=m.Observed_Extrap,
                                             ymin=m.Shannon_Extrap - sd.Shannon_Extrap,
                                             ymax=m.Shannon_Extrap + sd.Shannon_Extrap),
                color="grey", width=0) +
  geom_errorbarh(data = alpha.means[[4]], aes(y=m.Shannon_Extrap,
                                              xmin=m.Observed_Extrap - sd.Observed_Extrap,
                                              xmax=m.Observed_Extrap + sd.Observed_Extrap),
                 color="grey", height=0) +
  geom_point(data = sample_data.mb[[4]], aes(x=Observed_Extrap, y=Shannon_Extrap, color=Group), size=0.5) +
  geom_point(data = alpha.means[[4]], aes(x=m.Observed_Extrap, y=m.Shannon_Extrap,
                                          color=Group, shape=Group), size=3.5) +
  scale_shape_manual(values=c(16,17,15,18)) +
  theme_bw() + plot_theme +
  guides(shape=guide_legend(ncol=1, byrow=TRUE),
         color=guide_legend(ncol=1, byrow=TRUE, override.aes = list(size = 2, shape=c(16,17,15,18)))) +
  labs(x="ASV richness", y="Shannon diversity", tag="a") +
  scale_color_manual(values=c("chocolate1", "chartreuse3",
                              "darkorchid1", "gray62")) +
  theme(legend.position="none", 
        legend.title=element_blank(),
        legend.key.height=unit(0.75,"line"),
        legend.key.width=unit(0.75,"line"),
        legend.text=element_text(size=7),
        legend.spacing.y=unit(0.002,'cm'))

# (b) Series of ordinations excluding Aitchison ####
p.data <- list()
for(i in c(1,2,4,5)){
  temp <- bdiv.ordinations.all$PITS$points[[i]]
  temp$Herd <- factor(temp$Herd, levels=c("Hart Ranges", "Barkerville",
                                          "Central Selkirks",
                                          "Tonquin", "Brazeau", "Revelstoke pen", "LARS"))
  
  temp$GroupHerd <- paste0(temp$Group, temp$Herd)
  temp$GroupHerd <- factor(temp$GroupHerd,
                           levels=c("deep-snowHart Ranges", 
                                    "deep-snowBarkerville",
                                  "deep-snowCentral Selkirks",
                                    "shallow-snowTonquin", "shallow-snowBrazeau", "Revelstoke penRevelstoke pen",
                                    "LARSLARS"))
  
  if(i==1){temp$distance <- "Bray-Curtis"}
  if(i==2){temp$distance <- "Jaccard"}
  if(i==4){temp$distance <- "weighted UniFrac"}
  if(i==5){temp$distance <- "unweighted UniFrac"}
  
  p.data[[i]] <- temp
  rm(temp)
}

summary(bdiv.ordinations.all$PITS$models[[1]]) # 27.85%, 21.79%

p1 <- ggplot(data =  p.data[[1]], aes(x=MDS1, y=MDS2, color=GroupHerd, shape=GroupHerd)) + 
  geom_point(size=2) +
  facet_wrap(~distance) +
  stat_ellipse(aes(group=Group), level=0.95, linewidth=0.5, color="darkgrey") +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (27.9%)", y="PC2 (21.8%)", tag="b") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,4,8)], 
                                             brewer.pal(4, "Paired")[c(3:4)], 
                                             "darkorchid1", "gray62"),
                     labels=c("Hart Ranges", "Barkerville", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS")) +
  scale_shape_manual(name="Guide1", values=c(16,16,16,17,17,15,18),
                     labels=c("Hart Ranges", "Barkerville", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS"))

summary(bdiv.ordinations.all$PITS$models[[2]]) # 21.26%, 16.17%

p2 <- ggplot(data =  p.data[[2]], aes(x=MDS1, y=MDS2, color=GroupHerd, shape=GroupHerd)) + 
  geom_point(size=2) +
  facet_wrap(~distance) +
  stat_ellipse(aes(group=Group), level=0.95, linewidth=0.5, color="darkgrey") +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (21.3%)", y="PC2 (16.2%)") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,4,8)], 
                                             brewer.pal(4, "Paired")[c(3:4)], 
                                             "darkorchid1", "gray62"),
                     labels=c("Hart Ranges", "Barkerville", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS")) +
  scale_shape_manual(name="Guide1", values=c(16,16,16,17,17,15,18),
                     labels=c("Hart Ranges", "Barkerville", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS"))

summary(bdiv.ordinations.all$PITS$models[[4]]) # 37.93%, 22.77%

p3 <- ggplot(data =  p.data[[4]], aes(x=MDS1, y=MDS2, color=GroupHerd, shape=GroupHerd)) + 
  geom_point(size=2) +
  facet_wrap(~distance) +
  stat_ellipse(aes(group=Group), level=0.95, linewidth=0.5, color="darkgrey") +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (37.9%)", y="PC2 (22.9%)") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,4,8)], 
                                             brewer.pal(4, "Paired")[c(3:4)], 
                                             "darkorchid1", "gray62"),
                     labels=c("Hart Ranges", "Barkerville", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS")) +
  scale_shape_manual(name="Guide1", values=c(16,16,16,17,17,15,18),
                     labels=c("Hart Ranges", "Barkerville", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS"))

summary(bdiv.ordinations.all$PITS$models[[5]]) # 23.89%, 22.85%

p4 <- ggplot(data =  p.data[[5]], aes(x=MDS1, y=MDS2, color=GroupHerd, shape=GroupHerd)) + 
  geom_point(size=2) +
  facet_wrap(~distance) +
  stat_ellipse(aes(group=Group), level=0.95, linewidth=0.5, color="darkgrey") +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (23.9%)", y="PC2 (22.9%)") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,4,8)], 
                                             brewer.pal(4, "Paired")[c(3:4)], 
                                             "darkorchid1", "gray62"),
                     labels=c("Hart Ranges", "Barkerville", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS")) +
  scale_shape_manual(name="Guide1", values=c(16,16,16,17,17,15,18),
                     labels=c("Hart Ranges", "Barkerville", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS"))

legend <- cowplot::get_legend(
  ggplot(data =  p.data[[2]], aes(x=MDS1, y=MDS2, color=GroupHerd, shape=GroupHerd)) + 
    geom_point(size=2) +
    facet_wrap(~distance) +
    stat_ellipse(aes(group=Group), level=0.95, linewidth=0.5, color="darkgrey") +
    theme_bw() + plot_theme +
    guides(color=guide_legend(ncol=2)) +
    labs(x="PC1 (28.9%)", y="PC2 (22.5%)") +
    theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
          axis.line = element_line(colour = "black", size=0.5),
          legend.title=element_blank(),
          legend.position="bottom",
          legend.box="horizontal") +
    scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,4,8)], 
                                               brewer.pal(4, "Paired")[c(3:4)], 
                                               "darkorchid1", "gray62"),
                       labels=c("Hart Ranges", "Barkerville", "Central Selkirks",
                                "Tonquin", "Brazeau", "Revelstoke pen","LARS")) +
    scale_shape_manual(name="Guide1", values=c(16,16,16,17,17,15,18),
                       labels=c("Hart Ranges", "Barkerville", "Central Selkirks",
                                "Tonquin", "Brazeau", "Revelstoke pen","LARS"))
)

# Assemble figure ####
grid.arrange(p5, legend, p1, p2, p3, p4, ncol=2)

p1 <- ggplotGrob(p1)
p2 <- ggplotGrob(p2)
p3 <- ggplotGrob(p3)
p4 <- ggplotGrob(p4)
p5 <- ggplotGrob(p5)

p2$widths <- p1$widths
p3$widths <- p1$widths
p4$widths <- p1$widths
p2$heights <- p1$heights
p3$heights <- p1$heights
p4$heights <- p1$heights

plot <- arrangeGrob(p5, legend, p1, p2, p3, p4, ncol=2)

ggsave("~/FigS10_plant_population.jpg",
       plot, width=6.5, height=7.5, units="in", dpi=300)

rm(p1, p2, p3, p4, p5, p.data, plot) 

#### *FIG. S11: BACTERIAL RELATIVE ABUNDANCES + METHANOGENS ####
# (a) Relative abundance bar chart ####
y1 <- tax_glom(pseq.clr.sub$X16S, taxrank="Family", NArm = FALSE)
y1 <- subset_samples(y1, SampleID !="20-1348")

y5 <- as.data.frame(cbind(tax_table(y1)))
y5$Phylum <- as.character(y5$Phylum)
y5$Phylum[y5$Phylum %in% c("Actinobacteria",
                           "Spirochaetes", "Lentisphaerae",
                           "Fibrobacteres",
                           "Planctomycetes",
                           "Cyanobacteria/Chloroplast",
                           "Verrucomicrobia",
                           "Chloroflexi",
                           "Elusimicrobia",
                           "Synergistetes")] <- "Other"
y5$Family <- as.character(y5$Family)
y5$Family[y5$Phylum=="Other"] <- "Taxa <1% Abundance"

y5$Family[y5$Phylum=="Bacteroidetes" & !y5$Family %in% c("Bacteroidaceae", "Muribaculaceae")] <- 
  "Other Bacteroidetes"

y5$Family[y5$Phylum=="Firmicutes" & !y5$Family %in% c("Ruminococcaceae", "Lachnospiraceae")] <- "Other Firmicutes"

y5$Family[y5$Phylum=="Proteobacteria"] <- "Proteobacteria"

y5$Phylum <- as.factor(y5$Phylum)
y5$Class <- as.factor(y5$Class)
seq <- rownames(y5)
y5 <- tax_table(y5)
rownames(y5) <- seq
colnames(y5) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
tax_table(y1) <- y5

y2 <- tax_glom(y1, taxrank="Family")
y2 <- merge_samples(y2, group="Herd")
y2 <- prune_taxa(taxa_sums(y2) > 0, y2)
y2 <- transform_sample_counts(y2, function(x) 100*x/sum(x))

yt <- tax_glom(y1, taxrank="Family")
yt <- merge_samples(yt, group="Group")
yt <- prune_taxa(taxa_sums(yt) > 0, yt)
yt <- transform_sample_counts(yt, function(x) 100*x/sum(x))

y4 <- psmelt(y2)
yt <- psmelt(yt)

yt <- yt %>% dplyr::group_by(Sample, Phylum, Family) %>% dplyr::summarise(Abundance=sum(Abundance))
y4 <- y4 %>% dplyr::group_by(Group, Sample, Phylum, Family) %>% dplyr::summarise(Abundance=sum(Abundance))

yt$Facet <- rep("population")
y4 <- subset(y4, Group %in% c(1, 2))

y4$Facet <- rep("deep-snow herds")
y4$Facet[y4$Group==2] <- "shallow-snow\nherds"
y4$Group <- NULL

y4 <- rbind(yt, y4)

y4$Phylum <- factor(as.character(y4$Phylum), levels=c("Other","Proteobacteria","Bacteroidetes", "Firmicutes"))
y4$Family <- factor(y4$Family, levels=c("Taxa <1% Abundance",
                                        "Proteobacteria",
                                        "Muribaculaceae", "Bacteroidaceae", "Other Bacteroidetes", 
                                        "Lachnospiraceae", "Ruminococcaceae", "Other Firmicutes"))
y4 <- y4[order(y4$Phylum, y4$Family),]
colorCount <- length(unique(y4$Family))

myColors <- c("#999999", #Other
              "#99CC99", #Other Proteobacteri
              "#FFCC00", ##Other Bacteroidetes
              "#FF9900", ##Muribaculaceae
              "#CC9933", ##Bacteroidaceae
              "#3399FF", ##Other Firmicutes
              "#99CCFF", #Lachnospiraceae
              "#6699CC") #Ruminococcaceae
names(myColors) <- levels(y4$Family)

y4$Sample <- factor(y4$Sample, levels=c("deep-snow", "shallow-snow", "Revelstoke pen", "LARS",
                                   "Hart Ranges", "North Cariboo", "Barkerville", "Wells Grey North",
                                    "Columbia North", "Central Selkirks", "Tonquin", "Brazeau"))
y4$Facet <- factor(y4$Facet, levels=c("population", "deep-snow herds", "shallow-snow\nherds"))

p1 <- ggplot(y4, aes(x=Sample, y=Abundance, fill=Family)) +
  geom_bar(stat="identity") +
  facet_grid(~Facet, scales="free_x", space="free") +
  scale_fill_manual(name="Class", values=myColors) +
  scale_y_continuous(expand=c(0,0), limits=c(0, 102)) +
  #scale_x_discrete(labels=function(Type) str_wrap(Location.Tag, width=3)) +
  guides(fill=guide_legend(ncol=1)) +
  theme_bw() + plot_theme +
  theme(#panel.spacing=unit(2, "lines"),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(color="black", size = 8, hjust=1, vjust=1, angle=45),
        axis.text.y=element_text(color="black", size = 8, vjust=1),
        #strip.background=element_blank(), 
        legend.title=element_blank(),
        legend.position="right",
        legend.key.size=unit(1, "line")) +
  labs(x="Population / Herd", y = "Mean relative abundance (%)", tag="a")

# (b) Methanogens ####
yt <- tax_glom(ctrl_data$archaea, taxrank="Genus")
yt <- psmelt(yt)
yt <- subset(yt, Abundance > 0)
x <- yt %>% dplyr::group_by(Genus, Herd) %>% dplyr::count()
yt <- yt %>%
  dplyr::group_by(Genus, Herd) %>%
  dplyr::summarise(meanAbund = mean(Abundance),
                   sdAbund = sd(Abundance))
yt <- merge(x, yt, by=c("Genus", "Herd"), all=TRUE)

x <- sample_data.mb$X16S %>% dplyr::group_by(Herd) %>% dplyr::count()
colnames(x) <- c("Herd", "Total")

yt <- merge(yt, x, by="Herd", all=TRUE)

yt$prev <- 100 * yt$n / yt$Total

yt$Herd <- factor(yt$Herd, levels=c("Hart Ranges", "North Cariboo", "Barkerville", "Wells Grey North",
                                        "Columbia North", "Central Selkirks", "Tonquin", "Brazeau",
                                        "Revelstoke pen", "LARS"))

yt$lowlimit <- yt$meanAbund - yt$sdAbund
yt$lowlimit[yt$lowlimit < 0] <- 0

yt$meanAbund[yt$meanAbund==0] <- NA
yt <- subset(yt, is.na(Genus)=="FALSE")

p2 <- ggplot(yt, aes(x=Herd, y=prev)) +
  geom_point(aes(size=meanAbund)) +
  facet_wrap(~Genus) +
  scale_x_discrete(limits=levels(yt$Herd)) +
  scale_y_continuous(limits=c(10,105)) +
  theme_bw() + plot_theme +
  guides(size = guide_legend(title="Mean relative\nabundance (%)",
                             title.hjust=0.5)) +
  geom_vline(xintercept=c(6.5, 8.5, 9.5), color="grey50") +
  theme(
    axis.text.x=element_text(color="black", size = 8, hjust=1, vjust=1, angle=45)
  ) +
  labs(x = "Herd", y="Prevalence in herd (%)", tag="b")

# (-) Assemble figure ####
grid.arrange(p1, p2, nrow=2, heights=c(0.6, 0.4))

plot <- arrangeGrob(p1, p2, nrow=2, heights=c(0.6, 0.4))

ggsave("~/raw_figS11_16S_bars_methanogens.jpg", plot, 
       width=7, height=8, units="in")

rm(y1, y5, seq, yt, y4, x, plot, p1, p2, y2)

#### *FIG. S12: BACTERIA MULTI-DISTANCE ORDINATION PANEL ####
# Prepare data ####
p.data <- list()
x.data <- list()
for(i in c(1,2,4,5)){
  temp <- bdiv.ordinations.all$protists$points[[i]]
  temp$Herd <- factor(temp$Herd, levels=c("Hart Ranges", "North Cariboo", "Barkerville",
                                          "Wells Grey North", "Columbia North", "Central Selkirks",
                                          "Tonquin", "Brazeau", "Revelstoke pen", "LARS"))
  
  temp$GroupHerd <- paste0(temp$Group, temp$Herd)
  temp$GroupHerd <- factor(temp$GroupHerd,
                           levels=c("deep-snowHart Ranges", "deep-snowNorth Cariboo",
                                    "deep-snowBarkerville", "deep-snowWells Grey North",
                                    "deep-snowColumbia North", "deep-snowCentral Selkirks",
                                    "shallow-snowTonquin", "shallow-snowBrazeau", "Revelstoke penRevelstoke pen",
                                    "LARSLARS"))
  
  if(i==1){temp$distance <- "Bray-Curtis"}
  if(i==2){temp$distance <- "Jaccard"}
  if(i==4){temp$distance <- "weighted UniFrac"}
  if(i==5){temp$distance <- "unweighted UniFrac"}
  
  p.data[[i]] <- temp
  
  temp <- bdiv.ordinations.all$protists$models[[i]]
  x.data[[i]] <- c(100*temp[["CA"]][["eig"]][[1]]/temp[["CA"]]$tot.chi, 
                   100*temp[["CA"]][["eig"]][[2]]/temp[["CA"]]$tot.chi)
  rm(temp)
}

p1 <- ggplot(data =  p.data[[1]], aes(x=MDS1, y=MDS2, color=GroupHerd, shape=GroupHerd)) + 
  geom_point(size=2) +
  facet_wrap(~distance) +
  stat_ellipse(aes(group=Group), level=0.95, linewidth=0.5, color="darkgrey") +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (28.5%)", y="PC2 (13.1%)") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)], 
                                             brewer.pal(4, "Paired")[c(3:4)], 
                                             "darkorchid1", "gray62"),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS")) +
  scale_shape_manual(name="Guide1", values=c(16,16,16,16,16,16,17,17,15,18),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS"))

p2 <- ggplot(data =  p.data[[2]], aes(x=MDS1, y=MDS2, color=GroupHerd, shape=GroupHerd)) + 
  geom_point(size=2) +
  facet_wrap(~distance) +
  stat_ellipse(aes(group=Group), level=0.95, linewidth=0.5, color="darkgrey") +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (19.4%)", y="PC2 (9.1%)") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)], 
                                             brewer.pal(4, "Paired")[c(3:4)], 
                                             "darkorchid1", "gray62"),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS")) +
  scale_shape_manual(name="Guide1", values=c(16,16,16,16,16,16,17,17,15,18),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS"))

p3 <- ggplot(data =  p.data[[4]], aes(x=MDS1, y=MDS2, color=GroupHerd, shape=GroupHerd)) + 
  geom_point(size=2) +
  facet_wrap(~distance) +
  stat_ellipse(aes(group=Group), level=0.95, linewidth=0.5, color="darkgrey") +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (46.4%)", y="PC2 (13.0%)") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)], 
                                             brewer.pal(4, "Paired")[c(3:4)], 
                                             "darkorchid1", "gray62"),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS")) +
  scale_shape_manual(name="Guide1", values=c(16,16,16,16,16,16,17,17,15,18),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS"))

p4 <- ggplot(data =  p.data[[5]], aes(x=MDS1, y=MDS2, color=GroupHerd, shape=GroupHerd)) + 
  geom_point(size=2) +
  facet_wrap(~distance) +
  stat_ellipse(aes(group=Group), level=0.95, linewidth=0.5, color="darkgrey") +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (27.2%)", y="PC2 (12.5%)") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)], 
                                             brewer.pal(4, "Paired")[c(3:4)], 
                                             "darkorchid1", "gray62"),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS")) +
  scale_shape_manual(name="Guide1", values=c(16,16,16,16,16,16,17,17,15,18),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS"))

legend <- cowplot::get_legend(
  ggplot(data =  p.data[[2]], aes(x=MDS1, y=MDS2, color=GroupHerd, shape=GroupHerd)) + 
    geom_point(size=2) +
    facet_wrap(~distance) +
    stat_ellipse(aes(group=Group), level=0.95, linewidth=0.5, color="darkgrey") +
    theme_bw() + plot_theme +
    guides(color=guide_legend(ncol=3)) +
    labs(x="PC1 (28.9%)", y="PC2 (22.5%)") +
    theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
          axis.line = element_line(colour = "black", size=0.5),
          legend.title=element_blank(),
          legend.position="bottom",
          legend.box="horizontal") +
    scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)], 
                                               brewer.pal(4, "Paired")[c(3:4)], 
                                               "darkorchid1", "gray62"),
                       labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                                "Wells Grey North", "Columbia North", "Central Selkirks",
                                "Tonquin", "Brazeau", "Revelstoke pen","LARS")) +
    scale_shape_manual(name="Guide1", values=c(16,16,16,16,16,16,17,17,15,18),
                       labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                                "Wells Grey North", "Columbia North", "Central Selkirks",
                                "Tonquin", "Brazeau", "Revelstoke pen","LARS"))
)

# Assemble figure ####
grid.arrange(p1, p2, p3, p4, legend,
             layout_matrix=rbind(c(1,2),c(3,4),c(5,5)), heights=c(0.4,0.4,0.15))

plot <- arrangeGrob(p1, p2, p3, p4, legend,
                    layout_matrix=rbind(c(1,2),c(3,4),c(5,5)), heights=c(0.4,0.4,0.15))

ggsave("~/caribou/figures/raw_figS12_bacteria_multidistance.jpg",
       plot, width=6.5, height=7.5, units="in", dpi=300)

rm(p1, p2, p3, p4, legend, p.data, x.data)

#### *FIG. S13: BACTERIA AMONG-HERD VARIATION ####
# (a) Richness box plot ####
p1.data <- subset(sample_data.mb$protists, Group=="deep-snow")
p1.data$Herd <- factor(p1.data$Herd)
p1.data$Herd <- factor(p1.data$Herd, levels=c("Hart Ranges", "North Cariboo", "Barkerville",
                                              "Wells Grey North", "Columbia North", "Central Selkirks"))

p1.data <- p1.data[,c("Herd", "Observed_Extrap", "Shannon_Extrap")]
p1.data <- reshape2::melt(p1.data)

levels(p1.data$variable)[levels(p1.data$variable)=="Observed_Extrap"] <- "ASV richness"
levels(p1.data$variable)[levels(p1.data$variable)=="Shannon_Extrap"] <- "Shannon diversity"

p1 <- ggplot(p1.data, aes(x=Herd, y=value, fill=Herd)) + geom_boxplot() +
  facet_wrap(~variable, scales="free_y") +
  theme_bw() + plot_theme +
  labs(tag="a") +
  scale_fill_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)])) +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        legend.position="none")


# (b) Series of ordinations ####
p.data <- list()
x.data <- list()

for(i in c(1:5)){
  temp <- bdiv.ordinations.deepsnow$protists$points[[i]]
  
  temp$Herd <- as.character(temp$Herd)
  temp$Herd <- factor(temp$Herd, levels=c("Hart Ranges", "North Cariboo", "Barkerville",
                                          "Wells Grey North", "Columbia North", "Central Selkirks"))
  if(i==1){temp$distance <- "Bray-Curtis"}
  if(i==2){temp$distance <- "Jaccard"}
  if(i==3){temp$distance <- "Aitchison"}
  if(i==4){temp$distance <- "weighted UniFrac"}
  if(i==5){temp$distance <- "unweighted UniFrac"}
  
  p.data[[i]] <- temp
  
  temp <- bdiv.ordinations.deepsnow$protists$models[[i]]
  x.data[[i]] <- c(100*temp[["CA"]][["eig"]][[1]]/temp[["CA"]]$tot.chi, 
                   100*temp[["CA"]][["eig"]][[2]]/temp[["CA"]]$tot.chi)
  
  rm(temp)
}

p2 <- ggplot(p.data[[1]], aes(x=MDS1, y=MDS2, color=Herd)) + 
  geom_point(size=2, shape=16) +
  facet_wrap(~distance) +
  ggforce::geom_mark_ellipse(aes(color = Herd)) +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (15.9%)", y="PC2 (12.5%)", tag="b") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)]),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks"))

p3 <- ggplot(p.data[[2]], aes(x=MDS1, y=MDS2, color=Herd)) + 
  geom_point(size=2, shape=16) +
  facet_wrap(~distance) +
  ggforce::geom_mark_ellipse(aes(color = Herd)) +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (11.2%)", y="PC2 (9.7%)") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)]),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks"))

p4 <- ggplot(p.data[[3]], aes(x=PC1, y=PC2, color=Herd)) + 
  geom_point(size=2, shape=16) +
  facet_wrap(~distance) +
  ggforce::geom_mark_ellipse(aes(color = Herd)) +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (12.7%)", y="PC2 (8.5%)") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)]),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks"))

p5 <- ggplot(p.data[[4]], aes(x=MDS1, y=MDS2, color=Herd)) + 
  geom_point(size=2, shape=16) +
  facet_wrap(~distance) +
  ggforce::geom_mark_ellipse(aes(color = Herd)) +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (42.9%)", y="PC2 (31.6%)") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)]),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks"))

p6 <- ggplot(p.data[[5]], aes(x=MDS1, y=MDS2, color=Herd)) + 
  geom_point(size=2, shape=16) +
  facet_wrap(~distance) +
  ggforce::geom_mark_ellipse(aes(color = Herd)) +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=5)) +
  labs(x="PC1 (15.9%)", y="PC2 (11.7%)") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.position="none",
        legend.box="horizontal") +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)]),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks"))

p7 <- cowplot::get_legend(
  ggplot(p.data[[5]], aes(x=MDS1, y=MDS2, color=Herd)) + 
    geom_point(size=2, shape=16) +
    facet_wrap(~distance) +
    ggforce::geom_mark_ellipse(aes(color = Herd)) +
    theme_bw() + plot_theme +
    guides(color=guide_legend(ncol=1), line=FALSE) +
    labs(x="PC1 (32.8%)", y="PC2 (20.3%)") +
    theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
          axis.line = element_line(colour = "black", size=0.5),
          legend.title=element_blank(),
          legend.position="right") +
    scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)]),
                       labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                                "Wells Grey North", "Columbia North", "Central Selkirks"))
)

# (-) Assemble figure ####
grid.arrange(p1, p2, p3, p4, p5, p6, p7,
             layout_matrix=rbind(c(1,1,1),c(2,3,4),c(5,6,7)),
             heights=c(0.4,0.3,0.3))

p2 <- ggplotGrob(p2)
p3 <- ggplotGrob(p3)
p4 <- ggplotGrob(p4)
p5 <- ggplotGrob(p5)
p6 <- ggplotGrob(p6)

p2$widths <- p5$widths
p3$widths <- p5$widths
p4$widths <- p5$widths
p6$widths <- p5$widths

p3$heights <- p2$heights
p4$heights <- p2$heights
p5$heights <- p2$heights
p6$heights <- p2$heights

plot <- arrangeGrob(p1, p2, p3, p4, p5, p6, p7,
                    layout_matrix=rbind(c(1,1,1),c(2,3,4),c(5,6,7)),
                    heights=c(0.4,0.3,0.3))

ggsave("~/caribou/figures/raw_figS13_bacteria_deepsnow.jpg",
       plot, width=6.5, height=7.5, units="in", dpi=300)

rm(p1, p2, p3, p4, p5, p6, p7, plot)

#### *FIG. S14: PROTIST PANELS ####
# (a) Protist alpha diversity scatter plot ####
p1 <- ggplot() +
  geom_errorbar(data = alpha.means[[5]], aes(x=m.Observed_Extrap,
                                             ymin=m.Shannon_Extrap - sd.Shannon_Extrap,
                                             ymax=m.Shannon_Extrap + sd.Shannon_Extrap),
                color="grey", width=0) +
  geom_errorbarh(data = alpha.means[[5]], aes(y=m.Shannon_Extrap,
                                              xmin=m.Observed_Extrap - sd.Observed_Extrap,
                                              xmax=m.Observed_Extrap + sd.Observed_Extrap),
                 color="grey", height=0) +
  geom_point(data = sample_data.mb[[5]], aes(x=Observed_Extrap, y=Shannon_Extrap, color=Group), size=0.5) +
  
  geom_point(data = alpha.means[[5]], aes(x=m.Observed_Extrap, y=m.Shannon_Extrap,
                                          color=Group, shape=Group), size=3.5) +
  scale_shape_manual(values=c(16,17,15,18)) +
  theme_bw() + plot_theme +
  guides(shape=guide_legend(ncol=1, byrow=TRUE),
         color=guide_legend(ncol=1, byrow=TRUE, override.aes = list(size = 2, shape=c(16,17,15,18)))) +
  labs(x="ASV richness", y="Shannon diversity", tag="a") +
  scale_color_manual(values=c("chocolate1", "chartreuse3",
                              "darkorchid1", "gray62")) +
  theme(legend.position=c(0.73,0.18), legend.title=element_blank(),
        legend.key.height=unit(0.75,"line"),
        legend.key.width=unit(0.75,"line"),
        legend.text=element_text(size=7),
        legend.spacing.y=unit(0.002,'cm'))

# (b) Protist ordination ####
temp <- bdiv.ordinations.all[[5]][["points"]][[3]]
temp$Herd <- factor(temp$Herd, levels=c("Hart Ranges", "North Cariboo", "Barkerville",
                                        "Wells Grey North", "Columbia North", "Central Selkirks",
                                        "Tonquin", "Brazeau", "Revelstoke pen", "LARS"))

temp$GroupHerd <- paste0(temp$Group, temp$Herd)
temp$GroupHerd <- factor(temp$GroupHerd,
                         levels=c("deep-snowHart Ranges", "deep-snowNorth Cariboo",
                                  "deep-snowBarkerville", "deep-snowWells Grey North",
                                  "deep-snowColumbia North", "deep-snowCentral Selkirks",
                                  "shallow-snowTonquin", "shallow-snowBrazeau", "Revelstoke penRevelstoke pen",
                                  "LARSLARS"))

summary(bdiv.ordinations.all$protists$models[[3]])

p2 <- ggplot(data =  temp, aes(PC1, PC2)) + 
  geom_point(aes(color = GroupHerd, shape = GroupHerd), size=1.25) +
  stat_ellipse(aes(group=Group), level=0.95, linewidth=0.5, color="darkgrey") +
  theme_bw() + plot_theme +
  guides(color=guide_legend(ncol=1)) +
  labs(x="PC1 (23.5%)", y="PC2 (15.7%)", tag="b") +
  theme(panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.line = element_line(colour = "black", size=0.5),
        legend.title=element_blank(),
        legend.text=element_text(size=7),
        legend.position="right",
        legend.box="horizontal",
        legend.key.height=unit(0.75,"line"),
        legend.key.width=unit(0.75,"line"),
        legend.spacing.y=unit(0.01,"cm")) +
  scale_color_manual(name="Guide1", values=c(brewer.pal(8, "Set1")[c(1,2,4,5,7,8)], 
                                             brewer.pal(4, "Paired")[c(3:4)], 
                                             "darkorchid1", "gray62"),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS")) +
  scale_shape_manual(name="Guide1", values=c(16,16,16,16,16,16,17,17,15,18),
                     labels=c("Hart Ranges", "North Cariboo", "Barkerville",
                              "Wells Grey North", "Columbia North", "Central Selkirks",
                              "Tonquin", "Brazeau", "Revelstoke pen","LARS"))

# (c) Heat map ####
a.df.prt <- subset(diff.abund.GROUP[["protists"]][["Genus"]], glm.eBH < 0.05)
a.df.prt.hrd <- subset(diff.abund.HERD[["protists"]][["Genus"]], glm.eBH < 0.05)

df.prt <- subset(diff.abund.GROUP[["protists"]][["Genus"]], glm.eBH < 0.05 | Genus %in% a.df.prt.hrd$Genus)
df.prt.hrd <- subset(diff.abund.HERD[["protists"]][["Genus"]], glm.eBH < 0.05 | Genus %in% a.df.prt$Genus)

rm(a.df.prt, a.df.prt.hrd, df.prt.hrd)

df.prt2 <- subset(prev.data$protists$Genus, Taxon %in% df.prt$Genus)
df.prt2$Genus <- df.prt2$Taxon
df.prt2$Taxon <- NULL
df.prt2 <- subset(df.prt2, Abundance > 0.5 | `deep-snow` > 1)

df.prt2[,c("Prevalence", "Abundance")] <- NULL
rownames(df.prt2) <- df.prt2$Genus
df.prt2$Genus <- NULL
df.prt2 <- as.data.frame(t(df.prt2))
df.prt2$Group <- rownames(df.prt2)
df.prt2 <- reshape2::melt(df.prt2)

df.prt2$Facet <- rep("population")
df.prt2$Facet[df.prt2$Group %in% c("Barkerville", "Central Selkirks", "Hart Ranges",
                                   "North Cariboo", "Columbia North", "Wells Grey North")] <- "deep-snow"
df.prt2$Facet[df.prt2$Group %in% c("Brazeau", "Tonquin")] <- "shallow-snow"
df.prt2$Facet <- factor(df.prt2$Facet, levels=c(
  "population", "deep-snow", "shallow-snow"
))

df.prt2$Group <- factor(df.prt2$Group, levels=c(
  "deep-snow", "shallow-snow", "Revelstoke pen", "LARS",
  "Hart Ranges", "North Cariboo", "Barkerville", "Wells Grey North", "Columbia North", "Central Selkirks",
  "Tonquin", "Brazeau"
))

df.prt2$rounded <- as.numeric(df.prt2$value)
df.prt2$rounded[df.prt2$rounded==0] <- "-"
df.prt2$rounded[df.prt2$rounded != "-"] <- NA

df.prt2$mean_color <- as.numeric(as.character(df.prt2$rounded))

df.prt2$Genus <- df.prt2$variable
df.prt2$variable <- NULL
df.prt2 <- merge(df.prt2,
                 df.prt[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")],
                 by="Genus", all.x=TRUE, all.y=FALSE)

indic.prt <- indicator.Gini.synthesis[[5]]
indic.prt <- subset(indic.prt, taxon %in% df.prt2$Genus)

indic.prt2 <- subset(indic.prt, indicatee=="bothdeepsnow")
indic.prt2$indicatee <- rep("deep-snow")
indic.prt3 <- subset(indic.prt, indicatee=="bothdeepsnow")
indic.prt3$indicatee <- rep("Revelstoke pen")
indic.prt <- subset(indic.prt, indicatee !="bothdeepsnow" | is.na(indicatee=="TRUE"))

indic.prt4 <- subset(indic.prt, indicatee=="wild")
indic.prt4$indicatee <- rep("deep-snow")
indic.prt5 <- subset(indic.prt, indicatee=="wild")
indic.prt5$indicatee <- rep("shallow-snow")
indic.prt <- subset(indic.prt, indicatee !="wild" | is.na(indicatee=="TRUE"))

indic.prt <- rbind(indic.prt, indic.prt2, indic.prt3, indic.prt4, indic.prt5)
rm(indic.prt2, indic.prt3, indic.prt4, indic.prt5)

indic.prt$symbol <- rep("B")
indic.prt$symbol[indic.prt$p.value < 0.05] <- "A"
indic.prt$Group <- indic.prt$indicatee

df.prt2 <- merge(df.prt2, indic.prt[,c("Genus","A","B","symbol","Group")],
                 by=c("Group","Genus"), all=TRUE)
df.prt2 <- subset(df.prt2, is.na(Group)=="FALSE")
df.prt2$symbol[is.na(df.prt2$symbol)=="TRUE"] <- "B"

levels(df.prt2$Facet)[levels(df.prt2$Facet)=="shallow-snow"] <- "shallow-\nsnow"

df.prt2$Phylum[df.prt2$Phylum=="Archamoebae"] <- "Amoeb."

p3 <- ggplot(df.prt2, aes(x=Genus, y=Group, fill=sqrt(mean_color))) +
  geom_tile() +
  geom_text(aes(label=rounded), size=2.25) +
  geom_point(data=subset(df.prt2, symbol=="A"), 
             aes(x=Genus, y=Group)) +
  facet_grid(Facet~Phylum, scales="free", space="free") +
  scale_fill_gradient(low="white", high="firebrick3") +
  scale_color_manual(values=c("black", "darkgrey")) +
  theme_bw() + plot_theme + 
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0), limits=rev) +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1, size=7),
        strip.text = element_text(size=7),
        axis.text.y = element_text(size=7),
        legend.position="none") +
  labs(tag="c", x="Genus", y="Population / Herd\n")

legend <- cowplot::get_legend(
  ggplot(df.prt2, aes(x=Genus, y=Group, fill=sqrt(mean_color))) +
    geom_tile() +
    geom_text(aes(label=rounded), size=2.25) +
    geom_point(data=subset(df.prt2, symbol=="A"), 
               aes(x=Genus, y=Group)) +
    facet_grid(Facet~Phylum, scales="free", space="free") +
    scale_fill_gradient(low="white", high="firebrick3", 
                        breaks=c(0, 1, 2.236068, 3.162278, 4.1833, 5),
                        labels=c("0", "1", "5", "10", "17.5", "25"),
                        name="Relative\nabundance (%)") +
    scale_color_manual(values=c("black", "darkgrey")) +
    theme_bw() + plot_theme + 
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0), limits=rev) +
    theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1, size=7),
          strip.text = element_text(size=7),
          axis.text.y = element_text(size=7),
          legend.title=element_text(hjust=0.5, size=7),
          legend.position="right") +
    labs(tag="c", x="Genus", y="Population / Herd\n")
)


# (-) Assemble ####
grid.arrange(p1, p2, p3, layout_matrix=rbind(c(1,1,1,2,2,2,2),c(3,3,3,3,3,3,3)), heights=c(0.4,0.6))

p1 <- ggplotGrob(p1)
p2 <- ggplotGrob(p2)
p3 <- ggplotGrob(p3)

plot <- arrangeGrob(p1, p2, p3, layout_matrix=rbind(c(1,1,1,2,2,2,2),c(3,3,3,3,3,3,3)), heights=c(0.45,0.55))

ggsave("~/caribou/figures/raw_figS14_protist_panels.jpg", plot, width=6.5, height=7, units="in", dpi=300)
ggsave("~/caribou/figures/raw_figS14_legend.jpg", legend, width=6.5, height=4, units="in", dpi=300)

rm(p1, p2, p3, plot, df.prt, df.prt2, indic.prt)

#### *FIG. S15: OTHER FUNGAL READS ####
a.df.FITS <- subset(diff.abund.GROUP[["FITS"]][["Genus"]], glm.eBH < 0.05)
a.df.FITS.hrd <- subset(diff.abund.HERD[["FITS"]][["Genus"]], glm.eBH < 0.05)

df.FITS <- subset(diff.abund.GROUP[["FITS"]][["Genus"]], glm.eBH < 0.05 | Genus %in% a.df.FITS.hrd$Genus)
df.FITS.hrd <- subset(diff.abund.HERD[["FITS"]][["Genus"]], glm.eBH < 0.05 | Genus %in% a.df.FITS$Genus)

rm(a.df.FITS, a.df.FITS.hrd, df.FITS.hrd)

df.FITS2 <- subset(prev.data$FITS$Genus, Taxon %in% df.FITS$Genus)
df.FITS2$Genus <- df.FITS2$Taxon
df.FITS2$Taxon <- NULL
df.FITS2 <- subset(df.FITS2, Abundance > 0.5 | `deep-snow` > 1)

df.FITS2[,c("Prevalence", "Abundance")] <- NULL
rownames(df.FITS2) <- df.FITS2$Genus
df.FITS2$Genus <- NULL
df.FITS2 <- as.data.frame(t(df.FITS2))
df.FITS2$Group <- rownames(df.FITS2)
df.FITS2 <- reshape2::melt(df.FITS2)

df.FITS2$Facet <- rep("group")
df.FITS2$Facet[df.FITS2$Group %in% c("Barkerville", "Central Selkirks", "Hart Ranges",
                                   "North Cariboo", "Columbia North", "Wells Grey North")] <- "deep-snow"
df.FITS2$Facet[df.FITS2$Group %in% c("Brazeau", "Tonquin")] <- "shallow-snow"

df.FITS2$Group <- factor(df.FITS2$Group, levels=c(
  "deep-snow", "shallow-snow", "Revelstoke pen", "LARS",
  "Hart Ranges", "North Cariboo", "Barkerville", "Wells Grey North", "Columbia North", "Central Selkirks",
  "Tonquin", "Brazeau"
))
df.FITS2$Facet <- factor(df.FITS2$Facet, levels=c(
  "population", "deep-snow", "shallow-snow"
))

df.FITS2$rounded <- as.numeric(df.FITS2$value)
df.FITS2$rounded[df.FITS2$rounded==0] <- "-"
df.FITS2$rounded[df.FITS2$rounded != "-"] <- NA

df.FITS2$mean_color <- as.numeric(as.character(df.FITS2$rounded))
df.FITS2$mean_color[df.FITS2$mean_color > 50] <- 50

df.FITS2$Genus <- df.FITS2$variable
df.FITS2$variable <- NULL
df.FITS2 <- merge(df.FITS2,
                 df.FITS[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")],
                 by="Genus", all.x=TRUE, all.y=FALSE)

indic.FITS <- indicator.Gini.synthesis[[3]]
indic.FITS <- subset(indic.FITS, taxon %in% df.FITS2$Genus)

indic.FITS2 <- subset(indic.FITS, indicatee=="bothdeepsnow")
indic.FITS2$indicatee <- rep("deep-snow")
indic.FITS3 <- subset(indic.FITS, indicatee=="bothdeepsnow")
indic.FITS3$indicatee <- rep("Revelstoke pen")
indic.FITS <- subset(indic.FITS, indicatee !="bothdeepsnow" | is.na(indicatee=="TRUE"))

indic.FITS4 <- subset(indic.FITS, indicatee=="wild")
indic.FITS4$indicatee <- rep("deep-snow")
indic.FITS5 <- subset(indic.FITS, indicatee=="wild")
indic.FITS5$indicatee <- rep("shallow-snow")
indic.FITS <- subset(indic.FITS, indicatee !="wild" | is.na(indicatee=="TRUE"))

indic.FITS <- rbind(indic.FITS, indic.FITS2, indic.FITS3, indic.FITS4, indic.FITS5)
rm(indic.FITS2, indic.FITS3, indic.FITS4, indic.FITS5)

indic.FITS$symbol <- rep("B")
indic.FITS$symbol[indic.FITS$p.value < 0.05] <- "A"
indic.FITS$Group <- indic.FITS$indicatee

indic.FITS$Genus <- indic.FITS$taxon

df.FITS2 <- merge(df.FITS2, indic.FITS[,c("Genus","A","B","symbol","Group")],
                 by=c("Group","Genus"), all=TRUE)
df.FITS2 <- subset(df.FITS2, is.na(Group)=="FALSE")
df.FITS2$symbol[is.na(df.FITS2$symbol)=="TRUE"] <- "B"

levels(df.FITS2$Facet)[levels(df.FITS2$Facet)=="shallow-snow"] <- "shallow-\nsnow"

df.FITS2 <- subset(df.FITS2, Genus !="Uncl_Fungi")
df.FITS2 <- subset(df.FITS2, Genus !="Uncl_Ascomycota")

df.FITS2$Phylum[df.FITS2$Phylum=="Ascomycota"] <- df.FITS2$Class[df.FITS2$Phylum=="Ascomycota"]

levels(factor(df.FITS2$Phylum))
df.FITS2$Phylum[df.FITS2$Phylum=="Basidiomycota"] <- "Basidio-\nmycota"
df.FITS2$Phylum[df.FITS2$Phylum=="Eurotiomycetes"] <- "Eurotio-\nmycetes"
df.FITS2$Phylum[df.FITS2$Phylum=="Lecanoromycetes"] <- "Lecanoro-\nmycetes"
df.FITS2$Phylum[df.FITS2$Phylum=="Sordariomycetes"] <- "Sordario-\nmycetes"

df.FITS2$Phylum <- factor(df.FITS2$Phylum, levels=c("Leotiomycetes", "Dothideomycetes", "Eurotio-\nmycetes",
                                                    "Lecanoro-\nmycetes", "Sordario-\nmycetes", "Basidio-\nmycota"))

plot <- ggplot(df.FITS2, aes(x=Group, y=Genus, fill=sqrt(mean_color))) +
  geom_tile() +
  geom_text(aes(label=rounded), size=2.25) +
  geom_point(data=subset(df.FITS2, symbol=="A"), 
             aes(x=Group, y=Genus)) +
  facet_grid(Phylum~Facet, scales="free", space="free") +
  scale_fill_gradient(low="white", high="firebrick3", 
                      breaks=c(0, 2.236068, 3.162278, 4.1833, 5, 7.061068),
                      labels=c("0", "5", "10", "17.5", "25", ">50"),
                      name="Relative\nabundance (%)") +
  scale_color_manual(values=c("black", "darkgrey")) +
  theme_bw() + plot_theme + 
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0), limits=rev) +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1, size=7),
        strip.text = element_text(size=7),
        axis.text.y = element_text(size=7),
        legend.title=element_text(hjust=0.5, size=7),
        legend.position="right") +
  labs(x="Population / Herd", y="Genus")

ggsave("~/caribou/figures/raw_figS16_other_fungi_hetmap.jpg", plot, width=6.5, height=7, units="in", dpi=300)

rm(plot, df.FITS2, df.FITS, df.FITS.hrd, indic.FITS)

#### FIG. S16: GENOME HEAT MAP ####
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

plot_tab <- subset(plot_tab, facet=="GH")

plot_tab$variable <- as.character(plot_tab$variable)

plot_tab$variable[plot_tab$variable=="med.massiliensis"] <- "Mediterranea massiliensis"

plot_tab$variable[plot_tab$variable=="p.dorei"] <- "Phocaeicola dorei"
plot_tab$variable[plot_tab$variable=="p.vulgatus"] <- "Phocaeicola vulgatus"

plot_tab$variable[plot_tab$variable=="d.dubosii"] <- "Duncaniella dubosii"
plot_tab$variable[plot_tab$variable=="d.freteri"] <- "Duncaniella freteri"
plot_tab$variable[plot_tab$variable=="d.muris"] <- "Duncaniella muris"
plot_tab$variable[plot_tab$variable=="d.muricolitica"] <- "Duncaniella muricolitica"

plot_tab$variable[plot_tab$variable=="p.intestinale"] <- "Paramuribaculum intestinale"

plot_tab$variable[plot_tab$variable=="m.gordoncarteri"] <- "Muribaculum gordoncarteri"
plot_tab$variable[plot_tab$variable=="m.caecicola"] <- "Muribaculum caecicola"
plot_tab$variable[plot_tab$variable=="m.intestinale"] <- "Muribaculum intestinale"

plot_tab$variable <- factor(plot_tab$variable, levels=c("Mediterranea massiliensis",
                                                        "Phocaeicola dorei",
                                                        "Phocaeicola vulgatus",
                                                        "Duncaniella dubosii",
                                                        "Duncaniella freteri",
                                                        "Duncaniella muris",
                                                        "Duncaniella muricolitica",
                                                        "Muribaculum gordoncarteri",
                                                        "Muribaculum caecicola",
                                                        "Muribaculum intestinale",
                                                        "Paramuribaculum intestinale"))

plot <- ggplot(subset(plot_tab, facet=="GH"), aes(x=cluster, y=variable, fill=sqrt(value))) + 
  geom_tile() + theme_bw() + plot_theme +
  scale_fill_gradient(low="white", high="firebrick3", 
                      breaks=c(0, 1, 2, 3, 4, 5),
                      labels=c("0", "1", "4", "9", "16", "25"),
                      name="Number of\nhits") +
  #facet_grid(~facet, scales="free_x", space="free_x") +
  scale_y_discrete(limits=rev) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=6, vjust=0.5),
        strip.text = element_text(size=7),
        axis.text.y = element_text(size=7),
        legend.title=element_text(hjust=0.5, size=7)) +
  labs(x="Glycoside hydrolase (GH) cluster", y="Species")

ggsave("~/FigS16_genome_map.jpg", plot, width=8, height=3.5, units="in", dpi=300)

#### *TABLE S3: LICHENS AT LARS ####
lars <- subset_samples(pseq.clr.sub$FITS, Group=="LARS")
lars <- prune_taxa(taxa_sums(lars) > 0, lars)

lars <- transform_sample_counts(lars, function(x) 100 * x/sum(x))
lars.lich <- subset_taxa(lars, taxa_names(lars) %in% taxa_names(pseq.clr.sub$lichens))

lars.lich.df <- as.data.frame(cbind(otu_table(lars.lich)))
lars.lich.df <- as.data.frame(t(lars.lich.df))

lars.lich.df$prev <- rowSums(lars.lich.df > 0)
lars.lich.df$mean <- rowMeans(lars.lich.df[,c(1:13)])

lars.lich.df <- cbind(as.data.frame(cbind(tax_table(lars.lich))),
                      lars.lich.df)

#### *TABLE S5: RANDOM FOREST CONFUSION MATRIX ####
p.data <- list()
for(i in c(1:7)){
  temp <- as.data.frame(rf.results.all[[i]]$Models$confusion)
  temp$amplicon <- rep(names(rf.results.all)[i])
  p.data[[i]] <- temp
  rm(temp)
}

p.data <- dplyr::bind_rows(p.data)

write.csv(p.data, "~/rf.confusion.csv")

#### *TABLE S6: LICHENS DIFF ABUND / GINI / INDICATOR ####
diff.abund <- diff.abund.GROUP$lichens$Genus %>%
  dplyr::select(., -kw.ep, -kw.eBH, -glm.ep)

temp.indic <- indicator.Gini.synthesis[[7]]

table <- merge(diff.abund, temp.indic, by="Genus", all=TRUE)
table$Taxon <- table$Genus
table <- merge(table, prev.data$lichens$Genus[,c("Taxon", "Prevalence")], by="Taxon", all=TRUE)

write.csv(table, "~/lichen.megatable.csv")

#### *TABLE S7: ALGAE DIFF ABUND / GINI / INDICATOR ####
diff.abund <- diff.abund.GROUP$algae$Genus %>%
  dplyr::select(., -kw.ep, -kw.eBH, -glm.ep)

temp.indic <- indicator.Gini.synthesis[[6]]

table <- merge(diff.abund, temp.indic, by="Genus", all=TRUE)
table$Taxon <- table$Genus
table <- merge(table, prev.data$algae$Genus[,c("Taxon", "Prevalence")], by="Taxon", all=TRUE)

write.csv(table, "~/algae.megatable.csv")

#### *TABLE S8: PLANTS DIFF ABUND / GINI / INDICATOR ####
diff.abund <- diff.abund.GROUP$PITS$Genus %>%
  dplyr::select(., -kw.ep, -kw.eBH, -glm.ep)

temp.indic <- indicator.Gini.synthesis[[4]]

table <- merge(diff.abund, temp.indic, by="Genus", all=TRUE)
table$Taxon <- table$Genus
table <- merge(table, prev.data$PITS$Genus[,c("Taxon", "Prevalence")], by="Taxon", all=TRUE)

write.csv(table, "~/PITS.megatable.csv")

#### *TABLE S9: 16S DIFF ABUND / GINI / INDICATOR ####
diff.abund <- diff.abund.GROUP[["X16S"]]$Genus %>%
  dplyr::select(., -kw.ep, -kw.eBH, -glm.ep)

temp.indic <- indicator.Gini.synthesis[[1]]

table <- merge(diff.abund, temp.indic, by="Genus", all=TRUE)
table$Taxon <- table$Genus
table <- merge(table, prev.data$X16S$Genus[,c("Taxon", "Prevalence")], by="Taxon", all=TRUE)

write.csv(table, "~/bact.megatable.csv")

#### *TABLE S10: PROTISTS DIFF ABUND / GINI / INDICATOR ####
diff.abund <- diff.abund.GROUP$protists$Genus %>%
  dplyr::select(., -kw.ep, -kw.eBH, -glm.ep)

temp.indic <- indicator.Gini.synthesis[[5]]

table <- merge(diff.abund, temp.indic, by="Genus", all=TRUE)
table$Taxon <- table$Genus
table <- merge(table, prev.data$protists$Genus[,c("Taxon", "Prevalence")], by="Taxon", all=TRUE)

write.csv(table, "~/protist.megatable.csv")

#### *TABLE S11: FUNGI DIFF ABUND / GINI / INDICATOR ####
diff.abund <- diff.abund.GROUP$FITS$Genus %>%
  dplyr::select(., -kw.ep, -kw.eBH, -glm.ep)

temp.indic <- indicator.Gini.synthesis[[3]]
temp.indic$Genus <- rownames(temp.indic)

table <- merge(diff.abund, temp.indic, by="Genus", all=TRUE)
table$Taxon <- table$Genus
table <- merge(table, prev.data$FITS$Genus[,c("Taxon", "Prevalence")], by="Taxon", all=TRUE)

write.csv(table, "~/fungal.megatable.csv")
